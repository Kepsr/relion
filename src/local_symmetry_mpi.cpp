#include "src/local_symmetry_mpi.h"
//#define DEBUG

void local_symmetry_parameters_mpi::read(int argc, char **argv) {
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    local_symmetry_parameters::read(argc, argv);

    // Don't put any output to screen for mpi followers
    verb = node->isLeader();

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
    printMpiNodesMachineNames(*node);
}

void local_symmetry_parameters_mpi::run() {
    RFLOAT aa = 0, bb = 0, gg = 0.0, tmp_binning_factor = 1.0;
    RFLOAT dx = 0.0, dy = 0.0, dz = 0.0, cc = 0.0;
    RFLOAT mask_sum = 0.0, mask_ctr = 0.0, mask2_sum = 0.0, mask2_ctr = 0.0;

    Image<RFLOAT> map, mask, mask2;
    map.clear(); mask.clear(); mask2.clear();
    Matrix1D<RFLOAT> op_search_ranges, op, com0_int, com1_int, com1_float, com1_diff, vecR3;
    op_search_ranges.clear(); op.clear(); com0_int.clear(); com1_int.clear(); com1_float.clear(); com1_diff.clear(); vecR3.clear();
    std::vector<std::vector<Matrix1D<RFLOAT> > > op_list;
    op_list.clear();
    std::vector<Matrix1D<RFLOAT> > op_samplings, op_samplings_batch;
    op_samplings.clear(); op_samplings_batch.clear();
    MultidimArray<RFLOAT> op_samplings_batch_packed, src_cropped, dest_cropped, mask_cropped;
    op_samplings_batch_packed.clear(); src_cropped.clear(); dest_cropped.clear(); mask_cropped.clear();

    // Check options
    if (
        do_apply_local_symmetry || do_duplicate_local_symmetry || 
        do_txt2rln || do_transform || do_debug || 
        !do_local_search_local_symmetry_ops
    ) REPORT_ERROR("ERROR: Please specify '--search' as the only option! For other options use non-parallel version (without '_mpi') instead!");

    // Leader writes out commands
    if (!show_usage_for_an_option && !do_debug && node->isLeader()) {
        local_symmetry_parameters::writeCommand("relion_localsym.log", "mpirun -n " + integerToString(node->size) + " `which relion_localsym_mpi`");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<FileName> fn_mask_list; fn_mask_list.clear();
    std::vector<std::vector<FileName> > op_mask_list; op_mask_list.clear();
    if (node->isLeader()) {
        displayEmptyLine();

        #ifdef DEBUG
        std::cout << " DEBUG: relion_localsym_mpi is running ..." << std::endl;
        #endif

        // Leader gets search ranges (in degrees and pixels), sets offset_step (in pixels).
        if (angpix_image < 0.001)
            REPORT_ERROR("Invalid pixel size!");
        if (fn_op_mask_info_in != "None") {
            if (
                ang_range      < Xmipp::epsilon && 
                ang_rot_range  < Xmipp::epsilon && 
                ang_tilt_range < Xmipp::epsilon && 
                ang_psi_range  < Xmipp::epsilon
            ) {
                ang_range = 180.0;
                std::cout << " Initial searches: reset searching ranges of all 3 Euler angles to +/-180 degrees." << std::endl;
            } else {
                std::cout << " User-defined initial searches: ";
                if (ang_range > Xmipp::epsilon) {
                    std::cout << "searching ranges of all 3 Euler angles are set to +/-" << ang_range << " degree(s).";
                } else {
                    std::cout << "(rot, tilt, psi) ranges are +/- (" << ang_rot_range << ", " << ang_tilt_range << ", " << ang_psi_range << ") degree(s).";
                }
                std::cout << std::endl;
            }
        }
        Localsym_composeOperator(
            op_search_ranges,
            ang_range    > Xmipp::epsilon ? ang_range    : ang_rot_range,
            ang_range    > Xmipp::epsilon ? ang_range    : ang_tilt_range,
            ang_range    > Xmipp::epsilon ? ang_range    : ang_psi_range,
            offset_range > Xmipp::epsilon ? offset_range : offset_x_range,
            offset_range > Xmipp::epsilon ? offset_range : offset_y_range,
            offset_range > Xmipp::epsilon ? offset_range : offset_z_range
        );
        Localsym_scaleTranslations(op_search_ranges, 1.0 / angpix_image);
        offset_step /= angpix_image;

        // Leader parses and reads mask info file
        // Local searches
        if (fn_op_mask_info_in == "None") {
            if (fn_info_in.getExtension() == "star") {
                readRelionFormatMasksAndOperators(
                    fn_info_in, fn_mask_list, op_list, angpix_image, true
                );
            } else {
                FileName fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
                parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
                readDMFormatMasksAndOperators(
                    fn_parsed, fn_mask_list, op_list, angpix_image, true
                );
            }
        } else {
            // Global searches
            std::cout << " Global searches: option --i_mask_info " << fn_info_in << " is ignored." << std::endl;
            readRelionFormatMasksWithoutOperators(
                fn_op_mask_info_in, fn_mask_list, 
                op_list, op_mask_list, ang_range > 179.99, true
            );
        }

        // Leader set total number of masks
        int nr_masks = fn_mask_list.size();

        // Leader reads input map
        std::cout << std::endl << " Pixel size = " << angpix_image << " Angstrom(s)" << std::endl;
        std::cout << " Read input map " << fn_unsym << " ..." << std::endl;
        map.read(fn_unsym);
        map().setXmippOrigin();
        if (!isMultidimArray3DCubic(map()))
            REPORT_ERROR("ERROR: Input map " + fn_unsym + " is not 3D cube!");
        #ifdef DEBUG
        std::cout << " I am leader. The nxyzdim of map() is " << map().size() << std::endl;
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Leader broadcasts total number of masks
    node->relion_MPI_Bcast(&nr_masks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // All nodes loop over all masks
    for (int imask = 0; imask < nr_masks; imask++) {
        MPI_Barrier(MPI_COMM_WORLD);

        long int newdim = 0, cropdim = 0;
        if (node->isLeader()) {
            displayEmptyLine();

            // Leader reads and checks the mask
            std::cout << " Read mask #" << imask + 1 << ": " << fn_mask_list[imask] << " ..." << std::endl;
            mask.read(fn_mask_list[imask]);
            mask().setXmippOrigin();
            if (!isMultidimArray3DCubic(mask()))
                REPORT_ERROR("ERROR: Input mask " + fn_mask_list[imask] + " is not 3D cube!");
            if (!map().sameShape(mask()))
                REPORT_ERROR("ERROR: Input map " + fn_unsym + " and mask " + fn_mask_list[imask] + " should have the same size!");
            sum3DCubicMask(mask(), mask_sum, mask_ctr);

            // Get com0 of this mask. Assume that com0 has all integer values!
            getMinCropSize(mask(), com0_int, cropdim, offset_range / angpix_image);
            if (cropdim < 2)
                REPORT_ERROR("ERROR: Mask " + fn_mask_list[imask] + " is too small!");
            XX(com0_int) = round(XX(com0_int));
            YY(com0_int) = round(YY(com0_int));
            ZZ(com0_int) = round(ZZ(com0_int));
            std::cout << " Mask #" << imask + 1 << " : center of mass XYZ = (" << XX(com0_int) << ", " << YY(com0_int) << ", " << ZZ(com0_int) << ") pixel(s)."<< std::endl;

            // Crop the mask and the corresponding region of the map
            long int z0 = round(ZZ(com0_int)) + Xmipp::init(cropdim);
            long int zf = round(ZZ(com0_int)) + Xmipp::last(cropdim);
            long int y0 = round(YY(com0_int)) + Xmipp::init(cropdim);
            long int yf = round(YY(com0_int)) + Xmipp::last(cropdim);
            long int x0 = round(XX(com0_int)) + Xmipp::init(cropdim);
            long int xf = round(XX(com0_int)) + Xmipp::last(cropdim);

            std::cout << " Mask #" << imask + 1 << " : cropped box size = " << cropdim << " pixels." << std::endl;
            #ifdef DEBUG
            std::cout << " Window: x0, y0, z0 = " << x0 << ", " << y0 << ", " << z0 << "; xf, yf, zf = " << xf << ", " << yf << ", " << zf << std::endl;
            #endif
            mask_cropped = mask().windowed(x0, xf, y0, yf, z0, zf).setXmippOrigin();
            src_cropped  = map ().windowed(x0, xf, y0, yf, z0, zf).setXmippOrigin();

            // Rescale the map and the mask (if binning_factor > 1), set 'newdim'.
            tmp_binning_factor = 1.0;
            newdim = cropdim;
            if (binning_factor - 1.0 > Xmipp::epsilon) {
                newdim = ceil(RFLOAT(cropdim) / binning_factor);
                if (newdim < 2)
                    REPORT_ERROR("ERROR: Binning factor is too large / Mask is too small!");
                if (newdim + 1 < cropdim) {
                    // Need rescaling
                    // Dimension should always be even
                    newdim += newdim % 2;
                    resizeMap(mask_cropped, newdim);
                    mask_cropped.setXmippOrigin();
                    resizeMap(src_cropped, newdim);
                    src_cropped.setXmippOrigin();
                    tmp_binning_factor = RFLOAT(cropdim) / RFLOAT(newdim);
                    std::cout << " + Rescale cropped box size from " << cropdim << " to " << newdim << " pixels. Binning factor = " << tmp_binning_factor << std::endl;

                    // Mask values might go out of range after rescaling. Fix it if it happens
                    truncateMultidimArray(mask_cropped, 0.0, 1.0);
                } else {
                    newdim = cropdim;
                }
            }
            #ifdef DEBUG
            std::cout << " newdim= " << newdim << ", cropdim= " << cropdim << std::endl;
            #endif
        }
        MPI_Barrier(MPI_COMM_WORLD);

        node->relion_MPI_Bcast(&newdim, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // Follower allocate space for MultidimArray
        // TODO: check whether the space is allocated and the map is read and successfully broadcast!!!
        if (!node->isLeader()) {
            mask_cropped.clear();
            mask_cropped.initZeros(1, newdim, newdim, newdim);
            mask_cropped.setXmippOrigin();
            src_cropped.clear();
            src_cropped.initZeros(1, newdim, newdim, newdim);
            src_cropped.setXmippOrigin();
            dest_cropped.clear();
            dest_cropped.initZeros(1, newdim, newdim, newdim);
            dest_cropped.setXmippOrigin();
            #ifdef DEBUG
            if (node->rank == 1) {
                std::cout << " I am node rank 1. The nxyzdim of mask_cropped is " << mask_cropped.size() << std::endl;
                std::cout << " I am node rank 1. The nxyzdim of  src_cropped is " <<  src_cropped.size() << std::endl;
                std::cout << " I am node rank 1. The nxyzdim of dest_cropped is " << dest_cropped.size() << std::endl;
            }
            #endif
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Leader broadcasts the mask to all nodes
        #ifdef DEBUG
        if (node->isLeader())
            std::cout << " Leader is broadcasting cropped masked region #" << (imask + 1) << " ..." << std::endl;
        #endif
        node->relion_MPI_Bcast(mask_cropped.data, mask_cropped.size(), relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(src_cropped.data,  src_cropped.size(),  relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
        #ifdef DEBUG
        if (node->isLeader())
            std::cout << " Leader has completed broadcasting cropped masked region #" << (imask + 1) << "." << std::endl;
        #endif

        // Leader reads total number of operators for this mask
        int nr_ops = 0;
        if (node->isLeader()) {
            nr_ops = op_list[imask].size();
            #ifdef DEBUG
            std::cout << " nr_ops= " << nr_ops << std::endl;
            #endif
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Leader broadcasts total number of operators for this mask to all followers
        node->relion_MPI_Bcast(&nr_ops, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // All nodes loop over all operators of this mask
        for (int iop = 0; iop < nr_ops; iop++) {
            MPI_Barrier(MPI_COMM_WORLD);

            // Leader gets sampling points
            int nr_total_samplings = 0;
            if (node->isLeader()) {
                std::cout << std::endl;

                com1_float.resize(3);
                std::fill(com1_float.begin(), com1_float.end(), 0);
                com1_int.resize(3);
                std::fill(com1_int.begin(), com1_int.end(), 0);
                com1_diff.resize(3);
                std::fill(com1_diff.begin(), com1_diff.end(), 0);

                Localsym_decomposeOperator(op_list[imask][iop], aa, bb, gg, dx, dy, dz, cc);

                if (fn_op_mask_info_in == "None") {
                    // Local searches
                    // Get com1_float. (floating point numbers)
                    // Com1f = R * Com0 + v
                    Matrix2D<RFLOAT> mat1 = Euler::angles2matrix(aa, bb, gg);
                    com1_float = matmul(mat1, com0_int) + vectorR3(dx, dy, dz);
                } else {
                    // Global searches
                    // Leader reads and checks the mask
                    std::cout << " Read mask #" << imask + 1 << " operator #" << iop + 1 << " : " << op_mask_list[imask][iop] << " ..." << std::endl;
                    mask2.read(op_mask_list[imask][iop]);
                    mask2().setXmippOrigin();
                    if (!isMultidimArray3DCubic(mask2()))
                        REPORT_ERROR("ERROR: Input mask " + op_mask_list[imask][iop] + " is not 3D cube!");
                    if (!map().sameShape(mask2()))
                        REPORT_ERROR("ERROR: Input map " + fn_unsym + " and mask " + op_mask_list[imask][iop] + " should have the same size!");
                    sum3DCubicMask(mask2(), mask2_sum, mask2_ctr);
                    if (!similar3DCubicMasks(mask_sum, mask_ctr, mask2_sum, mask2_ctr))
                        std::cerr << " WARNING: masks " << fn_mask_list[imask] << " and " << op_mask_list[imask][iop] << " seem different! Please check whether they are covering regions from the same set!" << std::endl;

                    // Calculate Com1f of this mask
                    mask2().centerOfMass(com1_float);
                    std::cout << " Mask #" << imask + 1 << " operator #" << iop + 1 << " : center of mass XYZ = (" << XX(com1_float) << ", " << YY(com1_float) << ", " << ZZ(com1_float) << ") pixel(s)."<< std::endl;
                }

                // Get com1_int and com1_diff
                // diff = Com1f - Com1i
                XX(com1_int) = round(XX(com1_float));
                YY(com1_int) = round(YY(com1_float));
                ZZ(com1_int) = round(ZZ(com1_float));
                XX(com1_diff) = XX(com1_float) - XX(com1_int);
                YY(com1_diff) = YY(com1_float) - YY(com1_int);
                ZZ(com1_diff) = ZZ(com1_float) - ZZ(com1_int);

                // Crop this region
                long int z0 = round(ZZ(com1_int)) + Xmipp::init(cropdim);
                long int zf = round(ZZ(com1_int)) + Xmipp::last(cropdim);
                long int y0 = round(YY(com1_int)) + Xmipp::init(cropdim);
                long int yf = round(YY(com1_int)) + Xmipp::last(cropdim);
                long int x0 = round(XX(com1_int)) + Xmipp::init(cropdim);
                long int xf = round(XX(com1_int)) + Xmipp::last(cropdim);
                #ifdef DEBUG
                std::cout << " Window: x0, y0, z0 = " << x0 << ", " << y0 << ", " << z0 << "; xf, yf, zf = " << xf << ", " << yf << ", " << zf << std::endl;
                #endif
                dest_cropped = map().windowed(x0, xf, y0, yf, z0, zf).setXmippOrigin();

                // Do the same rescaling
                if (newdim != cropdim) {
                    resizeMap(dest_cropped, newdim);
                    dest_cropped.setXmippOrigin();
                }

                // Leader gets sampling points
                // Get sampling points - Rescale translational search ranges and steps
                Localsym_composeOperator(op, aa, bb, gg, XX(com1_diff), YY(com1_diff), ZZ(com1_diff), cc);
                if (newdim != cropdim) {
                    Localsym_scaleTranslations(op_search_ranges, 1.0 / tmp_binning_factor);
                    offset_step *= 1.0 / tmp_binning_factor;
                    Localsym_scaleTranslations(op, 1.0 / tmp_binning_factor);
                }
                #ifdef __unix__
                std::cout << " + Refining " << "\e[1m" << "Mask #" << imask + 1 << " Operator #" << iop + 1 << "\e[0m" << ": " << std::flush;
                #else
                std::cout << " + Refining Mask #" << imask + 1 << " Operator #" << iop + 1 << ": " << std::flush;
                #endif
                Localsym_outputOperator(op_list[imask][iop], &std::cout, angpix_image);
                std::cout << std::endl;
                getLocalSearchOperatorSamplings(
                    op, op_search_ranges, op_samplings,
                    ang_step, offset_step,
                    use_healpix_sampling, true
                );
                if (newdim != cropdim) {
                    Localsym_scaleTranslations(op_search_ranges, tmp_binning_factor);
                    offset_step *= tmp_binning_factor;
                    Localsym_scaleTranslations(op, tmp_binning_factor);
                }

                if (op_samplings.size() <= node->size)
                    REPORT_ERROR("ERROR: Too few sampling points! Use non-parallel version (without '_mpi') instead!");

                nr_total_samplings = op_samplings.size();
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Leader sends this 'dest' cropped region to all followers
            node->relion_MPI_Bcast(dest_cropped.data, dest_cropped.size(), relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);

            // Leader sends the number of total samplings to all followers
            node->relion_MPI_Bcast(&nr_total_samplings, 1, MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

            // All nodes allocate space for op_samplings_batch_packed
            // Allow 10 more empty units to prevent segmentation fault
            op_samplings_batch_packed.initZeros(nr_total_samplings / node->size + 10, NR_LOCALSYM_PARAMETERS);

            MPI_Barrier(MPI_COMM_WORLD);

            // Leader distributes sampling points to all followers
            long int first = 0, last = 0;
            if (node->isLeader()) {
                for (int id_rank = (node->size) - 1; id_rank >= 0; id_rank--) {
                    divide_equally(nr_total_samplings, node->size, id_rank, first, last);

                    // Beware: Ysize(op_samplings_batch_packed) is larger than (last - first + 1)
                    for (long int j = 0; j < Ysize(op_samplings_batch_packed); j++)
                    for (long int i = 0; i < Xsize(op_samplings_batch_packed); i++) {
                        if (i >= 0 && i <= last - first)
                            direct::elem(op_samplings_batch_packed, i, j) = op_samplings[i + first][j];
                    }

                    // Leader distributes sampling points to all followers
                    if (id_rank > 0)
                        node->relion_MPI_Send(op_samplings_batch_packed.data, (last - first + 1) * NR_LOCALSYM_PARAMETERS, relion_MPI::DOUBLE, id_rank, MPITag::LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD);
                    // If id_rank == 0 (leader), just keep op_samplings_batch_packed to the leader itself
                }
            } else {
                MPI_Status status;
                // Followers receive sampling points from leader
                // Important: Followers calculate first and last subscripts!
                divide_equally(nr_total_samplings, node->size, node->rank, first, last);
                node->relion_MPI_Recv(op_samplings_batch_packed.data, (last - first + 1) * NR_LOCALSYM_PARAMETERS, relion_MPI::DOUBLE, 0, MPITag::LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD, status);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // All nodes unpack sampling points
            op_samplings_batch.clear();
            for (long int j = 0; j < Ysize(op_samplings_batch_packed); j++) {
                op.resize(NR_LOCALSYM_PARAMETERS);
                std::fill(op.begin(), op.end(), 0);
                for (long int i = 0; i < Xsize(op_samplings_batch_packed); i++) {
                    op[i] = direct::elem(op_samplings_batch_packed, i, j);
                }
                op_samplings_batch.push_back(op);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // All nodes calculate CC, with leader profiling (DONT SORT!)
            calculateOperatorCC(src_cropped, dest_cropped, mask_cropped, op_samplings_batch, false, node->isLeader());
            for (int op_id = 0; op_id < op_samplings_batch.size(); op_id++) {
                direct::elem(op_samplings_batch_packed, op_id, CC_POS) = op_samplings_batch[op_id][CC_POS];
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Followers send their results back to leader
            if (!node->isLeader()) {
                node->relion_MPI_Send(op_samplings_batch_packed.data, (last - first + 1) * NR_LOCALSYM_PARAMETERS, relion_MPI::DOUBLE, 0, MPITag::LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD);
            } else {
                MPI_Status status;

                for (int id_rank = 0; id_rank < node->size; id_rank++) {
                    divide_equally(op_samplings.size(), node->size, id_rank, first, last);

                    // Leader receives op_samplings_batch_packed from followers
                    if (id_rank > 0)
                        node->relion_MPI_Recv(op_samplings_batch_packed.data, (last - first + 1) * NR_LOCALSYM_PARAMETERS, relion_MPI::DOUBLE, id_rank, MPITag::LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD, status);

                    // Leader does something for itself if id_rank == 0
                    for (long int j = 0; j < Ysize(op_samplings_batch_packed); j++)
                    for (long int i = 0; i < Xsize(op_samplings_batch_packed); i++) {
                        // Beware: Ysize(op_samplings_batch_packed) is larger than (last - first + 1)
                        if (i >= 0 && i <= last - first)
                            op_samplings[i + first][CC_POS] = direct::elem(op_samplings_batch_packed, i, CC_POS);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if (node->isLeader()) {
                // TODO: For rescaled maps
                if (newdim != cropdim) {
                    for (int isamp = 0; isamp < op_samplings.size(); isamp++)
                        Localsym_scaleTranslations(op_samplings[isamp], tmp_binning_factor);
                }
                // Now translations are all unscaled.

                // TODO: add vectors together!!!
                // Update com1_float
                for (int isamp = 0; isamp < op_samplings.size(); isamp++) {
                    // Get new_com1
                    // newCom1f = Com1f + best_trans_samp - diff
                    Localsym_shiftTranslations(op_samplings[isamp], com1_float - com1_diff); // equivalently, com1_int

                    // Update v = newCom1f + ( - newR * com0)
                    Localsym_decomposeOperator(op_samplings[isamp], aa, bb, gg, dx, dy, dz, cc);
                    Matrix2D<RFLOAT> mat1 = Euler::angles2matrix(aa, bb, gg);
                    vecR3 = vectorR3(dx, dy, dz) - matmul(mat1, com0_int);
                    Localsym_composeOperator(op_samplings[isamp], aa, bb, gg, XX(vecR3), YY(vecR3), ZZ(vecR3), cc);
                }

                // Leader sorts the results
                std::stable_sort(op_samplings.begin(), op_samplings.end(), compareOperatorsByCC);

                // Leader outputs the local searches results
                // "*_cc_mask001.tmp" -> "*_cc_mask001"
                auto fn_tmp = FileName::compose(fn_info_out.withoutExtension() + "_cc_mask", imask + 1, "tmp", 3).withoutExtension();
                auto fn_searched_op_samplings = FileName::compose(fn_tmp + "_op", iop + 1, "star", 3); // "*_cc_mask001_op001.star"
                writeRelionFormatLocalSearchOperatorResults(fn_searched_op_samplings, op_samplings, angpix_image);
                std::cout << " + List of sampling points for this local symmetry operator: " << fn_searched_op_samplings << std::endl;

                // Leader updates this operator and do screen output
                op_list[imask][iop] = op_samplings[0];
                std::cout << " + Done! Refined operator: " << std::flush;
                Localsym_outputOperator(op_samplings[0], &std::cout, angpix_image);
                std::cout << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (node->isLeader()) {
        // Leader writes out new mask info file
        if (fn_info_out.getExtension() == "star") {
            writeRelionFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
        } else {
            writeDMFormatMasksAndOperators    (fn_info_out, fn_mask_list, op_list, angpix_image);
        }

        displayEmptyLine();
        #ifdef __unix__
        std::cout << " Done! New local symmetry description file: " << "\e[1m" << fn_info_out << "\e[0m" << std::endl;
        #else
        std::cout << " Done! New local symmetry description file: " << fn_info_out << std::endl;
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
