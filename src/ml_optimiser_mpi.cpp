/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#include "src/ml_optimiser_mpi.h"
#include "src/ml_optimiser.h"
#include "src/postprocessing.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_ml_optimiser.h"
#endif
#ifdef ALTCPU
#include <tbb/tbb.h>
#include "src/acc/cpu/cpu_ml_optimiser.h"
#endif
#include <stdio.h>
#include <stdlib.h>


// Like modulo, but return the base instead of zero.
// e.g. 0, 1, 2, 3, 4, 5, ... => 3, 1, 2, 3, 1, 2, ... (mod 3)
static inline int modulo_alt(int x, int y) {
    return (x - 1) % y + 1;
}

// #define PRINT_GPU_MEM_INFO
// #define DEBUG
// #define DEBUG_MPIEXP2

#ifdef TIMING
int TIMING_MPIPACK, TIMING_MPIWAIT, TIMING_MPICOMBINEDISC, TIMING_MPICOMBINENETW, TIMING_MPISLAVEWORK;
int TIMING_MPISLAVEWAIT1, TIMING_MPISLAVEWAIT2, TIMING_MPISLAVEWAIT3;
#endif


/** Read from a binary file.
 * The array must be previously resized to the correct size.
 */
template <typename T>
void readBinary(MultidimArray<T> &arr, const FileName &fn) {
    std::ifstream ifs(fn.c_str(), std::ios::in | std::ios::binary);
    if (!ifs)
        REPORT_ERROR(static_cast<std::string>("MultidimArray::read: File " + fn + " not found"));

    for (auto &x: arr)
        ifs.read(reinterpret_cast<char*>(&x), sizeof(T));

}

/** Read from a binary file, while summing to the existing array
 * The array must be previously resized to the correct size.
 */
template <typename T>
void readBinaryAndSum(MultidimArray<T> &arr, const FileName &fn) {
    std::ifstream ifs(fn.c_str(), std::ios::in | std::ios::binary);
    if (!ifs)
        REPORT_ERROR(static_cast<std::string>("MultidimArray::read: File " + fn + " not found"));

    for (auto &x : arr) {
        T y;
        ifs.read(reinterpret_cast<char*>(&y), sizeof(T));
        x += y;
    }
}

/* Write to a binary file
*/
template <typename T>
void writeBinary(const MultidimArray<T> &arr, const FileName &fn) {
    std::ofstream ofs(fn.c_str(), std::ios::out | std::ios::binary);
    if (!ofs)
        REPORT_ERROR(static_cast<std::string>("MultidimArray::write: File " + fn + " cannot be opened for output"));

    for (const auto &x : arr)
        ofs.write(reinterpret_cast<const char*>(&x), sizeof(T));
}

void determine_gpu_compatibility(int devCount) {
    // Determine how many GPUs are capable of CUDA >= 3.5
    int compatibleDevices (0);
    for (int i = 0; i < devCount; i++) {
        cudaDeviceProp deviceProp;
        HANDLE_ERROR(cudaGetDeviceProperties(&deviceProp, i));
        if (std::make_pair(deviceProp.major, deviceProp.minor) >=
            std::make_pair(CUDA_CC_MAJOR, CUDA_CC_MINOR))
        {
            compatibleDevices++;
        } else {
            // std::cout << "Rank " << node->rank << " found a " << deviceProp.name << " GPU with compute-capability " << deviceProp.major << "." << deviceProp.minor << std::endl;
        }
    }
    if (compatibleDevices == 0) {
        REPORT_ERROR("You have no GPUs compatible with RELION (CUDA-capable and compute-capability >= 3.5)");
    } else if (compatibleDevices != devCount) {
        std::cerr << "WARNING : at least one of your GPUs is not compatible with RELION (CUDA-capable and compute-capability >= 3.5)" << std::endl;
    }
}

/// FIXME: Calling this more than once will leak memory
void MlOptimiserMpi::read(int argc, char **argv) {

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::read Entering " << std::endl;
    #endif

    // Define a new MpiNode
    node = new MpiNode(argc, argv);  // Where is the delete?

    if (node->isLeader()) PRINT_VERSION_INFO();

    // First read in non-parallelisation-dependent variables
    MlOptimiser::read(argc, argv, node->rank);

    int mpi_section = parser.addSection("MPI options");
    halt_all_followers_except_this = textToInteger(parser.getOption("--halt_all_followers_except", "For debugging: keep all followers except this one waiting", "-1"));
    do_keep_debug_reconstruct_files = parser.checkOption("--keep_debug_reconstruct_files", "For debugging: keep temporary data and weight files for debug-reconstructions.");

    // Don't put any output to screen for mpi followers
    ori_verb = verb;
    if (verb) { verb = node->isLeader(); }

    // #define DEBUG_BODIES
    #ifdef DEBUG_BODIES
    verb = node->rank == 1 && ori_verb;
    #endif

    // TMP for debugging only
    // if (node->rank==1)
    //    verb = 1;
    // Possibly also read parallelisation-dependent variables here

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::read done" << std::endl;
    #endif

}

void MlOptimiserMpi::finalise() { delete node; }

void MlOptimiserMpi::initialise() {

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::initialise Entering" << std::endl;
    #endif

    // Print information about MPI nodes:
    printMpiNodesMachineNames(*node, nr_threads);
    #ifdef CUDA
    /************************************************************************/
    //Setup GPU related resources
    int devCount, deviceAffinity;
    bool is_split(false);

    if (do_gpu) {
        MPI_Status status;

        std::vector<std::string> deviceIdentifiers;
        std::vector<std::string> uniqueDeviceIdentifiers;
        std::vector<int> deviceCounts(node->size, 0);

        std::vector<std::string> uniqueHosts;
        std::vector<int>  hostId(node->size);
        std::vector<int>         ranksOnHost;
        std::vector<int>       devicesOnHost;


        // ------------------------------ FIGURE OUT GLOBAL DEVICE MAP------------------------------------------
        if (!node->isLeader()) {
            // Send device count seen by this follower
            HANDLE_ERROR(cudaGetDeviceCount(&devCount));
            determine_gpu_compatibility(devCount);
            node->relion_MPI_Send(&devCount, 1, MPI_INT, 0, MPITag::INT, MPI_COMM_WORLD);

            // Send host-name seen by this follower
            char buffer[BUFSIZ];
            int len;
            MPI_Get_processor_name(buffer, &len);
            node->relion_MPI_Send(buffer, len + 1, MPI_CHAR, 0, MPITag::IDENTIFIER, MPI_COMM_WORLD);

        } else {

            for (int follower = 1; follower < node->size; follower++) {
                // Receive device count seen by this follower
                node->relion_MPI_Recv(&devCount, 1, MPI_INT, follower, MPITag::INT, MPI_COMM_WORLD, status);
                deviceCounts[follower] = devCount;

                // Receive host-name seen by this follower
                char buffer[BUFSIZ];
                node->relion_MPI_Recv(&buffer, BUFSIZ, MPI_CHAR, follower, MPITag::IDENTIFIER, MPI_COMM_WORLD, status);
                const std::string hostName (buffer);

                // push_back unique host names and count followers on them

                // Check for novel host
                const auto search = std::find(uniqueHosts.begin(), uniqueHosts.end(), hostName);

                if (search == uniqueHosts.end()) {
                    // if unknown, push new name and counter, and assign host-id to follower
                    uniqueHosts.push_back(hostName);
                    ranksOnHost.push_back(1);
                    hostId[follower] = ranksOnHost.size() - 1;
                } else {
                    // if known, increment counter and assign host-id to follower
                    int idi = search - uniqueHosts.begin();
                    ranksOnHost[idi]++;
                    hostId[follower] = idi;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // ------------------------------ SET DEVICE INDICES BASED ON HOST DEVICE-COUNTS ------------------------------------------
        if (node->isLeader()) {
            for (int i = 0; i < uniqueHosts.size(); i++) {
                std::cout << " uniqueHost " << uniqueHosts[i] << " has " << ranksOnHost[i] << " ranks." << std::endl;
            }
            std::vector < std::vector < std::string > > allThreadIDs;
            untangleDeviceIDs(gpu_ids, allThreadIDs);

            /*if (allThreadIDs.size()==1) // if devices are specified for exactly one rank, use it for all ranks
            {
                allThreadIDs.resize(node->size - 1);
                for (int rank = 1; rank<(node->size-1); rank++)
                    allThreadIDs[rank] = allThreadIDs[0];
            }*/
            if (allThreadIDs.size() > 0 && allThreadIDs.size() < node->size - 1) {
                // If one or more devices are specified, but not for all ranks,
                // extend selection to all ranks by modulus
                int N = allThreadIDs.size();
                allThreadIDs.resize(node->size - 1);
                for (int rank = 1; rank < node->size - 1; rank++) {
                    allThreadIDs[rank] = allThreadIDs[rank % N];
                }
            }

            // Sequential initialisation of GPUs on all ranks
            for (int rank = 1; rank < node->size; rank++) {
                bool fullAutomaticMapping = true;
                bool semiAutomaticMapping = true; // possible to set fully manual for specific ranks

                const auto &thread_ids = allThreadIDs[rank - 1];

                if (
                    allThreadIDs.size() < rank ||
                    allThreadIDs.front().empty() ||
                    !std::isdigit(gpu_ids.front())
                ) {
                    std::cout << "GPU-ids not specified for this rank. "
                    << "Threads will be mapped to available devices automatically." << std::endl;
                } else {
                    fullAutomaticMapping = false;
                    if (thread_ids.size() != nr_threads) {
                        std::cout << " Follower " << rank << " will distribute threads over devices ";
                        for (int j = 0; j < thread_ids.size(); j++)
                            std::cout << " " << thread_ids[j];
                        std::cout << std::endl;
                    } else {
                        semiAutomaticMapping = false;
                        std::cout << " Using explicit indexing on follower " << rank << " to assign devices ";
                        /// FIXME: This is exactly the same as above!
                        for (int j = 0; j < thread_ids.size(); j++)
                            std::cout << " " << thread_ids[j];
                        std::cout << std::endl;
                    }
                }

                for (int i = 0; i < nr_threads; i ++) {
                    int dev_id;
                    if (semiAutomaticMapping) {
                        if (fullAutomaticMapping) {
                            int ranksOnThisHost = ranksOnHost[hostId[rank]];
                            // std::cout << rank << ", " << hostId[rank] << ", "  << hostId.size()  << std::endl;
                            dev_id = ranksOnThisHost > 1 ? deviceCounts[rank] * (((rank - 1) % ranksOnThisHost) * nr_threads + i)  / (ranksOnThisHost * nr_threads) :
                                deviceCounts[rank] * i / nr_threads;
                            // std::cout << deviceCounts[rank] << "," << (rank - 1) << "," << ranksOnThisHost << "," << nr_threads << "," << i << "," << dev_id << std::endl;
                        } else {
                            dev_id = textToInteger(thread_ids[i % thread_ids.size()].c_str());
                        }
                    } else {
                        dev_id = textToInteger(thread_ids[i].c_str());
                    }

                    std::cout << " Thread " << i << " on follower " << rank << " mapped to device " << dev_id << std::endl;
                    node->relion_MPI_Send(&dev_id, 1, MPI_INT, rank, MPITag::INT, MPI_COMM_WORLD);
                }
            }
        } else {
            for (int i = 0; i < nr_threads; i ++) {
                node->relion_MPI_Recv(&deviceAffinity, 1, MPI_INT, 0, MPITag::INT, MPI_COMM_WORLD, status);

                const auto search = std::find(cudaDevices.begin(), cudaDevices.end(), deviceAffinity);

                int bundleId;
                if (search == cudaDevices.end()) {
                    // If bundle doesn't exist on device, make a new bundle
                    cudaDevices.push_back(deviceAffinity);
                    cudaDeviceShares.push_back(1);
                    bundleId = cudaDevices.size();
                } else {
                    bundleId = search - cudaDevices.begin();
                }
                HANDLE_ERROR(cudaSetDevice(deviceAffinity));
                cudaOptimiserDeviceMap.push_back(bundleId);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (!node->isLeader()) {
            int devCount = cudaDevices.size();
            node->relion_MPI_Send(&devCount, 1, MPI_INT, 0, MPITag::INT, MPI_COMM_WORLD);

            for (int i = 0; i < devCount; i++) {
                char buffer[BUFSIZ];
                int len;
                MPI_Get_processor_name(buffer, &len);
                std::string didS (buffer, len);

                didS = "Device " + std::to_string(cudaDevices[i]) + " on " + didS;

				// std::cout << "SENDING: " << didS << std::endl;

                strcpy(buffer, didS.c_str());
                node->relion_MPI_Send(buffer, didS.length() + 1, MPI_CHAR, 0, MPITag::IDENTIFIER, MPI_COMM_WORLD);
            }

            for (int i = 0; i < devCount; i++) {
                node->relion_MPI_Recv(&cudaDeviceShares[i], 1, MPI_INT, 0, MPITag::INT, MPI_COMM_WORLD, status);
	        	// std::cout << "Received: " << bundle->rank_shared_count << std::endl;
            }
        } else {
            std::vector<std::string> deviceIdentifiers;
            std::vector<std::string> uniqueDeviceIdentifiers;
            std::vector<int> uniqueDeviceIdentifierCounts;
            std::vector<int> deviceCounts(node->size - 1,0);

            for (int follower = 1; follower < node->size; follower++) {
                int devCount;
                node->relion_MPI_Recv(&devCount, 1, MPI_INT, follower, MPITag::INT, MPI_COMM_WORLD, status);

                deviceCounts[follower - 1] = devCount;

                for (int i = 0; i < devCount; i++) {
                    char buffer[BUFSIZ];
                    node->relion_MPI_Recv(&buffer, BUFSIZ, MPI_CHAR, follower, MPITag::IDENTIFIER, MPI_COMM_WORLD, status);
                    std::string deviceIdentifier(buffer);

                    deviceIdentifiers.push_back(deviceIdentifier);

                    const auto search = std::find(uniqueDeviceIdentifiers.begin(), uniqueDeviceIdentifiers.end(), deviceIdentifier);

                    if (search == uniqueDeviceIdentifiers.end()) {
                        uniqueDeviceIdentifiers.push_back(deviceIdentifier);
                        uniqueDeviceIdentifierCounts.push_back(1);
                    } else {
                        uniqueDeviceIdentifierCounts[search - uniqueDeviceIdentifiers.begin()]++;
                    }

					// std::cout << "RECEIVING: " << deviceIdentifier << std::endl;
                }
            }

            for (int i = 0; i < uniqueDeviceIdentifiers.size(); i++) {
                if (uniqueDeviceIdentifierCounts[i] > 1) {
                    std::cout << uniqueDeviceIdentifiers[i] << " is split between " << uniqueDeviceIdentifierCounts[i] << " followers" << std::endl;
                    is_split = true;
                }
            }

            int global_did = 0;

            for (int follower = 1; follower < node->size; follower++) {
                int devCount = deviceCounts[follower - 1];

                for (int i = 0; i < devCount; i++) {

                    const auto search = std::find(
                        uniqueDeviceIdentifiers.begin(),
                        uniqueDeviceIdentifiers.end(),
                        deviceIdentifiers[global_did]
                    );

                    const int idi = search - uniqueDeviceIdentifiers.begin();

                    node->relion_MPI_Send(&uniqueDeviceIdentifierCounts[idi], 1, MPI_INT, follower, MPITag::INT, MPI_COMM_WORLD);

                    global_did++;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (do_auto_refine) {
            if (!node->isLeader()) {
                unsigned long long boxLim (10000);
                for (int i = 0; i < cudaDevices.size(); i ++) {
                    MlDeviceBundle b(this);
                    b.setDevice(cudaDevices[i]);
                    unsigned long long t = b.checkFixedSizedObjects(cudaDeviceShares[i]);
                    boxLim = std::min(t, boxLim);
                }
                node->relion_MPI_Send(&boxLim, 1, MPI_UNSIGNED_LONG_LONG, 0, MPITag::INT, MPI_COMM_WORLD);
            } else {
                unsigned long long boxLim, LowBoxLim(10000);
                for (int follower = 1; follower < node->size; follower++) {
                    node->relion_MPI_Recv(&boxLim, 1, MPI_UNSIGNED_LONG_LONG, follower, MPITag::INT, MPI_COMM_WORLD, status);
                    LowBoxLim = std::min(boxLim, LowBoxLim);
                }

                Experiment temp;
                temp.read(fn_data);

                int t_ori_size = temp.getOpticsImageSize(0);
                // t_ori_size = temp.MDopt.getValue<int>(EMDL::IMAGE_SIZE, 0);

                if (LowBoxLim < t_ori_size) {
                    anticipate_oom = true;
                    std::cerr << "\n\n                         ***WARNING***\n\n"
                    "With the current settings and hardware, you will be able to \n"
                    "use an estimated image-size of " << LowBoxLim << " pixels during the last iteration...\n\n"
                    "...but your input box-size (image_size) is however " << t_ori_size << ". This means that \n"
                    "you will likely run out of memory on the GPU(s), and will have to then re-start \n"
                    "from the last completed iteration (i.e. continue from it) *without* the use of GPUs.\n " << std::endl;

                    if (is_split) {
                        std::cerr << "You are also splitting at least one GPU between two or more mpi-followers, which \n"
                        "might be the limiting factor, since each mpi-follower that shares a GPU increases the \n"
                        "use of memory. In this case we recommend running a single mpi-follower per GPU, which \n"
                        "will still yield good performance and possibly a more stable execution. \n" << std::endl;
                    }
                    #ifdef ACC_DOUBLE_PRECISION
                    int sLowBoxLim = (float) LowBoxLim * pow(2.0, 1.0 / 3.0);
                    std::cerr << "You are also using double precison on the GPU. If you were using single precision\n"
                    "(which in all tested cases is perfectly fine), then you could use a box size of ~"  << sLowBoxLim << "." << std::endl;
                    #endif
                    std::cerr << std::endl << std::endl;
                }
            }
        }
        mymodel.do_gpu = do_gpu;
    }
    /************************************************************************/
    #endif // CUDA


    MlOptimiser::initialiseGeneral(node->rank);

    initialiseWorkLoad();

    #ifdef MKLFFT
    // Enable multi-threaded FFTW
    int success = fftw_init_threads();
    if (0 == success)
        REPORT_ERROR("Multithreaded FFTW failed to initialize");

    // And allow plans before expectation to run using allowed
    // number of threads
    fftw_plan_with_nthreads(nr_threads);
    #endif

    if (!fn_sigma.empty()) {
        // Read in sigma_noise spetrum from file DEVELOPMENTAL!!! FOR DEBUGGING ONLY...
        int idx;
        MetaDataTable MDsigma;
        MDsigma.read(fn_sigma);
        for (long int _  : MDsigma) {
            const RFLOAT val = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA2_NOISE);
            idx = MDsigma.getValue<int>(EMDL::SPECTRAL_IDX);
            if (idx < Xsize(mymodel.sigma2_noise[0])) {
                mymodel.sigma2_noise[0](idx) = val;
            }
        }
        if (idx < Xsize(mymodel.sigma2_noise[0]) - 1 && verb > 0) {
            std::cout << " WARNING: provided sigma2_noise-spectrum has fewer entries (" << idx + 1 << ")"
                " than needed (" << Xsize(mymodel.sigma2_noise[0]) << "). "
                "Set rest to zero..." << std::endl;
        }

        mydata.getNumberOfImagesPerGroup(mymodel.nr_particles_per_group);
        for (int igroup = 0; igroup< mymodel.nr_groups; igroup++) {
            // Use the same spectrum for all classes
            mymodel.sigma2_noise[igroup] =  mymodel.sigma2_noise[0];
            // We set wsum_model.sumw_group as in calculateSumOfPowerSpectraAndAverageImage
            wsum_model.sumw_group[igroup] = mymodel.nr_particles_per_group[igroup];
        }
    } else if (do_calculate_initial_sigma_noise || do_average_unaligned) {
        MultidimArray<RFLOAT> Mavg;
        // Calculate initial sigma noise model from power_class spectra of the individual images
        // This is done in parallel
        // std::cout << " Hello world1! I am node " << node->rank << " out of " << node->size <<" and my hostname= " << getenv("HOSTNAME")<< std::endl;
        calculateSumOfPowerSpectraAndAverageImage(Mavg);

        // Set sigma2_noise and Iref from averaged poser spectra and Mavg
        if (!node->isLeader())
            MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage(Mavg);
        // std::cout << " Hello world3! I am node " << node->rank << " out of " << node->size <<" and my hostname= " << getenv("HOSTNAME")<< std::endl;
    }

    MlOptimiser::initialLowPassFilterReferences();

    // Initialise the data_versus_prior ratio to get the initial current_size right
    if (iter == 0 && !do_initialise_bodies && !node->isLeader())
        mymodel.initialiseDataVersusPrior(fix_tau); // fix_tau was set in initialiseGeneral

    // std::cout << " Hello world! I am node " << node->rank << " out of " << node->size <<" and my hostname= " << getenv("HOSTNAME")<< std::endl;

    // Only leader writes out initial mymodel (do not gather metadata yet)
    const int my_nr_subsets = do_split_random_halves ? 2 : 1;
    if (node->isLeader()) {
        MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);
    } else if (node->rank <= my_nr_subsets) {
        // Only the first_follower of each subset writes model to disc
        MlOptimiser::write(DO_WRITE_SAMPLING, DONT_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, node->rank);

        bool do_warn = false;
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            long int nr_particles = mymodel.nr_particles_per_group[igroup];
            if (nr_particles < 5 && node->rank == 1) {
                std:: cout << "WARNING: There are only " << nr_particles << " particles in group " << igroup + 1;
                // Only warn for half1 to avoid messy output
                if (do_split_random_halves) {
                    std::cout << " of half-set " << node->rank;
                }
                std::cout << std::endl;
                do_warn = true;
            }
        }
        if (do_warn) {
            std:: cout << "WARNING: You may want to consider joining some micrographs into larger groups to obtain more robust noise estimates. " << std::endl;
            std:: cout << "         You can do so by using the same rlnMicrographName for particles from multiple different micrographs in the input STAR file. " << std::endl;
            std:: cout << "         It is then best to join micrographs with similar defocus values and similar apparent signal-to-noise ratios. " << std::endl;
        }
    }

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::initialise Done" << std::endl;
    #endif
}

void MlOptimiserMpi::initialiseWorkLoad() {
    if (do_split_random_halves) {
        if (node->size <= 2)
            REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 3 MPI processes are required when splitting data into random halves");
        if (node->size % 2 == 0)
            REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: the number of MPI processes must be an odd number when gold-standard seperation is applied.");
    } else if (node->size <= 1)
        REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 2 MPI processes are required, otherwise use the sequential program");

    // Get the same random number generator seed for all mpi nodes
    if (random_seed == -1) {
        if (node->isLeader()) {
            random_seed = time(nullptr);
            for (int follower = 1; follower < node->size; follower++)
                node->relion_MPI_Send(&random_seed, 1, MPI_INT, follower, MPITag::RANDOMSEED, MPI_COMM_WORLD);
        } else {
            MPI_Status status;
            node->relion_MPI_Recv(&random_seed, 1, MPI_INT, 0, MPITag::RANDOMSEED, MPI_COMM_WORLD, status);
        }
    }

    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

    // Split the data into two random halves
    if (do_split_random_halves) {
        mydata.divideParticlesInRandomHalves(random_seed, do_helical_refine);
        my_halfset = node->myRandomSubset();
    }

    if (node->isLeader()) {
        // The leader never participates in any actual work
        my_first_particle_id = 0;
        my_last_particle_id = -1;
    } else {
        if (do_split_random_halves) {
            const div_t div = std::div(node->size - 1, 2);
            const int nr_followers_halfset1 = div.quot + div.rem;
            const int nr_followers_halfset2 = div.quot;
            if (node->myRandomSubset() == 1) {
                 // Divide first half of the images
                divide_equally(mydata.numberOfParticles(1), nr_followers_halfset1, node->rank / 2, my_first_particle_id, my_last_particle_id);
            } else {
                // Divide second half of the images
                divide_equally(mydata.numberOfParticles(2), nr_followers_halfset2, node->rank / 2 - 1, my_first_particle_id, my_last_particle_id);
                my_first_particle_id += mydata.numberOfParticles(1);
                my_last_particle_id  += mydata.numberOfParticles(1);
            }
        } else {
            const int nr_followers = node->size - 1;
            divide_equally(mydata.numberOfParticles(), nr_followers, node->rank - 1, my_first_particle_id, my_last_particle_id);
        }
    }

    // Now copy particle stacks to scratch if needed
    if (!fn_scratch.empty() && !do_preread_images) {
        mydata.setScratchDirectory(fn_scratch, do_reuse_scratch, verb);

        if (!do_reuse_scratch) {
            bool also_do_ctfimage = mymodel.data_dim == 3 && do_ctf_correction;
            if (do_parallel_disc_io) {
                FileName fn_lock = mydata.initialiseScratchLock(fn_scratch, fn_out);
                // One rank after the other, all followers pass through mydata.prepareScratchDirectory()
                // This way, only the first rank on each hostname will actually copy the particle stacks
                // The rest will just update the filenames in exp_model
                bool need_to_copy = false;
                for (int inode = 0; inode < node->size; inode++) {
                    if (inode == 0) { need_to_copy = false; }
                    if (inode > 0 && inode == node->rank) {
                        // The leader removes the lock if it existed
                        need_to_copy = mydata.prepareScratchDirectory(fn_scratch, fn_lock);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }

                int myverb = node->rank == 1 && ori_verb; // Only the first follower
                if (need_to_copy)
                mydata.copyParticlesToScratch(myverb, true, also_do_ctfimage, keep_free_scratch_Gb);

                MPI_Barrier(MPI_COMM_WORLD);
                if (!need_to_copy) {
                    // This initialises nr_parts_on_scratch on non-first ranks by pretending --reuse_scratch
                    mydata.setScratchDirectory(fn_scratch, true, verb);
                    keep_scratch = true; // Setting keep_scratch for non-first ranks, to ensure that only first rank on each node deletes scratch during cleanup
                }
            } else {
                // Only the leader needs to copy the data,
                // as only it will be reading in images.
                if (node->isLeader()) {
                    mydata.prepareScratchDirectory(fn_scratch);
                    mydata.copyParticlesToScratch(1, true,  also_do_ctfimage, keep_free_scratch_Gb);
                } else {
                    mydata.copyParticlesToScratch(0, false, also_do_ctfimage, keep_free_scratch_Gb);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (!do_split_random_halves) {
        if (!node->isLeader()) {
            /* Set up a bool-array with reference responsibilities for each rank. That is;
             * if (PPrefRank[i]==true)  //on this rank
             * 		(set up reference vol and MPI_Bcast)
             * else()
             * 		(prepare to receive from MPI_Bcast)
             */
            mymodel.PPrefRank.assign(mymodel.PPref.size(), true);

            for (int i = 0; i < mymodel.PPref.size(); i++)
                mymodel.PPrefRank[i] = i % (node->size - 1) == node->rank - 1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // #define DEBUG_WORKLOAD
    #ifdef DEBUG_WORKLOAD
    std::cerr << " node->rank= " << node->rank << " my_first_particle_id= " << my_first_particle_id << " my_last_particle_id= " << my_last_particle_id << std::endl;
    #endif
}

void MlOptimiserMpi::calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg) {

    // First calculate the sum of all individual power spectra on each subset
    MlOptimiser::calculateSumOfPowerSpectraAndAverageImage(Mavg, node->rank == 1);

    if (pipeline_control_check_abort_job())
        MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

    // Now combine all weighted sums
    // Leave the option of both for a while. Then, if there are no problems with the system via files keep that one and remove the MPI version from the code
    if (combine_weights_thru_disc) {
        combineAllWeightedSumsViaFile();
    } else {
        combineAllWeightedSums();
    }
    // After introducing SGD code in Dec 2016: no longer calculate Mavg for the 2 halves separately...
    // Just calculate Mavg from AllReduce, and divide by 2 * the accumulated wsum_group
    MultidimArray<RFLOAT> Msum = MultidimArray<RFLOAT>::zeros(Mavg);
    MPI_Allreduce(Mavg.data, Msum.data, Msum.size(), relion_MPI::DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Mavg = Msum;
    // When doing random halves, the wsum_model.sumw_group[igroup], which will be used to divide Mavg by is only calculated over half the particles!
    if (do_split_random_halves) { Mavg /= 2.0; }
}

void MlOptimiserMpi::expectation() {
    const int first_follower = 1;
    int n_trials_acc;
    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_1);)
    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::expectation: Entering " << std::endl;
    #endif

    // Use maximum of 100 particles for 3D and 10 particles for 2D estimations
    n_trials_acc = std::min(
        mymodel.ref_dim == 3 && mymodel.data_dim != 3 ? 100 : 10,
        (int) mydata.numberOfParticles()
    );

    #ifdef MKLFFT
    // Allow parallel FFTW execution
    fftw_plan_with_nthreads(nr_threads);
    #endif

    // Initialise some stuff
    // A. Update current size (may have been changed to ori_size in autoAdjustAngularSampling) and resolution pointers
    updateImageSizeAndResolutionPointers();

    // B. Set the PPref Fourier transforms, initialise wsum_model, etc.
    // The leader only holds metadata, it does not set up the wsum_model (to save memory)
    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_1a);)
    if (!node->isLeader()) {
        MlOptimiser::expectationSetup();

        mydata.MDimg.clear();

        #if !defined(__APPLE__)
        malloc_trim(0);
        #endif
    }

    MPI_Barrier(MPI_COMM_WORLD);
    }

    if (!do_split_random_halves) {
        if (!node->isLeader()) {
            for (int i = 0; i < mymodel.PPref.size(); i++) {
                /* NOTE: the first follower has rank 0 on the follower communicator node->followerC,
                 *       that's why we don't have to add 1, like this;
                 *       int sender = (i)%(node->size - 1)+1 ;
                 */

                int sender = i % (node->size - 1);
                // which rank did the heavy lifting? -> sender of information
                {
                    #ifdef MPI_DEBUG
                    std::cout << "relion_MPI_Bcast debug: rank = " << node->rank << " i = " << i << " mymodel.PPref[i].data.size() = " << mymodel.PPref[i].data.size() << " sender = " << sender << " followerC = " << node->followerC << std::endl;
                    #endif
                    // Communicating over all followers means we don't have to allocate on the leader.
                    node->relion_MPI_Bcast(
                        mymodel.PPref[i].data.data,
                        mymodel.PPref[0].data.size(), relion_MPI::COMPLEX, sender, node->followerC
                    );
                    // For multibody refinement with overlapping bodies, there may be more PPrefs than bodies!
                    if (i < mymodel.nr_classes * mymodel.nr_bodies) {
                        node->relion_MPI_Bcast(
                            mymodel.tau2_class[i].data,
                            mymodel.tau2_class[0].size(), relion_MPI::DOUBLE, sender, node->followerC
                        );
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    #ifdef DEBUG
    if (node->rank == 2) {
        for (int i = 0; i < mymodel.PPref.size(); i++) {
            FileName fn_tmp;
            fn_tmp.compose("PPref_", i, "dat");
            std::ofstream ofs(fn_tmp.c_str());
            for (auto &x : mymodel.PPref[i].data)
                ofs << x.real << std::endl;
        }
    }
    #endif

    }

    MPI_Status status;

    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_2);)
    // C. Calculate expected angular errors
    // Do not do this for maxCC
    // Only the first (reconstructing) follower (i.e. from half1) calculates expected angular errors
    if (!(iter == 1 && do_firstiter_cc) && !do_skip_align && !do_skip_rotate && !do_sgd) {
        int my_nr_images, length_fn_ctf;
        if (node->isLeader()) {
            // Leader sends metadata (but not imagedata) for first 100 particles to first_follower (for calculateExpectedAngularErrors)
            MlOptimiser::getMetaAndImageDataSubset(0, n_trials_acc - 1, false);
            my_nr_images = Ysize(exp_metadata);
            node->relion_MPI_Send(&my_nr_images, 1, MPI_INT, first_follower, MPITag::JOB_REQUEST, MPI_COMM_WORLD);
            node->relion_MPI_Send(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, first_follower, MPITag::METADATA, MPI_COMM_WORLD);
            // Also send exp_fn_ctfs if necessary
            length_fn_ctf = exp_fn_ctf.length() + 1; // +1 to include \0 at the end of the string
            node->relion_MPI_Send(&length_fn_ctf, 1, MPI_INT, first_follower, MPITag::JOB_REQUEST, MPI_COMM_WORLD);
            if (length_fn_ctf > 1)
                node->relion_MPI_Send((void*) exp_fn_ctf.c_str(), length_fn_ctf, MPI_CHAR, first_follower, MPITag::METADATA, MPI_COMM_WORLD);
        } else if (node->rank == first_follower) {
            // Follower has to receive all metadata from the leader!
            node->relion_MPI_Recv(&my_nr_images, 1, MPI_INT, 0, MPITag::JOB_REQUEST, MPI_COMM_WORLD, status);
            exp_metadata.resize(my_nr_images, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
            node->relion_MPI_Recv(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, 0, MPITag::METADATA, MPI_COMM_WORLD, status);
            node->relion_MPI_Recv(&length_fn_ctf, 1, MPI_INT, 0, MPITag::JOB_REQUEST, MPI_COMM_WORLD, status);
            if (length_fn_ctf > 1) {
                char *rec_buf2;
                rec_buf2 = (char *) malloc(length_fn_ctf);
                node->relion_MPI_Recv(rec_buf2, length_fn_ctf, MPI_CHAR, 0, MPITag::METADATA, MPI_COMM_WORLD, status);
                exp_fn_ctf = rec_buf2;
                free(rec_buf2);
            }
            calculateExpectedAngularErrors(0, n_trials_acc - 1);
        }

        // The reconstructing follower Bcast acc_rottilt, acc_psi, acc_trans to all other nodes!
        node->relion_MPI_Bcast(&acc_rot,   1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&acc_trans, 1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
    }
    }

    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_3);)

    // D. Update the angular sampling (all nodes except leader)
    if (
        !node->isLeader() && (do_auto_refine || do_sgd) && iter > 1 ||
        mymodel.nr_classes > 1 && allow_coarser_samplings
    ) {
        updateAngularSampling(node->rank == 1);
    }

    // The leader needs to know about the updated parameters from updateAngularSampling
    node->relion_MPI_Bcast(&has_fine_enough_angular_sampling, 1, MPI_INT, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&nr_iter_wo_resol_gain, 1, MPI_INT, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&nr_iter_wo_large_hidden_variable_changes, 1, MPI_INT, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_classes, 1, MPI_INT, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_offsets, 1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_orientations, 1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
    if (mymodel.nr_bodies > 1) {
        for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
            node->relion_MPI_Bcast(&mymodel.keep_fixed_bodies[ibody], 1, MPI_INT, first_follower, MPI_COMM_WORLD);
        }
    }

    // Feb15,2016 - Shaoda - copy the following variables to the leader
    if (do_helical_refine && helical_sigma_distance >= 0.0) {
        node->relion_MPI_Bcast(&sampling.healpix_order,           1, MPI_INT,       first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&mymodel.sigma2_rot,               1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&mymodel.sigma2_tilt,              1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&mymodel.sigma2_psi,               1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&mymodel.sigma2_offset,            1, relion_MPI::DOUBLE, first_follower, MPI_COMM_WORLD);
        node->relion_MPI_Bcast(&mymodel.orientational_prior_mode, 1, MPI_INT,       first_follower, MPI_COMM_WORLD);
    }

    // For multi-body refinement: check if all bodies are fixed.
    // If so, don't loop over all particles, but just return.
    if (mymodel.nr_bodies > 1) {
        int all_fixed = 1;
        for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
            all_fixed *= mymodel.keep_fixed_bodies[ibody];
        }
        if (all_fixed > 0) return;
    }

    // E. All nodes, except the leader, check memory and precalculate AB-matrices for on-the-fly shifts
    if (!node->isLeader()) {
        // Check whether everything fits into memory
        MlOptimiser::expectationSetupCheckMemory(node->rank == first_follower);

        // F. Precalculate AB-matrices for on-the-fly shifts
        // Use tabulated sine and cosine values instead for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts
        if (
            do_shifts_onthefly &&
            !(do_helical_refine && !ignore_helical_symmetry) &&
            !(do_sgd && iter > 1)
        ) {
            precalculateABMatrices();
        }
    }
    // Follower 1 sends has_converged to everyone else (in particular the leader needs it!)
    node->relion_MPI_Bcast(&has_converged,         1, MPI_INT, first_follower, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&do_join_random_halves, 1, MPI_INT, first_follower, MPI_COMM_WORLD);

    }

    MultidimArray<long int> first_last_nr_images (6);

    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_4);)

    {
    ifdefTIMING(TicToc tt (timer, TIMING_EXP_4a);)
    // Wait until expected angular errors have been calculated
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    }

    // Now perform real expectation step in parallel, use an on-demand leader-follower system
    #define JOB_FIRST         (first_last_nr_images(0))
    #define JOB_LAST          (first_last_nr_images(1))
    #define JOB_NIMG          (first_last_nr_images(2))
    #define JOB_LEN_FN_IMG    (first_last_nr_images(3))
    #define JOB_LEN_FN_CTF    (first_last_nr_images(4))
    #define JOB_LEN_FN_RECIMG (first_last_nr_images(5))
    #define JOB_NPAR (JOB_LAST - JOB_FIRST + 1)

    #ifdef CUDA
    /************************************************************************/
    //GPU memory setup

    std::vector<size_t> allocationSizes;
    MPI_Barrier(MPI_COMM_WORLD);

    if (do_gpu && !node->isLeader()) {
        for (int i = 0; i < cudaDevices.size(); i ++) {
            ifdefTIMING(TicToc tt (timer, TIMING_EXP_4b);)

            MlDeviceBundle *b = new MlDeviceBundle(this);
            b->setDevice(cudaDevices[i]);
            b->setupFixedSizedObjects();
            accDataBundles.push_back((void*) b);
        }

        std::vector<unsigned> threadcountOnDevice(accDataBundles.size(), 0);

        for (int i = 0; i < cudaOptimiserDeviceMap.size(); i++) {
            std::stringstream didSs;
            didSs << "RRr" << node->rank << "t" << i;
            MlOptimiserCuda *b = new MlOptimiserCuda(this, (MlDeviceBundle*) accDataBundles[cudaOptimiserDeviceMap[i]],didSs.str().c_str());
            b->resetData();
            cudaOptimisers.push_back((void*) b);
            threadcountOnDevice[cudaOptimiserDeviceMap[i]]++;
        }

        int devCount;
        HANDLE_ERROR(cudaGetDeviceCount(&devCount));
        HANDLE_ERROR(cudaDeviceSynchronize());

        for (int i = 0; i < accDataBundles.size(); i++) {
            if (
                ((MlDeviceBundle*) accDataBundles[i])->device_id < 0 ||
                ((MlDeviceBundle*) accDataBundles[i])->device_id >= devCount
            ) {
                //std::cerr << " using device_id=" << ((MlDeviceBundle*)accDataBundles[i])->device_id << " (device no. " << ((MlDeviceBundle*)accDataBundles[i])->device_id+1 << ") which is not within the available device range" << devCount << std::endl;
                CRITICAL(ERR_GPUID);
            } else {
                HANDLE_ERROR(cudaSetDevice(((MlDeviceBundle*)accDataBundles[i])->device_id));
            }

            size_t free, total, allocationSize;
            HANDLE_ERROR(cudaMemGetInfo(&free, &total));

            free = (float) free / (float) cudaDeviceShares[i];
            size_t required_free = requested_free_gpu_memory + GPU_THREAD_MEMORY_OVERHEAD_MB * 1000 * 1000 * threadcountOnDevice[i];

            if (free < required_free) {
                printf("WARNING: Ignoring required free GPU memory amount of %zu MB, due to space insufficiency.\n", required_free / 1000000);
                allocationSize = (double) free * 0.7;
            } else {
                allocationSize = free - required_free;
            }

            if (allocationSize < 200000000) {
                printf("WARNING: The available space on the GPU after initialization (%zu MB) might be insufficient for the expectation step.\n", allocationSize/1000000);
            }
            #ifdef PRINT_GPU_MEM_INFO
            printf("INFO: Projector model size %dx%dx%d\n", (int)mymodel.PPref[0].data.xdim, (int)mymodel.PPref[0].data.ydim, (int)mymodel.PPref[0].data.zdim );
            printf("INFO: Free memory for Custom Allocator of device bundle %d of rank %d is %d MB\n", i, node->rank, (int) ( ((float)allocationSize)/1000000.0 ) );
            #endif
            allocationSizes.push_back(allocationSize);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (do_gpu && !node->isLeader()) {
        for (int i = 0; i < accDataBundles.size(); i ++)
            ((MlDeviceBundle*) accDataBundles[i])->setupTunableSizedObjects(allocationSizes[i]);
    }
    #endif // CUDA
    #ifdef ALTCPU
    /************************************************************************/
    // CPU memory setup
    MPI_Barrier(MPI_COMM_WORLD);   // Is this really necessary?
    if (do_cpu  && ! node->isLeader()) {
        unsigned nr_classes = mymodel.PPref.size();
        // Allocate Array of complex arrays for this class
        if (posix_memalign((void **) &mdlClassComplex, MEM_ALIGN, nr_classes * sizeof (std::complex<XFLOAT> *)))
            CRITICAL(RAMERR);

        // Set up XFLOAT complex array shared by all threads for each class
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            int mdlX = mymodel.PPref[iclass].data.xdim;
            int mdlY = mymodel.PPref[iclass].data.ydim;
            int mdlZ = mymodel.PPref[iclass].data.zdim;
            size_t mdlXYZ;
            if (mdlZ == 0) {
                mdlXYZ = (size_t) mdlX * (size_t) mdlY;
            } else {
                mdlXYZ = (size_t) mdlX * (size_t) mdlY * (size_t) mdlZ;
            }

            try {
                mdlClassComplex[iclass] = new std::complex<XFLOAT>[mdlXYZ];
            } catch (std::bad_alloc& ba) {
                CRITICAL(RAMERR);
            }

            std::complex<XFLOAT> *pData = mdlClassComplex[iclass];

            // Copy results into complex number array
            for (size_t i = 0; i < mdlXYZ; i ++) {
                std::complex<XFLOAT> arrayval(
                    (XFLOAT) mymodel.PPref[iclass].data.data[i].real,
                    (XFLOAT) mymodel.PPref[iclass].data.data[i].imag
                );
                pData[i] = arrayval;
            }
        }

        MlDataBundle *b = new MlDataBundle();
        b->setup(this);
        accDataBundles.push_back((void*) b);
    }  // do_cpu
    #endif // ALTCPU
    /************************************************************************/

    #ifdef MKLFFT
    // Single-threaded FFTW execution for code inside parallel processing loop
    fftw_plan_with_nthreads(1);
    #endif

    }

    long int my_nr_particles = subset_size > 0 ? subset_size : mydata.numberOfParticles();
    if (node->isLeader()) {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP_5);)

        try {
            long int progress_bar_step_size = std::max(1l, my_nr_particles / 60);
            long int prev_barstep = 0;
            long int my_first_particle = 0.0;
            long int my_last_particle = my_nr_particles - 1;
            long int my_first_particle_halfset1 = 0;
            long int my_last_particle_halfset1  = mydata.numberOfParticles(1) - 1;
            long int my_first_particle_halfset2 = mydata.numberOfParticles(1);
            long int my_last_particle_halfset2  = mydata.numberOfParticles() - 1;

            if (subset_size > 0) {
                my_last_particle_halfset1 = my_nr_particles;
                my_last_particle_halfset2 = mydata.numberOfParticles(1) + my_nr_particles;

                if (do_split_random_halves)
                    progress_bar_step_size = std::max(1l, my_nr_particles * 2 / 60);
            }

            if (verb > 0) {
                if (do_sgd) {
                    std::cout << " Stochastic Gradient Descent iteration " << iter << " of " << nr_iter;
                } else {
                    std::cout << " Expectation iteration " << iter;
                    if (!do_auto_refine)
                        std::cout << " of " << nr_iter;
                    if (my_nr_particles < mydata.numberOfParticles())
                        std::cout << " (with " << my_nr_particles << " particles)";
                }
                std::cout << std::endl;
                init_progress_bar(my_nr_particles);
            }

            // Leader distributes all packages of SomeParticles
            int nr_followers_done = 0;
            int random_halfset = 0;
            long int nr_particles_todo, nr_particles_done = 0;
            long int nr_particles_done_halfset1 = 0;
            long int nr_particles_done_halfset2 = 0;
            long int my_nr_particles_done = 0;

            while (nr_followers_done < node->size - 1) {

                pipeline_control_check_abort_job();

                // Receive a job request from a follower
                node->relion_MPI_Recv(first_last_nr_images.data, first_last_nr_images.size(), MPI_LONG, MPI_ANY_SOURCE, MPITag::JOB_REQUEST, MPI_COMM_WORLD, status);
                // Which follower sent this request?
                int this_follower = status.MPI_SOURCE;

                // #define DEBUG_MPIEXP2
                #ifdef DEBUG_MPIEXP2
                std::cerr << " MASTER RECEIVING from follower= " << this_follower
                    << " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST
                    << " JOB_NIMG= "  << JOB_NIMG << " JOB_NPAR= "  << JOB_NPAR << std::endl;
                #endif
                // The first time a follower reports it only asks for input, but does not send output of a previous processing task. In that case JOB_NIMG==0
                // Otherwise, the leader needs to receive and handle the updated metadata from the followers
                if (JOB_NIMG > 0) {
                    exp_metadata.resize(JOB_NIMG, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
                    node->relion_MPI_Recv(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, this_follower, MPITag::METADATA, MPI_COMM_WORLD, status);

                    // The leader monitors the changes in the optimal orientations and classes
                    monitorHiddenVariableChanges(JOB_FIRST, JOB_LAST);

                    // The leader then updates the mydata.MDimg table
                    MlOptimiser::setMetaDataSubset(JOB_FIRST, JOB_LAST);
                    if (verb > 0 && nr_particles_done - prev_barstep > progress_bar_step_size) {
                        prev_barstep = nr_particles_done;
                        progress_bar(nr_particles_done + JOB_NPAR);
                    }
                }

                // See which random_halfset this follower belongs to, and keep track of the number of particles that have been processed already
                if (do_split_random_halves) {
                    random_halfset = modulo_alt(this_follower, 2);
                    if (random_halfset == 1) {
                        my_nr_particles_done = nr_particles_done_halfset1;
                        nr_particles_todo = my_last_particle_halfset1 - my_first_particle_halfset1 + 1;
                        JOB_FIRST = nr_particles_done_halfset1;
                        JOB_LAST  = std::min(my_last_particle_halfset1, JOB_FIRST + nr_pool - 1);
                    } else {
                        my_nr_particles_done = nr_particles_done_halfset2;
                        nr_particles_todo = my_last_particle_halfset2 - my_first_particle_halfset2 + 1;
                        JOB_FIRST = mydata.numberOfParticles(1) + nr_particles_done_halfset2;
                        JOB_LAST  = std::min(my_last_particle_halfset2, JOB_FIRST + nr_pool - 1);
                    }
                } else {
                    random_halfset = 0;
                    my_nr_particles_done = nr_particles_done;
                    nr_particles_todo =  my_last_particle - my_first_particle + 1;
                    JOB_FIRST = nr_particles_done;
                    JOB_LAST  = std::min(my_last_particle, JOB_FIRST + nr_pool - 1);
                }

                // Now send out a new job
                if (my_nr_particles_done < nr_particles_todo) {
                    MlOptimiser::getMetaAndImageDataSubset(JOB_FIRST, JOB_LAST, !do_parallel_disc_io);
                    JOB_NIMG = Ysize(exp_metadata);
                    JOB_LEN_FN_IMG = exp_fn_img.length() + 1; // +1 to include \0 at the end of the string
                    JOB_LEN_FN_CTF = exp_fn_ctf.length() + 1;
                    JOB_LEN_FN_RECIMG = exp_fn_recimg.length() + 1;
                } else {
                    // There are no more particles in the list
                    JOB_FIRST = -1;
                    JOB_LAST = -1;
                    JOB_NIMG = 0;
                    JOB_LEN_FN_IMG = 0;
                    JOB_LEN_FN_CTF = 0;
                    JOB_LEN_FN_RECIMG = 0;
                    exp_metadata.clear();
                    exp_imagedata.clear();

                    // No more particles, this follower is done now
                    nr_followers_done++;
                }

                // std::cerr << "subset= " << subset << " half-set= " << random_halfset
                //		<< " nr_subset_particles_done= " << nr_subset_particles_done
                //		<< " nr_particles_todo=" << nr_particles_todo << " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST << std::endl;

                #ifdef DEBUG_MPIEXP2
                std::cerr << " MASTER SENDING to follower= " << this_follower
                    << " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST
                    << " JOB_NIMG= " << JOB_NIMG << " JOB_NPAR= " << JOB_NPAR << std::endl;
                #endif
                node->relion_MPI_Send(first_last_nr_images.data, first_last_nr_images.size(), MPI_LONG, this_follower, MPITag::JOB_REPLY, MPI_COMM_WORLD);

                //806 Leader also sends the required metadata and imagedata for this job
                if (JOB_NIMG > 0) {
                    node->relion_MPI_Send(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, this_follower, MPITag::METADATA, MPI_COMM_WORLD);
                    if (do_parallel_disc_io) {
                        node->relion_MPI_Send((void*) exp_fn_img.c_str(), JOB_LEN_FN_IMG, MPI_CHAR, this_follower, MPITag::METADATA, MPI_COMM_WORLD);
                        // Send filenames of images to the followers
                        if (JOB_LEN_FN_CTF > 1)
                            node->relion_MPI_Send((void*) exp_fn_ctf.c_str(), JOB_LEN_FN_CTF, MPI_CHAR, this_follower, MPITag::METADATA, MPI_COMM_WORLD);
                        if (JOB_LEN_FN_RECIMG > 1)
                            node->relion_MPI_Send((void*) exp_fn_recimg.c_str(), JOB_LEN_FN_RECIMG, MPI_CHAR, this_follower, MPITag::METADATA, MPI_COMM_WORLD);
                    } else {
                        // new in 3.1: first send the image_size of these particles (as no longer necessarily the same as mymodel.ori_size...)
                        int my_image_size = mydata.getOpticsImageSize(mydata.getOpticsGroup(JOB_FIRST, 0));
                        node->relion_MPI_Send(&my_image_size, 1, MPI_INT, this_follower, MPITag::IMAGE_SIZE, MPI_COMM_WORLD);

                        // Send imagedata to the followers
                        node->relion_MPI_Send(exp_imagedata.data, exp_imagedata.size(), relion_MPI::DOUBLE, this_follower, MPITag::IMAGE, MPI_COMM_WORLD);
                    }
                }

                // Update the total number of particles that has been done already
                nr_particles_done += JOB_NPAR;
                if (do_split_random_halves)
                    // Also update the number of particles that has been done for each subset
                    (random_halfset == 1 ?
                        nr_particles_done_halfset1 :
                        nr_particles_done_halfset2) += JOB_NPAR;
            }
        } catch (RelionError XE) {
            std::cerr << "leader encountered error: " << XE;
            MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);
        }
    } else {  // if not Leader
        ifdefTIMING(TicToc tt (timer, TIMING_EXP_6);)

        try {
            if (halt_all_followers_except_this > 0 && node->rank != halt_all_followers_except_this) {
                // Let all followers except this one sleep forever
                while (true) sleep(1000);
            }

            // Followers do the real work (The follower does not need to know to which random_halfset he belongs)
            // Start off with an empty job request
            JOB_FIRST = 0;
            JOB_LAST = -1; // So that initial nr_particles (=JOB_LAST-JOB_FIRST+1) is zero!
            JOB_NIMG = 0;
            JOB_LEN_FN_IMG = 0;
            JOB_LEN_FN_CTF = 0;
            JOB_LEN_FN_RECIMG = 0;
            node->relion_MPI_Send(first_last_nr_images.data, first_last_nr_images.size(), MPI_LONG, 0, MPITag::JOB_REQUEST, MPI_COMM_WORLD);

            while (true) {
                {
                ifdefTIMING(TicToc tt (timer, TIMING_MPISLAVEWAIT1);)
                //Receive a new bunch of particles
                node->relion_MPI_Recv(first_last_nr_images.data, first_last_nr_images.size(), MPI_LONG, 0, MPITag::JOB_REPLY, MPI_COMM_WORLD, status);
                }

                //Check whether I am done
                if (JOB_NIMG <= 0) {
                    #ifdef DEBUG
                    std::cerr <<" follower " << node->rank << " has finished expectation.." <<std::endl;
                    #endif
                    exp_imagedata.clear();
                    exp_metadata.clear();
                    break;
                } else {
                    {
                    ifdefTIMING(TicToc tt (timer, TIMING_MPISLAVEWAIT2);)

                    // Also receive the imagedata and the metadata for these images from the leader
                    exp_metadata.resize(JOB_NIMG, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
                    node->relion_MPI_Recv(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, 0, MPITag::METADATA, MPI_COMM_WORLD, status);

                    // Receive the image filenames or the exp_imagedata
                    if (do_parallel_disc_io) {
                        // Resize the exp_fn_img strings
                        char *rec_buf;
                        rec_buf = (char *) malloc(JOB_LEN_FN_IMG);
                        node->relion_MPI_Recv(rec_buf, JOB_LEN_FN_IMG, MPI_CHAR, 0, MPITag::METADATA, MPI_COMM_WORLD, status);
                        exp_fn_img = rec_buf;
                        free(rec_buf);
                        if (JOB_LEN_FN_CTF > 1) {
                            char *rec_buf2;
                            rec_buf2 = (char *) malloc(JOB_LEN_FN_CTF);
                            node->relion_MPI_Recv(rec_buf2, JOB_LEN_FN_CTF, MPI_CHAR, 0, MPITag::METADATA, MPI_COMM_WORLD, status);
                            exp_fn_ctf = rec_buf2;
                            free(rec_buf2);
                        }
                        if (JOB_LEN_FN_RECIMG > 1) {
                            char *rec_buf3;
                            rec_buf3 = (char *) malloc(JOB_LEN_FN_RECIMG);
                            node->relion_MPI_Recv(rec_buf3, JOB_LEN_FN_RECIMG, MPI_CHAR, 0, MPITag::METADATA, MPI_COMM_WORLD, status);
                            exp_fn_recimg = rec_buf3;
                            free(rec_buf3);
                        }
                    } else {
                        int mysize;
                        node->relion_MPI_Recv(&mysize, 1, MPI_INT, 0, MPITag::IMAGE_SIZE, MPI_COMM_WORLD, status);
                        // resize the exp_imagedata array
                        if (mymodel.data_dim == 3) {
                            const int n = 1 + (has_converged && do_use_reconstruct_images) + do_ctf_correction;
                            exp_imagedata.resize(n * mysize, mysize, mysize);
                        } else {
                            const int n = 1 + (has_converged && do_use_reconstruct_images);
                            exp_imagedata.resize(n * JOB_NIMG, mysize, mysize);
                        }
                        node->relion_MPI_Recv(exp_imagedata.data, exp_imagedata.size(), relion_MPI::DOUBLE, 0, MPITag::IMAGE, MPI_COMM_WORLD, status);
                    }

                    if (pipeline_control_check_abort_job())
                        MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

                    // Now process these images
                    #ifdef DEBUG_MPIEXP
                    std::cerr << " SLAVE EXECUTING node->rank= " << node->rank << " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST << std::endl;
                    #endif
                    }

                    {
                    ifdefTIMING(TicToc tt (timer, TIMING_MPISLAVEWORK);)
                    expectationSomeParticles(JOB_FIRST, JOB_LAST);
                    }

                    {
                    ifdefTIMING(TicToc tt (timer, TIMING_MPISLAVEWAIT3);)
                    // Report to the leader how many particles I have processed
                    node->relion_MPI_Send(first_last_nr_images.data, first_last_nr_images.size(), MPI_LONG, 0, MPITag::JOB_REQUEST, MPI_COMM_WORLD);
                    // Also send the metadata belonging to those
                    node->relion_MPI_Send(exp_metadata.data, exp_metadata.size(), relion_MPI::DOUBLE, 0, MPITag::METADATA, MPI_COMM_WORLD);
                    }
                }
            }
		// TODO: define MPI_COMM_SLAVES!!!!	MPI_Barrier(node->MPI_COMM_SLAVES);

            #ifdef CUDA
            if (do_gpu) {
                for (const auto &bundle : accDataBundles) {
                    ifdefTIMING(TicToc tt (timer, TIMING_EXP_7);)

                    MlDeviceBundle* b = (MlDeviceBundle*) bundle;
                    b->syncAllBackprojects();

                    for (int j = 0; j < b->backprojectors.size(); j++) {
                        const unsigned long s = wsum_model.BPref[j].data.size();
                        XFLOAT *reals   = new XFLOAT[s];
                        XFLOAT *imags   = new XFLOAT[s];
                        XFLOAT *weights = new XFLOAT[s];

                        b->backprojectors[j].getMdlData(reals, imags, weights);

                        for (unsigned long n = 0; n < s; n++) {
                            wsum_model.BPref[j].data.data[n].real += (RFLOAT) reals[n];
                            wsum_model.BPref[j].data.data[n].imag += (RFLOAT) imags[n];
                            wsum_model.BPref[j].weight.data[n]    += (RFLOAT) weights[n];
                        }

                        delete [] reals;
                        delete [] imags;
                        delete [] weights;

                        b->projectors[j].clear();
                        b->backprojectors[j].clear();
                    }

                    for (auto &plan : b->coarseProjectionPlans) plan.clear();

                }

                {
                ifdefTIMING(TicToc tt (timer, TIMING_EXP_8);)

                for (const auto &optimiser : cudaOptimisers)
                    delete (MlOptimiserCuda*) optimiser;

                cudaOptimisers.clear();

                for (const auto &bundle : accDataBundles) {

                    ((MlDeviceBundle*) bundle)->allocator->syncReadyEvents();
                    ((MlDeviceBundle*) bundle)->allocator->freeReadyAllocs();

                    #ifdef DEBUG_CUDA
                    if (((MlDeviceBundle*) bundle)->allocator->getNumberOfAllocs() != 0) {
                        printf("DEBUG_ERROR: Non-zero allocation count encountered in custom allocator between iterations.\n");
                        ((MlDeviceBundle*) bundle)->allocator->printState();
                        fflush(stdout);
                        CRITICAL(ERR_CANZ);
                    }
                    #endif

                    delete (MlDeviceBundle*) bundle;
                }

                accDataBundles.clear();

                }
            }
            #endif // CUDA
            #ifdef ALTCPU
            if (do_cpu) {
                MlDataBundle* b = (MlDataBundle*) accDataBundles[0];

                #ifdef DEBUG
                std::cerr << "Faux thread id: " << b->thread_id << std::endl;
                #endif

                for (int j = 0; j < b->backprojectors.size(); j++) {
                    const unsigned long s = wsum_model.BPref[j].data.size();
                    XFLOAT *reals   = nullptr;
                    XFLOAT *imags   = nullptr;
                    XFLOAT *weights = nullptr;

                    b->backprojectors[j].getMdlDataPtrs(reals, imags, weights);

                    for (unsigned long n = 0; n < s; n++) {
                        wsum_model.BPref[j].data.data[n].real += (RFLOAT) reals[n];
                        wsum_model.BPref[j].data.data[n].imag += (RFLOAT) imags[n];
                        wsum_model.BPref[j].weight.data[n] += (RFLOAT) weights[n];
                    }

                    b->projectors[j].clear();
                    b->backprojectors[j].clear();
                }

                for (auto &plan : b->coarseProjectionPlans) plan.clear();

                delete b;
                accDataBundles.clear();

                // Now clean up
                unsigned nr_classes = mymodel.nr_classes;
                for (int iclass = 0; iclass < nr_classes; iclass++) {
                    delete [] mdlClassComplex[iclass];
                }
                free(mdlClassComplex);

                tbbCpuOptimiser.clear();
            }
            #endif  // ALTCPU
        } catch (RelionError XE) {
            std::cerr << "follower " << node->rank << " encountered error: " << XE;
            MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);
        }
    }

    #ifdef  MKLFFT
    // Allow parallel FFTW execution to continue now that we are outside the parallel
    // portion of expectation
    fftw_plan_with_nthreads(nr_threads);
    #endif

    // Just make sure the temporary arrays are empty...
    exp_imagedata.clear();
    exp_metadata.clear();

    if (verb > 0) {
        progress_bar(my_nr_particles);
    }

    // Measure how long I have to wait for the rest
    {
    ifdefTIMING(TicToc tt (timer, TIMING_MPIWAIT);)
    node->barrierWait();
    }

    // Wait until expected angular errors have been calculated
    MPI_Barrier(MPI_COMM_WORLD);

    // All followers reset the size of their projector to zero to save memory
    if (!node->isLeader()) {
        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
            mymodel.PPref[iclass].initialiseData(0);
    }


    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::expectation: done" << std::endl;
    #endif
}

void MlOptimiserMpi::combineAllWeightedSumsViaFile() {
    {
    ifdefTIMING(TicToc tt (timer, TIMING_MPICOMBINEDISC);)
    MultidimArray<RFLOAT> Mpack;

    int nr_halfsets = do_split_random_halves ? 2 : 1;

    // Only need to combine if there are more than one followers per subset!
    if ((node->size - 1) / nr_halfsets > 1) {
        // A. First all followers pack up their wsum_model (this is done simultaneously)
        if (!node->isLeader()) {
            wsum_model.pack(Mpack); // use negative piece and nr_pieces to only make a single Mpack, i.e. do not split into multiple pieces
        }

        // B. All followers write their Mpack to disc. Do this SEQUENTIALLY to prevent heavy load on disc I/O
        for (int this_follower = 1; this_follower < node->size; this_follower++ ) {
            if (this_follower == node->rank) {
                FileName fn_pack;
                fn_pack.compose(fn_out + "_rank", node->rank, "tmp");
                writeBinary(Mpack, fn_pack);
                // std::cerr << "Rank " << node->rank <<" has written: " <<fn_pack << " sum= " <<Mpack.sum()<< std::endl;
            }
            if (!do_parallel_disc_io)
                MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // C. First follower of each subset reads all other followers' Mpack; sum; and write sum to disc
        // Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
        for (int first_follower = 1; first_follower <= nr_halfsets; first_follower++) {
            if (node->rank == first_follower) {
                for (int other_follower = first_follower + nr_halfsets; other_follower < node->size; other_follower += nr_halfsets) {
                    FileName fn_pack;
                    fn_pack.compose(fn_out + "_rank", other_follower, "tmp");
                    readBinaryAndSum(Mpack, fn_pack);
                    // std::cerr << "Follower " <<node->rank<<" has read " <<fn_pack<< " sum= " <<Mpack.sum() << std::endl;
                }
                // After adding all Mpacks together: write the sum to disc
                FileName fn_pack;
                fn_pack.compose(fn_out + "_rank", node->rank, "tmp");
                writeBinary(Mpack, fn_pack);
                //std::cerr << "Follower " <<node->rank<<" is writing total SUM in " <<fn_pack << std::endl;
            }
            if (!do_parallel_disc_io)
                MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // D. All other followers read the summed Mpack from their first_follower
        // Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
        for (int this_follower = nr_halfsets + 1; this_follower < node->size; this_follower++ ) {
            if (this_follower == node->rank) {
                // Find out who is the first follower in this subset
                int first_follower = do_split_random_halves ? modulo_alt(this_follower, 2) : 1;

                // Read the corresponding Mpack (which now contains the sum of all Mpacks)
                FileName fn_pack;
                fn_pack.compose(fn_out + "_rank", first_follower, "tmp");
                readBinary(Mpack, fn_pack);
                // std::cerr << "Rank " << node->rank << " has read: " << fn_pack << " sum= " << Mpack.sum() << std::endl;
            }
            if (!do_parallel_disc_io)
                MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // E. All followers delete their own temporary file
        // Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
        for (int this_follower = 1; this_follower < node->size; this_follower++ ) {
            if (this_follower == node->rank) {
                FileName fn_pack;
                fn_pack.compose(fn_out + "_rank", node->rank, "tmp");
                remove(fn_pack.c_str());
                //std::cerr << "Rank " << node->rank <<" has deleted: " <<fn_pack << std::endl;
            }
            if (!do_parallel_disc_io)
                MPI_Barrier(MPI_COMM_WORLD);
        }

        // F. Finally all followers unpack Msum into their wsum_model (do this simultaneously)
        if (!node->isLeader())
            wsum_model.unpack(Mpack);

    } // end if ((node->size - 1)/nr_halfsets > 1)
    }
}

void MlOptimiserMpi::combineAllWeightedSums() {
    {
    ifdefTIMING(TicToc tt (timer, TIMING_MPICOMBINENETW);)

    // Pack all weighted sums in Mpack
    MultidimArray<RFLOAT> Mpack, Msum;
    MPI_Status status;

    // First follower manually sums over all other followers of it's subset
    // And then sends results back to all those followers
    // When splitting the data into two random halves, perform two passes: one for each subset
    int nr_halfsets = do_split_random_halves ? 2 : 1;
    #ifdef DEBUG
    std::cerr << " starting combineAllWeightedSums..." << std::endl;
    #endif
    // Only combine weighted sums if there are more than one followers per subset!
    if ((node->size - 1) / nr_halfsets > 1) {
        // Loop over possibly multiple instances of Mpack of maximum size
        int piece = 0;
        int nr_pieces = 1;
        while (piece < nr_pieces) {
            // All nodes except those who will reset nr_pieces piece will pass while loop in next pass
            nr_pieces = 0;

            // First all followers pack up their wsum_model
            if (!node->isLeader()) {
                wsum_model.pack(Mpack, piece, nr_pieces);
                // The first follower(s) set Msum equal to Mpack, the others initialise to zero
                Msum = node->rank <= nr_halfsets ? Mpack :
                    MultidimArray<RFLOAT>::zeros(Mpack.xdim, Mpack.ydim, Mpack.zdim, Mpack.ndim);
            }

            // Loop through all followers: each follower sends its Msum to the next follower for its subset.
            // Each next follower sums its own Mpack to the received Msum and sends it on to the next follower
            for (int this_follower = 1; this_follower < node->size; this_follower++) {

                // Find out who is the first follower in this subset
                int first_follower = do_split_random_halves ? modulo_alt(this_follower, 2) : 1;

                // Find out who is the next follower in this subset
                int other_follower = this_follower + nr_halfsets;

                if (other_follower < node->size) {
                    if (node->rank == this_follower) {
                        #ifdef DEBUG
                        std::cerr << " AA SEND node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " other_follower= " << other_follower << std::endl;
                        #endif
                        node->relion_MPI_Send(Msum.data, Msum.size(), relion_MPI::DOUBLE, other_follower, MPITag::PACK, MPI_COMM_WORLD);
                    } else if (node->rank == other_follower) {
                        node->relion_MPI_Recv(Msum.data, Msum.size(), relion_MPI::DOUBLE, this_follower, MPITag::PACK, MPI_COMM_WORLD, status);
                        #ifdef DEBUG
                        std::cerr << " AA RECV node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " other_follower= " << other_follower << std::endl;
                        #endif
                        // Add my own Mpack to send onwards in the next step
                        Msum += Mpack;
                    }
                } else {
                    // Now this_follower has reached the last follower, which passes the final Msum to the first one (i.e. first_follower)
                    if (node->rank == this_follower) {
                        #ifdef DEBUG
                        std::cerr << " BB SEND node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " first_follower= " << first_follower << std::endl;
                        #endif
                        node->relion_MPI_Send(Msum.data, Msum.size(), relion_MPI::DOUBLE, first_follower, MPITag::PACK, MPI_COMM_WORLD);
                    } else if (node->rank == first_follower) {
                        node->relion_MPI_Recv(Msum.data, Msum.size(), relion_MPI::DOUBLE, this_follower, MPITag::PACK, MPI_COMM_WORLD, status);
                        #ifdef DEBUG
                        std::cerr << " BB RECV node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " first_follower= " << first_follower << std::endl;
                        #endif
                    }
                }
            }

            // Now loop through all followers again to pass around the Msum
            for (int this_follower = 1; this_follower < node->size; this_follower++) {
                // Find out who is the next follower in this subset
                int other_follower = this_follower + nr_halfsets;

                // Do not send to the last follower, because it already had its Msum from the cycle above, therefore subtract nr_halfsets from node->size
                if (other_follower < node->size - nr_halfsets) {
                    if (node->rank == this_follower) {
                        #ifdef DEBUG
                        std::cerr << " CC SEND node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " other_follower= " << other_follower << std::endl;
                        #endif
                        node->relion_MPI_Send(Msum.data, Msum.size(), relion_MPI::DOUBLE, other_follower, MPITag::PACK, MPI_COMM_WORLD);
                    } else if (node->rank == other_follower) {
                        node->relion_MPI_Recv(Msum.data, Msum.size(), relion_MPI::DOUBLE, this_follower, MPITag::PACK, MPI_COMM_WORLD, status);
                        #ifdef DEBUG
                        std::cerr << " CC RECV node->rank= " << node->rank << " Msum.size()= " << Msum.size()
                            << " this_follower= " << this_follower << " other_follower= " << other_follower << std::endl;
                        #endif
                    }
                }
            }

            // Finally all followers unpack Msum into their wsum_model
            if (!node->isLeader()) {
                // Subtract 1 from piece because it was incremented already...
                wsum_model.unpack(Msum, piece - 1);
            }

        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    }
}

void MlOptimiserMpi::combineWeightedSumsTwoRandomHalvesViaFile() {
    // Just sum the weighted halves from follower 1 and follower 2 and Bcast to everyone else
    if (!do_split_random_halves)
        REPORT_ERROR("MlOptimiserMpi::combineWeightedSumsTwoRandomHalvesViaFile BUG: you cannot combineWeightedSumsTwoRandomHalves if you have not split random halves");

    MultidimArray<RFLOAT> Mpack;
    FileName fn_pack = fn_out + ".tmp";

    // Everyone packs up his wsum_model (simultaneously)
    // The followers from 3 and onwards also need this in order to have the correct Mpack size to be able to read in the summed Mpack
    if (!node->isLeader())
        wsum_model.pack(Mpack);

    // Rank 2 writes it Mpack to file
    if (node->rank == 2) {
        writeBinary(Mpack, fn_pack);
        // std::cerr << "Rank " << node->rank <<" has written: " <<fn_pack << " sum= " <<Mpack.sum()<< std::endl;
    }

    // Wait until rank2 is ready
    MPI_Barrier(MPI_COMM_WORLD);

    // Now rank1 reads the file of rank2 and adds it to its own Mpack and (over)write the sum to disc
    if (node->rank == 1) {
        readBinaryAndSum(Mpack, fn_pack);
        writeBinary(Mpack, fn_pack);
    }

    // Now all followers read the sum and unpack
    // Do this sequentially in order not to have very heavy disc I/O
    for (int this_follower = 2; this_follower < node->size; this_follower++) {
        if (node->rank == this_follower)
            readBinary(Mpack, fn_pack);
        if (!do_parallel_disc_io)
            MPI_Barrier(MPI_COMM_WORLD);
    }

    // in case we're doing do_parallel_disc_io
    MPI_Barrier(MPI_COMM_WORLD);

    // Remove temporary file
    if (node->rank == 1)
        remove(fn_pack.c_str());

    // Then everyone except the leader unpacks
    if (!node->isLeader())
        wsum_model.unpack(Mpack);
}

void MlOptimiserMpi::combineWeightedSumsTwoRandomHalves() {
    // Just sum the weighted halves from follower 1 and follower 2 and Bcast to everyone else
    if (!do_split_random_halves)
        REPORT_ERROR("MlOptimiserMpi::combineWeightedSumsTwoRandomHalves BUG: you cannot combineWeightedSumsTwoRandomHalves if you have not split random halves");

    MultidimArray<RFLOAT> Mpack;
    MPI_Status status;

    int piece = 0;
    int nr_pieces = 1;
    long int pack_size;
    while (piece < nr_pieces) {
        // All nodes except those who will reset nr_pieces will pass while next time
        nr_pieces = 0;

        if (node->rank == 2) {
            wsum_model.pack(Mpack, piece, nr_pieces, false); // do not clear the model!
            node->relion_MPI_Send(Mpack.data, Mpack.size(), relion_MPI::DOUBLE, 1, MPITag::PACK, MPI_COMM_WORLD);
            Mpack.clear();
        } else if (node->rank == 1) {
            if (verb > 0) std::cout << " Combining two random halves ..." << std::endl;
            wsum_model.pack(Mpack, piece, nr_pieces);
            auto Msum = MultidimArray<RFLOAT>::zeros(Mpack.xdim, Mpack.ydim, Mpack.zdim, Mpack.ndim);
            node->relion_MPI_Recv(Msum.data, Msum.size(), relion_MPI::DOUBLE, 2, MPITag::PACK, MPI_COMM_WORLD, status);
            Msum += Mpack;
            // Unpack the sum (subtract 1 from piece because it was incremented already...)
            wsum_model.unpack(Msum, piece - 1);
            Mpack.clear();
        }
    }

    // Now follower 1 has the sum of the two halves
    MPI_Barrier(MPI_COMM_WORLD);

    // Bcast to everyone
    piece = 0;
    nr_pieces = 1;
    while (piece < nr_pieces) {

        // All nodes except those who will reset nr_pieces will pass while next time
        nr_pieces = 0;

        // The leader does not have a wsum_model!
        if (!node->isLeader()) {
            // Let's have everyone repack their Mpack, so we know the size etc
            wsum_model.pack(Mpack, piece, nr_pieces);

            // rank one sends Mpack to everyone else
            if (node->rank == 1) {
                for (int other_follower = 2; other_follower < node->size; other_follower++)
                    node->relion_MPI_Send(Mpack.data, Mpack.size(), relion_MPI::DOUBLE, other_follower, MPITag::PACK, MPI_COMM_WORLD);
            } else {
                node->relion_MPI_Recv(Mpack.data, Mpack.size(), relion_MPI::DOUBLE, 1, MPITag::PACK, MPI_COMM_WORLD, status);
            }

            // Everyone unpacks the new Mpack
            wsum_model.unpack(Mpack, piece - 1);
            Mpack.clear();
        }
    }
}

void MlOptimiserMpi::maximization() {
    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::maximization: Entering " << std::endl;
    #endif

    RFLOAT helical_twist_half1, helical_rise_half1, helical_twist_half2, helical_rise_half2;
    {
    ifdefTIMING(TicToc tt (timer, TIMING_RECONS);)

    // For multi-body refinement: check if all bodies are fixed. If so, just return
    if (mymodel.nr_bodies > 1) {
        int all_fixed = 1;
        for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
            all_fixed *= mymodel.keep_fixed_bodies[ibody];
        }
        if (all_fixed > 0) return;
    }

    if (verb > 0) {
        std::cout << " Maximization ..." << std::endl;
        init_progress_bar(mymodel.nr_classes);
    }

    helical_twist_half1 = helical_twist_half2 = helical_twist_initial;
    helical_rise_half1  = helical_rise_half2  = helical_rise_initial;

    // First reconstruct all classes in parallel
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        auto &fsc_curve = mymodel.fsc_halves_class[ibody];

        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            {
            ifdefTIMING(TicToc tt (timer, RCT_1);)
            // Either ibody or iclass can be larger than 0, never 2 at the same time!
            const int ith_recons = mymodel.nr_bodies > 1 ? ibody : iclass;
            auto &reference = mymodel.Iref[ith_recons];
            auto &gradient       = mymodel.Igrad[ith_recons];

            if (wsum_model.pdf_class[iclass] > 0.0) {
                // Parallelise: each MPI-node has a different reference
                int reconstruct_rank1;
                if (do_split_random_halves) {
                    reconstruct_rank1 = 2 * (ith_recons % ((node->size - 1) / 2)) + 1;
                } else {
                    reconstruct_rank1 = ith_recons % (node->size - 1) + 1;
                }
                if (node->rank == reconstruct_rank1) {

                    if (wsum_model.BPref[ith_recons].weight.sum() > Xmipp::epsilon) {

                        // Ideally, if do_sgd is false, old_reference should not even be in scope.
                        const auto old_reference = do_sgd ? reference : MultidimArray<RFLOAT>();

                        wsum_model.BPref[ith_recons].updateSSNRarrays(
                            mymodel.tau2_fudge_factor,
                            mymodel.tau2_class[ith_recons],
                            mymodel.sigma2_class[ith_recons],
                            mymodel.data_vs_prior_class[ith_recons],
                            mymodel.fourier_coverage_class[ith_recons],
                            fsc_curve,
                            do_split_random_halves,
                            do_join_random_halves || do_always_join_random_halves
                        );

                        if (do_external_reconstruct) {
                            FileName fn_ext_root;
                            if (iter >= 0) {
                                fn_ext_root.compose(fn_out+"_it", iter, "", 3);
                            } else {
                                fn_ext_root = fn_out;
                            }
                            if (do_split_random_halves && !do_join_random_halves) {
                                fn_ext_root += "_half1";
                            }
                            if (mymodel.nr_bodies > 1) {
                                fn_ext_root.compose(fn_ext_root + "_body",  ibody + 1,  "", 3);
                            } else {
                                fn_ext_root.compose(fn_ext_root + "_class", iclass + 1, "", 3);
                            }
                            wsum_model.BPref[ith_recons].externalReconstruct(
                                reference,
                                fn_ext_root,
                                mymodel.fsc_halves_class[ith_recons],
                                mymodel.tau2_class[ith_recons],
                                mymodel.sigma2_class[ith_recons],
                                mymodel.data_vs_prior_class[ith_recons],
                                do_join_random_halves || do_always_join_random_halves,
                                mymodel.tau2_fudge_factor,
                                node->rank == 1  // only first follower is verbose
                            );
                        } else {
                            wsum_model.BPref[ith_recons].reconstruct(
                                reference,
                                gridding_nr_iter,
                                do_map,
                                mymodel.tau2_class[ith_recons],
                                mymodel.tau2_fudge_factor,
                                wsum_model.pdf_class[iclass],
                                minres_map,
                                false
                            );
                        }

                        if (do_sgd) {

                            // Use stochastic expectation maximisation instead of SGD.
                            if (do_avoid_sgd) {
                                if (iter < sgd_ini_iter) {
                                    for (auto &x : reference) {
                                        x = std::max(0.0, x);
                                    }
                                }
                                reference -= old_reference;
                            }

                            // Now update formula: dV_kl^(n) = (mu) * dV_kl^(n-1) + (1-mu)*step_size*G_kl^(n)
                            // where G_kl^(n) is now in mymodel.Iref[iclass]!!!
                            for (long int n = 0; n < gradient.size(); n++) {
                                gradient[n] = mu * gradient[n]
                                    + (1.0 - mu) * sgd_stepsize * reference[n];
                            }

                            // update formula: V_kl^(n+1) = V_kl^(n) + dV_kl^(n)
                            reference = old_reference + gradient;
                        }
                    }

                    // Apply the body mask
                    if (mymodel.nr_bodies > 1) {
                        // 19may2015 translate the reconstruction back to its C.O.M.
                        mymodel.Iref[ibody] = translate(mymodel.Iref[ibody], mymodel.com_bodies[ibody], DONT_WRAP);

                        //#define DEBUG_BODIES_SPI
                        #ifdef DEBUG_BODIES_SPI
                        // Also write out unmasked body reconstruction
                        FileName fn_tmp;
                        fn_tmp.compose(fn_out + "_unmasked_half1_body", ibody + 1, "spi");
                        Image<RFLOAT> Itmp;
                        Itmp() = mymodel.Iref[ibody];
                        Itmp.write(fn_tmp);
                        #endif

                    }

                    // Apply local symmetry according to a list of masks and their operators
                    if (
                        fn_local_symmetry_masks    .size() >= 1 &&
                        fn_local_symmetry_operators.size() >= 1 &&
                        !has_converged
                    ) {
                        applyLocalSymmetry(reference, fn_local_symmetry_masks, fn_local_symmetry_operators);
                    }
                    // Shaoda 26 Jul 2015 - Helical symmetry local refinement
                    if (
                        do_helical_refine && !ignore_helical_symmetry &&
                        mymodel.ref_dim != 2 &&
                        iter > 1 && do_helical_symmetry_local_refinement
                    ) {
                        localSearchHelicalSymmetry(
                            reference,
                            mymodel.pixel_size,
                            particle_diameter / 2.0,
                            helical_tube_inner_diameter / 2.0,
                            helical_tube_outer_diameter / 2.0,
                            helical_z_percentage,
                            mymodel.helical_rise_min,
                            mymodel.helical_rise_max,
                            mymodel.helical_rise_inistep,
                            mymodel.helical_rise[ith_recons],
                            mymodel.helical_twist_min,
                            mymodel.helical_twist_max,
                            mymodel.helical_twist_inistep,
                            mymodel.helical_twist[ith_recons]
                        );
                    }
                    // Sjors & Shaoda Apr 2015 - Apply real space helical symmetry and real space Z axis expansion.
                    if (
                        do_helical_refine && !ignore_helical_symmetry &&
                        mymodel.ref_dim != 2 &&
                        !has_converged
                    ) {
                        imposeHelicalSymmetryInRealSpace(
                            reference,
                            mymodel.pixel_size,
                            particle_diameter / 2.0,
                            helical_tube_inner_diameter / 2.0,
                            helical_tube_outer_diameter / 2.0,
                            helical_z_percentage,
                            mymodel.helical_rise[ith_recons],
                            mymodel.helical_twist[ith_recons],
                            width_mask_edge
                        );
                    }
                    helical_rise_half1  = mymodel.helical_rise [ith_recons];
                    helical_twist_half1 = mymodel.helical_twist[ith_recons];

                    // Also perform the unregularized reconstruction
                    if (do_auto_refine && has_converged)
                        readTemporaryDataAndWeightArraysAndReconstruct(ith_recons, 1);

                }

                // In some cases there is not enough memory to reconstruct two random halves in parallel
                // Therefore the following option exists to perform them sequentially
                if (do_sequential_halves_recons)
                    MPI_Barrier(MPI_COMM_WORLD);

                // When splitting the data into two random halves, perform two reconstructions in parallel: one for each subset
                if (do_split_random_halves) {
                    int reconstruct_rank2 = 2 * (ith_recons % ((node->size - 1) / 2)) + 2;

                    if (node->rank == reconstruct_rank2) {
                        // Rank 2 does not need to do the joined reconstruction
                        if (!do_join_random_halves) {

                            const auto old_reference = do_sgd ? reference : MultidimArray<RFLOAT>();

                            wsum_model.BPref[ith_recons].updateSSNRarrays(
                                mymodel.tau2_fudge_factor,
                                mymodel.tau2_class[ith_recons],
                                mymodel.sigma2_class[ith_recons],
                                mymodel.data_vs_prior_class[ith_recons],
                                mymodel.fourier_coverage_class[ith_recons],
                                fsc_curve,
                                do_split_random_halves,
                                do_join_random_halves || do_always_join_random_halves
                            );

                            if (do_external_reconstruct) {
                                FileName fn_ext_root;
                                if (iter >= 0) {
                                    fn_ext_root.compose(fn_out+"_it", iter, "", 3);
                                } else {
                                    fn_ext_root = fn_out;
                                }
                                if (do_split_random_halves && !do_join_random_halves) {
                                    fn_ext_root += "_half2";
                                }
                                if (mymodel.nr_bodies > 1) {
                                    fn_ext_root.compose(fn_ext_root + "_body",  ibody + 1, "", 3);
                                } else {
                                    fn_ext_root.compose(fn_ext_root + "_class", iclass + 1, "", 3);
                                }
                                wsum_model.BPref[ith_recons].externalReconstruct(
                                    reference,
                                    fn_ext_root,
                                    mymodel.fsc_halves_class[ith_recons],
                                    mymodel.tau2_class[ith_recons],
                                    mymodel.sigma2_class[ith_recons],
                                    mymodel.data_vs_prior_class[ith_recons],
                                    do_join_random_halves || do_always_join_random_halves,
                                    mymodel.tau2_fudge_factor);
                            } else {
                                wsum_model.BPref[ith_recons].reconstruct(
                                    reference,
                                    gridding_nr_iter,
                                    do_map,
                                    mymodel.tau2_class[ith_recons],
                                    mymodel.tau2_fudge_factor,
                                    wsum_model.pdf_class[iclass],
                                    minres_map,
                                    false
                                );
                            }

                            if (do_sgd) {
                                // Use stochastic expectation maximisation, instead of SGD.
                                if (do_avoid_sgd) {
                                    if (iter < sgd_ini_iter) {
                                        for (auto &x : reference) {
                                            x = std::max(0.0, x);
                                        }
                                    }
                                    reference -= old_reference;
                                }

                                // Now update formula: dV_kl^(n) = (mu) * dV_kl^(n-1) + (1-mu)*step_size*G_kl^(n)
                                // where G_kl^(n) is now in mymodel.Iref[iclass]!!!
                                for (long int n = 0; n < gradient.size(); n++) {
                                    gradient[n] = mu * gradient[n] +
                                        (1.0 - mu) * sgd_stepsize * reference[n];
                                }

                                // update formula: V_kl^(n+1) = V_kl^(n) + dV_kl^(n)
                                reference = old_reference + gradient;
                            }

                            // Apply the body mask
                            if (mymodel.nr_bodies > 1) {
                                // 19 May 2015 Translate the reconstruction back to its C.O.M.
                                mymodel.Iref[ibody] = translate(mymodel.Iref[ibody], mymodel.com_bodies[ibody], DONT_WRAP);

                                #ifdef DEBUG_BODIES_SPI
                                FileName fn_tmp;
                                fn_tmp.compose(fn_out + "_unmasked_half2_body", ibody + 1, "spi");
                                Image<RFLOAT> Itmp (mymodel.Iref[ibody]);
                                Itmp.write(fn_tmp);
                                #endif
                            }

                            // Apply local symmetry according to a list of masks and their operators
                            if (fn_local_symmetry_masks.size() >= 1 && fn_local_symmetry_operators.size() >= 1 && !has_converged)
                                applyLocalSymmetry(reference, fn_local_symmetry_masks, fn_local_symmetry_operators);

                            // Shaoda 26 Jul 2015 - Helical symmetry local refinement
                            if (
                                do_helical_refine && !ignore_helical_symmetry &&
                                mymodel.ref_dim != 2 &&
                                iter > 1 && do_helical_symmetry_local_refinement
                            ) {
                                localSearchHelicalSymmetry(
                                    reference,
                                    mymodel.pixel_size,
                                    particle_diameter / 2.0,
                                    helical_tube_inner_diameter / 2.0,
                                    helical_tube_outer_diameter / 2.0,
                                    helical_z_percentage,
                                    mymodel.helical_rise_min,
                                    mymodel.helical_rise_max,
                                    mymodel.helical_rise_inistep,
                                    mymodel.helical_rise[ith_recons],
                                    mymodel.helical_twist_min,
                                    mymodel.helical_twist_max,
                                    mymodel.helical_twist_inistep,
                                    mymodel.helical_twist[ith_recons]
                                );
                            }
                            // Sjors & Shaoda Apr 2015 - Apply real space helical symmetry and real space Z axis expansion.
                            if (
                                do_helical_refine && !ignore_helical_symmetry &&
                                mymodel.ref_dim != 2 &&
                                !has_converged
                            ) {
                                imposeHelicalSymmetryInRealSpace(
                                    reference,
                                    mymodel.pixel_size,
                                    particle_diameter / 2.0,
                                    helical_tube_inner_diameter / 2.0,
                                    helical_tube_outer_diameter / 2.0,
                                    helical_z_percentage,
                                    mymodel.helical_rise[ith_recons],
                                    mymodel.helical_twist[ith_recons],
                                    width_mask_edge
                                );
                            }
                            helical_rise_half2  = mymodel.helical_rise [ith_recons];
                            helical_twist_half2 = mymodel.helical_twist[ith_recons];
                        }

                        // But rank 2 always does the unfiltered reconstruction
                        if (do_auto_refine && has_converged) {
                            readTemporaryDataAndWeightArraysAndReconstruct(ith_recons, 2);
                        }
                    }
                }
            } else {
                // When doing SGD, keep the previous reference but re-initialise the gradient to zero.
                // When not doing SGD, re-initialise the reference to zero.
                (do_sgd ? gradient : reference).initZeros();
            }
            }
            // #define DEBUG_RECONSTRUCT
            #ifdef DEBUG_RECONSTRUCT
            MPI_Barrier(MPI_COMM_WORLD);
            #endif
        }
    }

    #ifdef DEBUG
    std::cerr << "rank= " << node->rank << " has reached barrier of reconstruction" << std::endl;
    #endif
    MPI_Barrier(MPI_COMM_WORLD);

    #ifdef DEBUG
    std::cerr << "All classes have been reconstructed" << std::endl;
    #endif

    {
    ifdefTIMING(TicToc tt (timer, RCT_2);)
    // Once reconstructed, broadcast new models to all other nodes
    // This cannot be done in the reconstruction loop itself because then it will be executed sequentially
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            // either ibody or iclass can be larger than 0, never 2 at the same time!
            int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;
            auto &gradient = mymodel.Igrad[ith_recons];

            if (do_split_random_halves) {
                if (!do_join_random_halves) {
                    MPI_Status status;
                    // Make sure I am sending from the rank where the reconstruction was done (see above) to all other followers of this subset
                    // Loop twice through this, as each class was reconstructed by two different followers!!
                    int nr_halfsets = 2;
                    for (int ihalfset = 1; ihalfset <= nr_halfsets; ihalfset++) {
                        if (node->myRandomSubset() == ihalfset) {
                            int reconstruct_rank = 2 * (ith_recons % ((node->size - 1) / 2)) + ihalfset; // first pass halfset1, second pass halfset2
                            int my_first_recv = node->myRandomSubset();

                            for (int recv_node = my_first_recv; recv_node < node->size; recv_node += nr_halfsets) {
                                if (node->rank == reconstruct_rank && recv_node != node->rank) {
                                    #ifdef DEBUG
                                    std::cerr << "ihalfset= " <<ihalfset<<" Sending iclass=" <<iclass<<" Sending ibody=" <<ibody<<" from node " <<reconstruct_rank<<" to node " <<recv_node << std::endl;
                                    #endif
                                    node->relion_MPI_Send(mymodel.Iref                  [ith_recons].data, mymodel.Iref                  [ith_recons].size(), relion_MPI::DOUBLE, recv_node, MPITag::IMAGE,    MPI_COMM_WORLD);
                                    node->relion_MPI_Send(mymodel.data_vs_prior_class   [ith_recons].data, mymodel.data_vs_prior_class   [ith_recons].size(), relion_MPI::DOUBLE, recv_node, MPITag::METADATA, MPI_COMM_WORLD);
                                    node->relion_MPI_Send(mymodel.fourier_coverage_class[ith_recons].data, mymodel.fourier_coverage_class[ith_recons].size(), relion_MPI::DOUBLE, recv_node, MPITag::METADATA, MPI_COMM_WORLD);
                                    node->relion_MPI_Send(mymodel.sigma2_class          [ith_recons].data, mymodel.sigma2_class          [ith_recons].size(), relion_MPI::DOUBLE, recv_node, MPITag::RFLOAT,   MPI_COMM_WORLD);
                                    // node->relion_MPI_Send(fsc_curve.data, fsc_curve.size(), relion_MPI::DOUBLE, recv_node, MPITag::RANDOMSEED, MPI_COMM_WORLD);
                                } else if (node->rank != reconstruct_rank && node->rank == recv_node) {
                                    node->relion_MPI_Recv(mymodel.Iref                  [ith_recons].data, mymodel.Iref                  [ith_recons].size(), relion_MPI::DOUBLE, reconstruct_rank, MPITag::IMAGE,    MPI_COMM_WORLD, status);
                                    node->relion_MPI_Recv(mymodel.data_vs_prior_class   [ith_recons].data, mymodel.data_vs_prior_class   [ith_recons].size(), relion_MPI::DOUBLE, reconstruct_rank, MPITag::METADATA, MPI_COMM_WORLD, status);
                                    node->relion_MPI_Recv(mymodel.fourier_coverage_class[ith_recons].data, mymodel.fourier_coverage_class[ith_recons].size(), relion_MPI::DOUBLE, reconstruct_rank, MPITag::METADATA, MPI_COMM_WORLD, status);
                                    node->relion_MPI_Recv(mymodel.sigma2_class          [ith_recons].data, mymodel.sigma2_class          [ith_recons].size(), relion_MPI::DOUBLE, reconstruct_rank, MPITag::RFLOAT,   MPI_COMM_WORLD, status);
                                    // node->relion_MPI_Recv(fsc_curve.data, fsc_curve.size(), relion_MPI::DOUBLE, reconstruct_rank, MPITag::RANDOMSEED, MPI_COMM_WORLD, status);
                                    #ifdef DEBUG
                                    std::cerr << "ihalfset= " <<ihalfset<< " Received!!!=" <<iclass<<" ibody=" <<ibody<<" from node " <<reconstruct_rank<<" at node " <<node->rank<< std::endl;
                                    #endif
                                }
                            }
                        }
                    }
                    // No one should continue until we're all here
                    MPI_Barrier(MPI_COMM_WORLD);

                    // Now all followers have all relevant reconstructions to continue
                }
            } else {
                int reconstruct_rank = ith_recons % (node->size - 1) + 1;
                // Broadcast the reconstructed references to all other MPI nodes
                node->relion_MPI_Bcast(
                    mymodel.Iref[ith_recons].data,
                    mymodel.Iref[ith_recons].size(),
                    relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD
                );
                if (do_sgd)
                    node->relion_MPI_Bcast(
                        gradient.data,
                        gradient.size(),
                        relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD
                    );
                // Broadcast the data_vs_prior spectra to all other MPI nodes
                node->relion_MPI_Bcast(
                    mymodel.data_vs_prior_class[ith_recons].data,
                    mymodel.data_vs_prior_class[ith_recons].size(),
                    relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD
                );
                // Broadcast the fourier_coverage spectra to all other MPI nodes
                node->relion_MPI_Bcast(
                    mymodel.fourier_coverage_class[ith_recons].data,
                    mymodel.fourier_coverage_class[ith_recons].size(),
                    relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD
                    );
                // Broadcast the sigma2_class spectra to all other MPI nodes
                node->relion_MPI_Bcast(
                    mymodel.sigma2_class[ith_recons].data,
                    mymodel.sigma2_class[ith_recons].size(),
                    relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD
                );
                // Broadcast helical rise and twist of this 3D class
                if (do_helical_refine && !ignore_helical_symmetry) {
                    node->relion_MPI_Bcast(&mymodel.helical_rise [iclass], 1, relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD);
                    node->relion_MPI_Bcast(&mymodel.helical_twist[iclass], 1, relion_MPI::DOUBLE, reconstruct_rank, MPI_COMM_WORLD);
                }
            }

            // Re-set the origin (this may be lost in some cases??)
            mymodel.Iref[ith_recons].setXmippOrigin();
            if (do_sgd)
                gradient.setXmippOrigin();

            // 5 Aug 2015 - Shaoda, helical symmetry refinement, broadcast refined helical parameters
            if (iter > 1 && do_helical_refine && !ignore_helical_symmetry && do_helical_symmetry_local_refinement) {
                int reconstruct_rank1 = do_split_random_halves ?
                    2 * (ith_recons % ((node->size - 1) / 2)) + 1 :
                         ith_recons %  (node->size - 1)       + 1;
                node->relion_MPI_Bcast(&helical_twist_half1, 1, relion_MPI::DOUBLE, reconstruct_rank1, MPI_COMM_WORLD);
                node->relion_MPI_Bcast(&helical_rise_half1,  1, relion_MPI::DOUBLE, reconstruct_rank1, MPI_COMM_WORLD);

                // When splitting the data into two random halves, perform two reconstructions in parallel: one for each subset
                if (do_split_random_halves) {
                    int reconstruct_rank2 = 2 * (ith_recons % ((node->size - 1) / 2)) + 2;
                    node->relion_MPI_Bcast(&helical_twist_half2, 1, relion_MPI::DOUBLE, reconstruct_rank2, MPI_COMM_WORLD);
                    node->relion_MPI_Bcast(&helical_rise_half2,  1, relion_MPI::DOUBLE, reconstruct_rank2, MPI_COMM_WORLD);
                }
            }
        }
    }
    }

    }

    if (node->isLeader()) {
        // The leader also updates the changes in hidden variables
        { updateOverallChangesInHiddenVariables(); }
        ifdefTIMING(TicToc tt (timer, RCT_4);)
    } else {
        // Now do the maximisation of all other parameters (and calculate the tau2_class-spectra of the reconstructions
        // The lazy leader never does this: it only handles metadata and does not have the weighted sums
        { maximizationOtherParameters(); }
        ifdefTIMING(TicToc tt (timer, RCT_3);)
    }

    // The leader broadcasts the changes in hidden variables to all other nodes
    node->relion_MPI_Bcast(&current_changes_optimal_classes,          1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&current_changes_optimal_orientations,     1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&current_changes_optimal_offsets,          1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&nr_iter_wo_large_hidden_variable_changes, 1, MPI_INT,       0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_classes,         1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_offsets,         1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);
    node->relion_MPI_Bcast(&smallest_changes_optimal_orientations,    1, relion_MPI::DOUBLE, 0, MPI_COMM_WORLD);

    if (verb > 0)
        progress_bar(mymodel.nr_classes);

    if (
        do_helical_refine && !ignore_helical_symmetry &&
        mymodel.ref_dim != 2 &&
        verb > 0
    ) {
        outputHelicalSymmetryStatus(
            iter,
            helical_rise_initial,
            mymodel.helical_rise_min,
            mymodel.helical_rise_max,
            helical_twist_initial,
            mymodel.helical_twist_min,
            mymodel.helical_twist_max,
            do_helical_symmetry_local_refinement,
            mymodel.helical_rise,
            mymodel.helical_twist,
            helical_rise_half1,
            helical_rise_half2,
            helical_twist_half1,
            helical_twist_half2,
            do_split_random_halves, // TODO: && !join_random_halves ???
            std::cout
        );
    }
    if (
        do_helical_refine && !ignore_helical_symmetry &&
        do_split_random_halves
    ) {
        mymodel.helical_rise [0] = (helical_rise_half1  + helical_rise_half2)  / 2.0;
        mymodel.helical_twist[0] = (helical_twist_half1 + helical_twist_half2) / 2.0;
    }

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::maximization: done" << std::endl;
    #endif
}

void MlOptimiserMpi::joinTwoHalvesAtLowResolution() {
    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::joinTwoHalvesAtLowResolution: Entering " << std::endl;
    #endif

    if (!do_split_random_halves)
        REPORT_ERROR("BUG: you should not be in MlOptimiserMpi::joinTwoHalvesAtLowResolution!");

    // Loop over all classes (this will be just one class for now...)
    RFLOAT myres = std::max(low_resol_join_halves, 1.0 / mymodel.current_resolution);
    int lowres_r_max = ceil(mymodel.ori_size * mymodel.pixel_size / myres);

    for (int ibody = 0; ibody< mymodel.nr_bodies; ibody++) {
        #ifdef DEBUG
        std::cerr << " ibody= " << ibody << " node->rank= " << node->rank << " mymodel.keep_fixed_bodies[ibody]= " << mymodel.keep_fixed_bodies[ibody] << std::endl;
        #endif

        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        int reconstruct_rank1 = 2 * (ibody % ((node->size - 1) / 2)) + 1;
        int reconstruct_rank2 = 2 * (ibody % ((node->size - 1) / 2)) + 2;
        #ifdef DEBUG
        std::cerr << " ibody= " << ibody << " node->rank= " << node->rank << " reconstruct_rank1= " << reconstruct_rank1 << " reconstruct_rank2= " << reconstruct_rank2 << std::endl;
        #endif

        if (node->rank == reconstruct_rank1 || node->rank == reconstruct_rank2) {
            MultidimArray<Complex> lowres_data;
            MultidimArray<RFLOAT>  lowres_weight;
            wsum_model.BPref[ibody].getLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);

            if (node->rank == reconstruct_rank2) {
                MPI_Status status;

                #ifdef DEBUG
                std::cerr << " RANK2A: node->rank= " << node->rank << std::endl;
                std::cerr << "AAArank=2 lowresdata: "; lowres_data.printShape();
                #endif
                // The second follower sends its lowres_data and lowres_weight to the first follower
                node->relion_MPI_Send(lowres_data.data, 2 * lowres_data.size(),   relion_MPI::DOUBLE, reconstruct_rank1, MPITag::IMAGE,  MPI_COMM_WORLD);
                node->relion_MPI_Send(lowres_weight.data,   lowres_weight.size(), relion_MPI::DOUBLE, reconstruct_rank1, MPITag::RFLOAT, MPI_COMM_WORLD);

                // Now the first follower is calculating the average...
                #ifdef DEBUG
                std::cerr << " RANK2B: node->rank= " << node->rank << std::endl;
                std::cerr << "BBBrank=2 lowresdata: "; lowres_data.printShape();
                #endif

                // Then the second follower receives the average back from the first follower
                node->relion_MPI_Recv(lowres_data.data, 2 * lowres_data.size(),   relion_MPI::DOUBLE, reconstruct_rank1, MPITag::IMAGE,  MPI_COMM_WORLD, status);
                node->relion_MPI_Recv(lowres_weight.data,   lowres_weight.size(), relion_MPI::DOUBLE, reconstruct_rank1, MPITag::RFLOAT, MPI_COMM_WORLD, status);


            } else if (node->rank == reconstruct_rank1) {

                #ifdef DEBUG
                std::cerr << " RANK1A: node->rank= " << node->rank << std::endl;
                #endif
                std::cout << " Averaging half-reconstructions up to " << myres << " Angstrom resolution to prevent diverging orientations ..." << std::endl;
                std::cout << " Note that only for higher resolutions the FSC-values are according to the gold-standard!" << std::endl;
                MPI_Status status;
                MultidimArray<Complex> lowres_data_half2;
                MultidimArray<RFLOAT>  lowres_weight_half2;
                lowres_data_half2.resize(lowres_data);
                lowres_weight_half2.resize(lowres_weight);
                #ifdef DEBUG
                std::cerr << "AAArank=1 lowresdata: "; lowres_data.printShape();
                std::cerr << "AAArank=1 lowresdata_half2: "; lowres_data_half2.printShape();
                std::cerr << "RANK1B: node->rank= " << node->rank << std::endl;
                #endif
                // The first follower receives the average from the second follower
                node->relion_MPI_Recv(lowres_data_half2.data, 2 * lowres_data_half2.size(), relion_MPI::DOUBLE, reconstruct_rank2, MPITag::IMAGE, MPI_COMM_WORLD, status);
                node->relion_MPI_Recv(lowres_weight_half2.data, lowres_weight_half2.size(), relion_MPI::DOUBLE, reconstruct_rank2, MPITag::RFLOAT, MPI_COMM_WORLD, status);

                // The first follower calculates the average of the two lowres_data and lowres_weight arrays
                #ifdef DEBUG
                std::cerr << "BBBrank=1 lowresdata: ";       lowres_data.printShape();
                std::cerr << "BBBrank=1 lowresdata_half2: "; lowres_data_half2.printShape();
                #endif
                for (long int n = 0; n < lowres_data.size(); n++) {
                    lowres_data[n] += lowres_data_half2[n];
                    lowres_data[n] /= 2.0;
                    lowres_weight[n] += lowres_weight_half2[n];
                    lowres_weight[n] /= 2.0;
                }

                // The first follower sends the average lowres_data and lowres_weight also back to the second follower
                node->relion_MPI_Send(lowres_data.data, 2 * lowres_data.size(),   relion_MPI::DOUBLE, reconstruct_rank2, MPITag::IMAGE,  MPI_COMM_WORLD);
                node->relion_MPI_Send(lowres_weight.data,   lowres_weight.size(), relion_MPI::DOUBLE, reconstruct_rank2, MPITag::RFLOAT, MPI_COMM_WORLD);

            }

            // Now that both followers have the average lowres arrays, set them back into the backprojector
            wsum_model.BPref[ibody].setLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);
        }
    }

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::joinTwoHalvesAtLowResolution: done" << std::endl;
    #endif
}

void MlOptimiserMpi::reconstructUnregularisedMapAndCalculateSolventCorrectedFSC() {
    if (do_sgd || subset_size > 0)
        REPORT_ERROR("BUG! You cannot do solvent-corrected FSCs and subsets!");

    if (fn_mask.empty()) return;

    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        auto &fsc_curve = mymodel.fsc_halves_class[ibody];

        int reconstruct_rank1 = 2 * (ibody % ((node->size - 1) / 2)) + 1;
        int reconstruct_rank2 = 2 * (ibody % ((node->size - 1) / 2)) + 2;  // reconstruct_rank1 + 1
        if (
            mymodel.ref_dim == 3 && (
                 node->rank == reconstruct_rank1 ||
                (node->rank == reconstruct_rank2 && do_split_random_halves)
            )
        ) {
            FileName fn_root;
            if (iter >= 0) {
                fn_root.compose(fn_out + "_it", iter, "", 3);
            } else {
                fn_root = fn_out;
            }
            int random_halfset = modulo_alt(node->rank, 2);
            fn_root += "_half" + integerToString(random_halfset);
            if (mymodel.nr_bodies > 1) {
                fn_root.compose(fn_root + "_body", ibody + 1, "", 3);
            } else {
                fn_root.compose(fn_root + "_class", 1, "", 3);
            }
            BackProjector BPextra(wsum_model.BPref[ibody]);

            Image<RFLOAT> Iunreg;
            MultidimArray<RFLOAT> dummy;
            BPextra.reconstruct(Iunreg(), gridding_nr_iter, false, dummy);

            if (mymodel.nr_bodies > 1) {
                // 19may2015 translate the reconstruction back to its C.O.M.
                Iunreg() = translate(Iunreg(), mymodel.com_bodies[ibody], DONT_WRAP);
            }

            // Update header information
            Iunreg().setXmippOrigin();
            Iunreg.setStatisticsInHeader();
            Iunreg.setSamplingRateInHeader(mymodel.pixel_size);
            // And write the resulting model to disc
            Iunreg.write(fn_root + "_unfil.mrc");
        }

        // reconstruct_rank1 also sends the current_size to the leader, so that it knows where to cut the FSC to zero
        MPI_Status status;
        if (node->rank == reconstruct_rank1)
            node->relion_MPI_Send(&mymodel.current_size, 1, MPI_INT, 0, MPITag::INT, MPI_COMM_WORLD);
        if (node->rank == 0)
            node->relion_MPI_Recv(&mymodel.current_size, 1, MPI_INT, reconstruct_rank1, MPITag::INT, MPI_COMM_WORLD, status);

        MPI_Barrier(MPI_COMM_WORLD);

        if (node->rank == 0) {
            // Let's do this on the leader (hopefully it has more memory)
            if (mymodel.nr_bodies > 1) {
                std::cout << " Calculating solvent-corrected gold-standard FSC for " << ibody + 1 << "th body ..." << std::endl;
            } else {
                std::cout << " Calculating solvent-corrected gold-standard FSC ..." << std::endl;
            }

            // Read in the half-reconstruction from rank2 and perform the postprocessing-like FSC correction
            FileName fn_root1, fn_root2;
            if (iter >= 0) {
                fn_root1.compose(fn_out + "_it", iter, "", 3);
            } else {
                fn_root1 = fn_out;
            }
            if (mymodel.nr_bodies > 1) {
                fn_root1.compose(fn_root1 + "_half1_body", ibody + 1, "", 3);
                fn_root2.compose(fn_root1 + "_half2_body", ibody + 1, "", 3);
            } else {
                fn_root1.compose(fn_root1 + "_half1_class", 1, "", 3);
                fn_root2.compose(fn_root1 + "_half2_class", 1, "", 3);
            }
            auto Iunreg1 = Image<RFLOAT>::from_filename(fn_root1 + "_unfil.mrc");
            auto Iunreg2 = Image<RFLOAT>::from_filename(fn_root2 + "_unfil.mrc");
            Iunreg1().setXmippOrigin();
            Iunreg2().setXmippOrigin();

            // Calculate FSC of the unmasked maps
            const auto fsc_unmasked = getFSC(Iunreg1(), Iunreg2());

            // At what resolution does the FSC drop below 0.8 (if at all)?
            const auto search = std::find_if(fsc_unmasked.begin() + 1, fsc_unmasked.end(),
                [] (RFLOAT fsc) { return fsc < 0.8; });
            if (search == fsc_unmasked.end()) {
                std::cerr << " WARNING: FSC curve between unmasked maps never drops below 0.8. Using unmasked FSC as FSC_true... " << std::endl;
                std::cerr << " WARNING: This message should go away during the later stages of refinement!" << std::endl;

                fsc_curve = fsc_unmasked;
            } else {

                Image<RFLOAT> Imask = mymodel.nr_bodies > 1 ?
                    Image<RFLOAT>(mymodel.masks_bodies[ibody]) :
                    Image<RFLOAT>::from_filename(fn_mask);
                Imask().setXmippOrigin();

                // Calculate FSC of the masked maps
                const auto fsc_masked = getFSC(Iunreg1() * Imask(), Iunreg2() * Imask());

                // Do phase-randomised FSC-correction for the solvent mask
                const int randomize_at = search - fsc_unmasked.begin();
                if (verb > 0) {
                    std::cout.width(35);
                    std::cout << std::left << "  + randomize phases beyond: ";
                    std::cout << Xsize(Iunreg1()) * mymodel.pixel_size / randomize_at << " Angstroms" << std::endl;
                }
                // Mask randomized phases maps and calculated fsc_random_masked
                const auto fsc_random_masked = getFSC(
                    randomizePhasesBeyond(Iunreg1(), randomize_at) * Imask(),
                    randomizePhasesBeyond(Iunreg2(), randomize_at) * Imask()
                );

                fsc_curve = Postprocessing::calculateFSCtrue(
                    fsc_masked, fsc_random_masked, randomize_at
                );
            }

            // Set fsc_halves_class explicitly to zero beyond the current_size
            for (int i = mymodel.current_size / 2 + 1; i < fsc_curve.size(); i++)
                direct::elem(fsc_curve, i) = 0.0;
        }

        // Now the leader sends the FSC curve to everyone else
        node->relion_MPI_Bcast(
            fsc_curve.data, fsc_curve.size(),
            relion_MPI::DOUBLE, 0, MPI_COMM_WORLD
        );

    }
}

void MlOptimiserMpi::writeTemporaryDataAndWeightArrays() {
    if (node->rank == 1 || do_split_random_halves && node->rank == 2) {
        //#define DEBUG_RECONSTRUCTION
        #ifdef DEBUG_RECONSTRUCTION
        FileName fn_root = fn_out + "_it" + integerToString(iter, 3) + "_half" + integerToString(node->rank);
        #else
        FileName fn_root = fn_out + "_half" + integerToString(node->rank);
        #endif

        // Write out temporary arrays for all classes
        for (int ibody  = 0; ibody  < mymodel.nr_bodies;  ibody++ ) {
        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            int ith_recons = mymodel.nr_bodies > 1 ? ibody : iclass;

            FileName fn_tmp;
            if (mymodel.nr_bodies > 1) {
                fn_tmp.compose(fn_root + "_body",  ibody + 1,  "", 3);
            } else {
                fn_tmp.compose(fn_root + "_class", iclass + 1, "", 3);
            }
            if (mymodel.pdf_class[iclass] > 0.0) {
                const auto &data = wsum_model.BPref[ith_recons].data;
                Image<RFLOAT> Ireal;
                Ireal().resize(data);
                for (long int n = 0; n < data.size(); n++) {
                    Ireal()[n] = data[n].real;
                }
                Ireal.write(fn_tmp + "_data_real.mrc");
                Image<RFLOAT> Iimag;
                Iimag().resize(data);
                for (long int n = 0; n < data.size(); n++) {
                    Iimag()[n] = data[n].imag;
                }
                Iimag.write(fn_tmp + "_data_imag.mrc");
                Image<RFLOAT>(wsum_model.BPref[ith_recons].weight).write(fn_tmp + "_weight.mrc");
            }
        }
        }
     }

    MPI_Barrier(MPI_COMM_WORLD);
}

void MlOptimiserMpi::readTemporaryDataAndWeightArraysAndReconstruct(int iclass, int ihalf) {

    #ifdef DEBUG_RECONSTRUCTION
    FileName fn_root = fn_out + "_it" + integerToString(iter, 3) + "_half" + integerToString(node->rank);
    #else
    FileName fn_root = fn_out + "_half" + integerToString(ihalf);
    #endif
    if (mymodel.nr_bodies > 1) {
        fn_root.compose(fn_root + "_body",  iclass + 1, "", 3);  // Surely ibody?
    } else {
        fn_root.compose(fn_root + "_class", iclass + 1, "", 3);
    }

    // Read temporary arrays back in
    auto Ireal = Image<RFLOAT>::from_filename(fn_root + "_data_real.mrc");
    Ireal().setXmippOrigin();
    Ireal().xinit = 0;
    if (!Ireal().sameShape(wsum_model.BPref[iclass].data)) {
        wsum_model.BPref[iclass].data.printShape(std::cerr);
        Ireal().printShape(std::cerr);
        REPORT_ERROR("Incompatible size of " + fn_root + "_data_real.mrc");
    }
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Ireal()) {
        wsum_model.BPref[iclass].data.elem(i, j, k).real = Ireal().elem(i, j, k);
    }

    auto Iimag = Image<RFLOAT>::from_filename(fn_root + "_data_imag.mrc");
    Iimag().setXmippOrigin();
    Iimag().xinit = 0;
    if (!Iimag().sameShape(wsum_model.BPref[iclass].data)) {
        wsum_model.BPref[iclass].data.printShape(std::cerr);
        Iimag().printShape(std::cerr);
        REPORT_ERROR("Incompatible size of " + fn_root + "_data_imag.mrc");
    }
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Iimag()) {
        wsum_model.BPref[iclass].data.elem(i, j, k).imag = Iimag().elem(i, j, k);
    }

    auto Iweight = Image<RFLOAT>::from_filename(fn_root + "_weight.mrc");
    Iweight().setXmippOrigin();
    Iweight().xinit = 0;
    if (!Iweight().sameShape(wsum_model.BPref[iclass].weight)) {
        wsum_model.BPref[iclass].weight.printShape(std::cerr);
        Iweight().printShape(std::cerr);
        REPORT_ERROR("Incompatible size of " + fn_root + "_weight.mrc");
    }
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Iweight()) {
        wsum_model.BPref[iclass].weight.elem(i, j, k) = Iweight().elem(i, j, k);
    }

    // Now perform the unregularized reconstruction
    Image<RFLOAT> Iunreg;
    MultidimArray<RFLOAT> dummy;
    wsum_model.BPref[iclass].reconstruct(Iunreg(), gridding_nr_iter, false, dummy);

    if (mymodel.nr_bodies > 1) {
        // 19may2015 translate the reconstruction back to its C.O.M.
        Iunreg() = translate(Iunreg(), mymodel.com_bodies[iclass], DONT_WRAP);
    }

    // Update header information
    Iunreg.setStatisticsInHeader();
    Iunreg.setSamplingRateInHeader(mymodel.pixel_size);
    // And write the resulting model to disc
    Iunreg.write(fn_root + "_unfil.mrc");

    // remove temporary arrays from the disc
    #ifndef DEBUG_RECONSTRUCTION
    if (!do_keep_debug_reconstruct_files) {
        remove((fn_root + "_data_real.mrc").c_str());
        remove((fn_root + "_data_imag.mrc").c_str());
        remove((fn_root +    "_weight.mrc").c_str());
    }
    #endif
}

void MlOptimiserMpi::compareTwoHalves() {
    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::compareTwoHalves: Entering " << std::endl;
    #endif

    if (!do_split_random_halves)
        REPORT_ERROR("ERROR: you should not be in MlOptimiserMpi::compareTwoHalves if !do_split_random_halves");

    if (mymodel.nr_classes > 1)
        REPORT_ERROR("ERROR: you should not be in MlOptimiserMpi::compareTwoHalves if mymodel.nr_classes > 1");

    // Only do gold-standard FSC comparisons for single-class refinements
    // TODO: Rank 0 and 1 do all bodies sequentially here... That parallelisation could be improved...
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;
        
        auto &fsc_curve = mymodel.fsc_halves_class[ibody];

        // The first two followers calculate the sum of the downsampled average of all bodies
        if (node->rank == 1 || node->rank == 2) {
            MultidimArray<Complex> avg1 = wsum_model.BPref[ibody].getDownsampledAverage();

            //#define DEBUG_FSC
            #ifdef DEBUG_FSC
            MultidimArray<Complex> avg;
            MultidimArray<RFLOAT> Mavg;
            if (mymodel.ref_dim == 2) {
                Mavg.resize(mymodel.ori_size, mymodel.ori_size);
            } else {
                Mavg.resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
            }

            FourierTransformer transformer_debug;
            transformer_debug.setReal(Mavg);
            avg.alias(transformer_debug.getFourier());
            wsum_model.BPref[0].decenter(avg1, avg, wsum_model.BPref[0].r_max * wsum_model.BPref[0].r_max);
            transformer_debug.inverseFourierTransform();
            FileName fnt;
            fnt.compose("downsampled_avg_half", node->rank, "spi");
            Image<RFLOAT> It;
            CenterFFT(Mavg, true);
            It() = Mavg;
            It.write(fnt);
            #endif

            if (node->rank == 2) {
                // The second follower sends its average to the first follower
                node->relion_MPI_Send(avg1.data, 2 * avg1.size(), relion_MPI::DOUBLE, 1, MPITag::IMAGE, MPI_COMM_WORLD);
            } else if (node->rank == 1) {

                std::cout << " Calculating gold-standard FSC";
                if (mymodel.nr_bodies > 1) {
                    std::cout << " for " << ibody + 1 << "th body";
                }
                std::cout << " ..." << std::endl;

                // The first follower receives the average from the second follower and calculates the FSC between them
                MPI_Status status;
                MultidimArray<Complex> avg2 (avg1.xdim, avg1.ydim, avg1.zdim, avg1.ndim);
                node->relion_MPI_Recv(avg2.data, 2 * avg2.size(), relion_MPI::DOUBLE, 2, MPITag::IMAGE, MPI_COMM_WORLD, status);
                fsc_curve = wsum_model.BPref[ibody].calculateDownSampledFourierShellCorrelation(avg1, avg2);
            }

        }

        // Now follower 1 sends the fsc curve to everyone else
        node->relion_MPI_Bcast(fsc_curve.data, fsc_curve.size(), relion_MPI::DOUBLE, 1, MPI_COMM_WORLD);

    }

    #ifdef DEBUG
    std::cerr << "MlOptimiserMpi::compareTwoHalves: done" << std::endl;
    #endif
}

void MlOptimiserMpi::iterate() {
    #ifdef TIMING
    // MPI-specific timing stuff goes here...
    TIMING_MPIWAIT        = timer.setNew("mpiWaitEndOfExpectation");
    TIMING_MPICOMBINEDISC = timer.setNew("mpiCombineThroughDisc");
    TIMING_MPICOMBINENETW = timer.setNew("mpiCombineThroughNetwork");
    TIMING_MPISLAVEWORK   = timer.setNew("mpiFollowerWorking");
    TIMING_MPISLAVEWAIT1  = timer.setNew("mpiFollowerWaiting1");
    TIMING_MPISLAVEWAIT2  = timer.setNew("mpiFollowerWaiting2");
    TIMING_MPISLAVEWAIT3  = timer.setNew("mpiFollowerWaiting3");
    #endif

    // Launch threads etc.
    MlOptimiser::iterateSetup();

    // Initialize the current resolution
    updateCurrentResolution();

    for (iter = iter + 1; iter <= nr_iter; iter++) {
        {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP);)

        // Update subset_size
        updateSubsetSize(node->isLeader());

        // Randomly take different subset of the particles each time we do a new "iteration" in SGD
        if (random_seed > 0) {
            mydata.randomiseParticlesOrder(random_seed+iter, do_split_random_halves,  subset_size < mydata.numberOfParticles() );
        } else if (verb > 0) {
            std::cerr << " WARNING: skipping randomisation of particle order because random_seed equals zero..." << std::endl;
        }

        // Nobody can start the next iteration until everyone has finished
        MPI_Barrier(MPI_COMM_WORLD);

        // Only first follower checks for convergence and prints stats to the stdout
        if (do_auto_refine) { checkConvergence(node->rank == 1); }

        expectation();
        #ifdef DEBUG
        std::cerr << " finished expectation..." << std::endl;
        #endif

        MPI_Barrier(MPI_COMM_WORLD);

        if (do_skip_maximization) {
            // Only write data.star file and break from the iteration loop
            if (node->isLeader()) {
                // The leader only writes the data file (he's the only one who has and manages these data!)
                iter = -1; // write output file without iteration number
                MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);
                if (verb > 0)
                    std::cout << " Auto-refine: Skipping maximization step, so stopping now... " << std::endl;
            }
            break;
        }

        // Now combine all weighted sums
        // Leave the option to both for a while. Then, if there are no problems with the system via files keep that one and remove the MPI version from the code
        #ifdef DEBUG
        std::cerr << " before combineAllWeightedSums..." << std::endl;
        #endif
        if (combine_weights_thru_disc) {
            combineAllWeightedSumsViaFile();
        } else {
            combineAllWeightedSums();
        }
        #ifdef DEBUG
        std::cerr << " after combineAllWeightedSums..." << std::endl;
        #endif

        MPI_Barrier(MPI_COMM_WORLD);

        // Sjors & Shaoda Apr 2015
        // This function does enforceHermitianSymmetry, applyHelicalSymmetry and applyPointGroupSymmetry sequentially.
        // First it enforces Hermitian symmetry to the back-projected Fourier 3D matrix.
        // Then helical symmetry is applied in Fourier space. It does rise and twist for all asymmetrical units in Fourier space.
        // Finally it applies point group symmetry (such as Cn, ...).
        // DEBUG
        if (node->isLeader() && verb > 0) {
            if (
                do_helical_refine && !ignore_helical_symmetry && mymodel.ref_dim != 2
            ) {
                if (mymodel.helical_nr_asu > 1)
                    std::cout << " Applying helical symmetry from the last iteration for all asymmetrical units in Fourier space..." << std::endl;
                if (iter > 1 && do_helical_symmetry_local_refinement) {
                    std::cout << " Refining helical symmetry in real space..." << std::endl;
                    std::cout << " Applying refined helical symmetry in real space..." << std::endl;
                } else {
                    std::cout << " Applying helical symmetry from the last iteration in real space..." << std::endl;
                }
            }
        }
        symmetriseReconstructions();

        if (verb > 0 && node->isLeader() && fn_local_symmetry_masks.size() >= 1 && fn_local_symmetry_operators.size() >= 1)
            std::cout << " Applying local symmetry in real space according to " << fn_local_symmetry_operators.size() << " operators..." << std::endl;

        // Write out data and weight arrays to disc in order to also do an unregularized reconstruction
        #ifndef DEBUG_RECONSTRUCTION
        if (do_auto_refine && has_converged || do_keep_debug_reconstruct_files)
        #endif
            writeTemporaryDataAndWeightArrays();

        // Inside iterative refinement: do FSC-calculation BEFORE the solvent flattening, otherwise over-estimation of resolution
        // anyway, now that this is done inside BPref, there would be no other way...
        if (do_split_random_halves) {

            // For asymmetric molecules, join 2 half-reconstructions at the lowest resolutions to prevent them from diverging orientations
            if (low_resol_join_halves > 0.0) { joinTwoHalvesAtLowResolution(); }

            #ifdef DEBUG
            std::cerr << " before compareHalves..." << std::endl;
            #endif
            // Sjors 27-oct-2015
            // Calculate gold-standard FSC curve
            if (do_phase_random_fsc && (fn_mask != "None" || mymodel.nr_bodies > 1)) {
                reconstructUnregularisedMapAndCalculateSolventCorrectedFSC();
            } else {
                compareTwoHalves();
            }
            #ifdef DEBUG
            std::cerr << " after compareHalves..." << std::endl;
            #endif

            // For automated sampling procedure
            if (!node->isLeader()) {
                // the leader does not have the correct mymodel.current_size, it only handles metadata!

                // Check that incr_size is at least the number of shells as between FSC=0.5 and FSC=0.143
                for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

                    if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
                        continue;

                    auto &fsc_curve = mymodel.fsc_halves_class[ibody];

                    int fsc05   = -1;
                    int fsc0143 = -1;
                    for (long int i = 0; i < Xsize(fsc_curve); i++) {
                        if (direct::elem(fsc_curve, i) < 0.5 && fsc05 < 0)
                        { fsc05 = i; }
                        if (direct::elem(fsc_curve, i) < 0.143 && fsc0143 < 0)
                        { fsc0143 = i; }
                    }

                    // At least fsc05 - fsc0143 + 5 shells as incr_size
                    incr_size = std::max(incr_size, fsc0143 - fsc05 + 5);
                    if (!has_high_fsc_at_limit)
                        has_high_fsc_at_limit = (direct::elem(fsc_curve, mymodel.current_size / 2 - 1) > 0.2);
                }
            }

            // Upon convergence join the two random halves
            if (do_join_random_halves || do_always_join_random_halves) {
                if (combine_weights_thru_disc) {
                    combineWeightedSumsTwoRandomHalvesViaFile();
                } else {
                    combineWeightedSumsTwoRandomHalves();
                }
            }
        }

        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_MAX);)

        maximization();

        // Make sure all nodes have the same resolution, set the data_vs_prior array from half1 also for half2
        // Because there is an if-statement on ave_Pmax to set the image size, also make sure this one is the same for both halves
        if (do_split_random_halves) {
            node->relion_MPI_Bcast(&mymodel.ave_Pmax, 1, relion_MPI::DOUBLE, 1, MPI_COMM_WORLD);
            if (mymodel.nr_bodies > 1) {
                // Multiple bodies may have been reconstructed on rank other than 1!
                for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
                    int reconstruct_rank1 = 2 * (ibody % ( (node->size - 1) / 2)) + 1;
                    node->relion_MPI_Bcast(
                        mymodel.data_vs_prior_class[ibody].data,
                        mymodel.data_vs_prior_class[ibody].size(),
                        relion_MPI::DOUBLE, reconstruct_rank1, MPI_COMM_WORLD
                    );
                }

            } else {
                for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
                    node->relion_MPI_Bcast(
                        mymodel.data_vs_prior_class[iclass].data,
                        mymodel.data_vs_prior_class[iclass].size(),
                        relion_MPI::DOUBLE, 1, MPI_COMM_WORLD
                    );
            }
        }

        }

        MPI_Barrier(MPI_COMM_WORLD);

        {
        ifdefTIMING(TicToc tt (timer, TIMING_ITER_HELICALREFINE);)
        if (node->isLeader()) {
            if (
                do_helical_refine && !do_skip_align && !do_skip_rotate &&
                mymodel.ref_dim == 3
            ) {
                updatePriorsForHelicalReconstruction(
                    mydata.MDimg,
                    helical_sigma_distance * ((RFLOAT) mymodel.ori_size),
                    mymodel.helical_rise, mymodel.helical_twist,
                    helical_nstart,
                    mymodel.data_dim == 3,
                    do_auto_refine,
                    mymodel.sigma2_rot, mymodel.sigma2_tilt, mymodel.sigma2_psi,
                    mymodel.sigma2_offset,
                    helical_keep_tilt_prior_fixed,
                    verb
                );
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Directly use fn_out, without "_it" specifier, so unmasked refs will be overwritten at every iteration
        if (do_write_unmasked_refs && node->rank == 1)
        mymodel.write(fn_out + "_unmasked", sampling, false, true);

        // Mask the reconstructions to get rid of noisy solvent areas
        // Skip masking upon convergence (for validation purposes)
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_SOLVFLAT);)
        if (do_solvent && !has_converged) solventFlatten();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_UPDATERES);)
        // Re-calculate the current resolution, do this before writing to get the correct values in the output files
        updateCurrentResolution();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_ITER_WRITE);)
        // If we are joining random halves, then do not write an optimiser file so that it cannot be restarted!
        bool do_write_optimiser = !do_join_random_halves;
        // Write out final map without iteration number in the filename
        if (do_join_random_halves) { iter = -1; }

        if (node->rank == 1 || (node->rank == 2 && do_split_random_halves && !do_join_random_halves)) {
            // Only the first_follower of each subset writes model to disc (do not write the data.star file, only leader will do this)
            MlOptimiser::write(DO_WRITE_SAMPLING, DONT_WRITE_DATA, do_write_optimiser, DO_WRITE_MODEL, node->rank);
        } else if (node->isLeader()) {
            // The leader only writes the data file (he's the only one who has and manages these data!)
            MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);
        }
        }

        if (do_auto_refine && has_converged) {
            if (verb > 0) {
                std::cout << " Auto-refine: Refinement has converged, stopping now... " << std::endl;

                if (mymodel.nr_bodies == 1) {
                    std::cout << " Auto-refine: + Final reconstruction from all particles is saved as: " <<  fn_out << "_class001.mrc" << std::endl;
                } else {
                    std::cout << " Auto-refine: + Final reconstructions of each body from all particles are saved as " <<  fn_out << "_bodyNNN.mrc, where NNN is the body number" << std::endl;
                }
                std::cout << " Auto-refine: + Final model parameters are stored in: " << fn_out << "_model.star" << std::endl;
                std::cout << " Auto-refine: + Final data parameters are stored in: "  << fn_out << "_data.star"  << std::endl;

                if (mymodel.tau2_fudge_factor > 1.0) {
                    std::cout << " Auto-refine: + SEVERE WARNING: Because you used a tau2_fudge of " << mymodel.tau2_fudge_factor << " your resolution during this refinement will be inflated!" << std::endl;
                    std::cout << " Auto-refine: + SEVERE WARNING: You have to run a postprocessing on the unfil.mrc maps to get a gold-standard resolution estimate!"  << std::endl;
                } else if (do_phase_random_fsc) {
                    std::cout << " Auto-refine: + Final resolution (already with masking) is: " << 1./mymodel.current_resolution << std::endl;
                } else {
                    std::cout << " Auto-refine: + Final resolution (without masking) is: " << 1./mymodel.current_resolution << std::endl;
                    std::cout << " Auto-refine: + But you may want to run relion_postprocess to mask the unfil.mrc maps and calculate a higher resolution FSC" << std::endl;
                }

                if (acc_rot > 10.0) {
                    std::cout << " Auto-refine: + WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles!" << std::endl;
                    std::cout << " Auto-refine: + WARNING: This has been observed to lead to spurious FSC curves, so be VERY wary of inflated resolution estimates..." << std::endl;
                    std::cout << " Auto-refine: + WARNING: You most probably do NOT want to publish these results!" << std::endl;
                    std::cout << " Auto-refine: + WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
                }

                if (do_use_reconstruct_images)
                    std::cout << " Auto-refine: + Used rlnReconstructImageName images for final reconstruction. Ignore filtered map, and only assess the unfiltered half-reconstructions!" << std::endl;
            }
            break;
        }

        #ifdef TIMING
        // Only first follower prints it timing information
        if (node->rank == 1) timer.printTimes(false);
        #endif

        if (do_auto_refine && has_converged) break;

    }

    // Hopefully this barrier will prevent some bus errors
    MPI_Barrier(MPI_COMM_WORLD);

    // delete threads etc.
    MlOptimiser::iterateWrapUp();
    MPI_Barrier(MPI_COMM_WORLD);
}
