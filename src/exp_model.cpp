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
#include "src/exp_model.h"
#include <sys/statvfs.h>

long int Experiment::numberOfParticles(int random_subset) {
    if (random_subset == 0) return particles.size();
    if (random_subset == 1) return nr_particles_subset1;
    if (random_subset == 2) return nr_particles_subset2;
    REPORT_ERROR("ERROR: Experiment::numberOfParticles invalid random_subset: " + integerToString(random_subset));
}

// Get the total number of images in a given particle
long int Experiment::numberOfImagesInParticle(long int part_id) {
    return particles[part_id].images.size();
}

long int Experiment::numberOfMicrographs() {
    return micrographs.size();
}

long int Experiment::numberOfGroups() {
    return groups.size();
}

int Experiment::numberOfOpticsGroups() {
    return obsModel.numberOfOpticsGroups();
}

bool Experiment::hasCtfPremultiplied() {
    for (int og = 0; og < numberOfOpticsGroups(); og++)
        if (obsModel.getCtfPremultiplied(og)) return true;

    return false;
}

RFLOAT Experiment::getOpticsPixelSize(int optics_group) {
    return obsModel.getPixelSize(optics_group);
}

int Experiment::getOpticsImageSize(int optics_group) {
    return obsModel.getBoxSize(optics_group);
}

long int Experiment::getMicrographId(long int part_id, int img_id) {
    return particles[part_id].images[img_id].micrograph_id;
}

long int Experiment::getGroupId(long int part_id, int img_id) {
    return particles[part_id].images[img_id].group_id;
}

int Experiment::getOpticsGroup(long part_id, int img_id) {
    return particles[part_id].images[img_id].optics_group;
}

int Experiment::getRandomSubset(long int part_id) {
    return particles[part_id].random_subset;
}

int Experiment::getOriginalImageId(long part_id, int img_id) {
    return particles[part_id].images[img_id].id;
}
RFLOAT Experiment::getImagePixelSize(long int part_id, int img_id) {
    const int optics_group = particles[part_id].images[img_id].optics_group;
    return obsModel.getPixelSize(optics_group);
}

void Experiment::getNumberOfImagesPerGroup(std::vector<long int> &nr_particles_per_group) {
    nr_particles_per_group.resize(groups.size());

    for (long int part_id = 0; part_id < particles.size(); part_id++)
        for (int img_id = 0; img_id < particles[part_id].images.size(); img_id++)
            nr_particles_per_group[particles[part_id].images[img_id].group_id] += 1;
}

MetaDataTable Experiment::getMetaDataImage(long int part_id, int img_id) {
    MetaDataTable result;
    result.addObject(MDimg.getObject(getOriginalImageId(part_id, img_id)));
    return result;
}

long int Experiment::addParticle(const std::string &part_name, int random_subset) {

    ExpParticle particle;
    particle.name = part_name;
    particle.random_subset = random_subset;

    // Push back this particle in the particles vector and its sorted index in sorted_idx
    sorted_idx.push_back(particles.size());
    particles.push_back(particle);

    // Return the current part_id in the particles vector
    return particles.size() - 1;
}

int Experiment::addImageToParticle(
    long int part_id, std::string img_name, long int ori_img_id, long int group_id, long int micrograph_id,
    int optics_group, bool unique
) {
    if (group_id >= groups.size())
        REPORT_ERROR("Experiment::addImageToParticle: group_id out of range");

    if (micrograph_id >= micrographs.size())
        REPORT_ERROR("Experiment::addImageToParticle: micrograph_id out of range");

    if (optics_group >= obsModel.numberOfOpticsGroups())
        REPORT_ERROR("Experiment::addImageToParticle: optics_group out of range");

    ExpImage img;
    img.name = img_name;
    img.id = ori_img_id;
    img.particle_id = part_id;
    img.group_id = group_id;
    img.micrograph_id = micrograph_id;
    img.optics_group = optics_group;
    if (unique)
        nr_images_per_optics_group[optics_group]++;
    img.optics_group_id = nr_images_per_optics_group[optics_group] - 1;

    if (img.optics_group_id < 0)
        REPORT_ERROR("Logic error in Experiment::addImageToParticle.");

    // Push back this particle in the particles vector
    particles[part_id].images.push_back(img);
    (micrographs[micrograph_id].image_ids).push_back(img.id);

    return particles[part_id].images.size() - 1;
}

long int Experiment::addGroup(const std::string &group_name, int _optics_group) {
    // Add new group to this Experiment
    ExpGroup group;
    group.id = groups.size(); // start counting groups at 0!
    group.optics_group = _optics_group;
    group.name = group_name;

    // Push back this micrograph
    groups.push_back(group);

    // Return the id in the micrographs vector
    return group.id;
}

long int Experiment::addMicrograph(const std::string &mic_name) {
    // Add new micrograph to this Experiment
    ExpMicrograph micrograph;
    micrograph.id = micrographs.size();
    micrograph.name = mic_name;

    // Push back this micrograph
    micrographs.push_back(micrograph);

    // Return the id in the micrographs vector
    return micrograph.id;
}

void Experiment::divideParticlesInRandomHalves(int seed, bool do_helical_refine) {
    // Only do this if the random_subset of all original_particles is zero
    bool anyzero = false, allzero = true;
    nr_particles_subset1 = 0;
    nr_particles_subset2 = 0;
    for (const auto &particle : particles) {
        if (particle.random_subset == 0) {
            anyzero = true;
        } else {
            allzero = false;
            // Keep track of how many particles there are in each subset
            if (particle.random_subset == 1) {
                nr_particles_subset1++;
            } else if (particle.random_subset == 2) {
                nr_particles_subset2++;
            } else {
                REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(particle.random_subset));
            }
        }

        if (!allzero && anyzero)
            REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: some random subset values are zero and others are not. They should all be zero, or all bigger than zero!");
    }

    if (allzero) {
        // Only randomise them if the random_subset values were not read in from the STAR file
        srand(seed);
        if (do_helical_refine) {
            std::string mic_name, img_name;
            int nr_swaps, nr_segments_subset1, nr_segments_subset2, helical_tube_id;
            std::vector<std::pair<std::string, int> > vec_mics;

            bool divide_according_to_helical_tube_id = false;
            if (MDimg.containsLabel(EMDL::PARTICLE_HELICAL_TUBE_ID))
                divide_according_to_helical_tube_id = true;

            // Count micrograph names
            std::map<std::string, int> map_mics;
            for (long int part_id = 0; part_id < particles.size(); part_id++) {
                // Get name of micrograph of the first image in this particle
                long int mic_id = particles[part_id].images[0].micrograph_id;
                mic_name = micrographs[mic_id].name;
                if (divide_according_to_helical_tube_id) {
                    long int ori_img_id = getOriginalImageId(part_id, 0);
                    helical_tube_id = MDimg.getValue<int>(EMDL::PARTICLE_HELICAL_TUBE_ID, ori_img_id);
                    if (helical_tube_id < 1)
                        REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Helical tube ID should be positive integer!");
                    mic_name += std::string("_TUBEID_");
                    mic_name += std::string(integerToString(helical_tube_id));
                }
                if (!map_mics.emplace(mic_name, 1).second)
                    map_mics[mic_name]++;
            }

            vec_mics.clear();
            for (const auto &mic : map_mics)
                vec_mics.push_back(mic);

            // NEW RANDOMISATION (better than the old one)
            nr_swaps = 0;
            for (int ptr_a = 0; ptr_a < vec_mics.size() - 1; ptr_a++) {
                std::pair<std::string, int> tmp;
                int ptr_b = round(rnd_unif(ptr_a, vec_mics.size() - 1));
                if (ptr_b <= ptr_a || ptr_b >= vec_mics.size())
                    continue;
                nr_swaps++;
                tmp = vec_mics[ptr_a];
                vec_mics[ptr_a] = vec_mics[ptr_b];
                vec_mics[ptr_b] = tmp;
            }

            // Divide micrographs into halves
            map_mics.clear();
            nr_segments_subset1 = nr_segments_subset2 = 0;
            for (auto &pair : vec_mics) {
                if (nr_segments_subset1 < nr_segments_subset2) {
                    nr_segments_subset1 += pair.second;
                    pair.second = 1;
                } else {
                    nr_segments_subset2 += pair.second;
                    pair.second = 2;
                }
                map_mics.insert(pair);
            }

            for (long int part_id = 0; part_id < particles.size(); part_id++) {
                // Get name of micrograph of the first image in this particle
                long int mic_id = particles[part_id].images[0].micrograph_id;
                mic_name = micrographs[mic_id].name;
                if (divide_according_to_helical_tube_id) {
                    long int ori_img_id = getOriginalImageId(part_id, 0);
                    helical_tube_id = MDimg.getValue<int>(EMDL::PARTICLE_HELICAL_TUBE_ID, ori_img_id);
                    if (helical_tube_id < 1)
                        REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Helical tube ID should be positive integer!");
                    mic_name += std::string("_TUBEID_") + std::string(integerToString(helical_tube_id));
                }
                particles[part_id].random_subset = map_mics[mic_name];
            }
        } else {
            for (auto &particle : particles) {
                particle.random_subset = rand() % 2 + 1;  // randomly 1 or 2
            }
        }

        // Now that random subsets have been assigned, count the number of particles in each subset and set new labels in entire MDimg
        for (long int part_id = 0; part_id < particles.size(); part_id++) {
            int random_subset = getRandomSubset(part_id);

            if (random_subset == 1) {
                nr_particles_subset1++;
            } else if (random_subset == 2) {
                nr_particles_subset2++;
            } else {
                REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
            }

            for (int img_id = 0; img_id < numberOfImagesInParticle(part_id); img_id++) {
                long int ori_img_id = getOriginalImageId(part_id, img_id);
                MDimg.setValue(EMDL::PARTICLE_RANDOM_SUBSET, random_subset, ori_img_id);
            }
        }
    }

    if (nr_particles_subset2 == 0 || nr_particles_subset1 == 0)
        REPORT_ERROR("ERROR: one of your half sets has no segments. Is rlnRandomSubset set to 1 or 2 in your particles STAR file? Or in case you're doing helical, half-sets are always per-filament, so provide at least 2 filaments.");

    std::stable_sort(sorted_idx.begin(), sorted_idx.end(), compareRandomSubsetParticles(particles));

}

void Experiment::randomiseParticlesOrder(int seed, bool do_split_random_halves, bool do_subsets) {
    //This static flag is for only randomize once
    static bool randomised = false;
    if (!randomised || do_subsets) {
        srand(seed);

        if (do_split_random_halves) {
            std::stable_sort(sorted_idx.begin(), sorted_idx.end(), compareRandomSubsetParticles(particles));

            // sanity check
            long int nr_half1 = 0, nr_half2 = 0;
            for (long int i = 0; i < particles.size(); i++) {
                const int random_subset = particles[i].random_subset;
                if (random_subset == 1) {
                    nr_half1++;
                } else if (random_subset == 2) {
                    nr_half2++;
                } else {
                    REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
                }
            }

            if (nr_half1 != nr_particles_subset1)
                REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid half1 size:" + integerToString(nr_half1) + " != " + integerToString(nr_particles_subset1));
            if (nr_half2 != nr_particles_subset2)
                REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid half2 size:" + integerToString(nr_half2) + " != " + integerToString(nr_particles_subset2));

            // Randomise the two particle lists
            std::random_shuffle(sorted_idx.begin(), sorted_idx.begin() + nr_half1);
            std::random_shuffle(sorted_idx.begin() + nr_half1, sorted_idx.end());

            // Make sure the particles are sorted on their optics_group.
            // Otherwise CudaFFT re-calculation of plans every time image size changes slows down things a lot!
            std::stable_sort(sorted_idx.begin(), sorted_idx.begin() + nr_half1, compareOpticsGroupsParticles(particles));
            std::stable_sort(sorted_idx.begin() + nr_half1, sorted_idx.end(), compareOpticsGroupsParticles(particles));

        } else {
            // Just randomise the entire vector
            std::random_shuffle(sorted_idx.begin(), sorted_idx.end());

            // Make sure the particles are sorted on their optics_group.
            // Otherwise CudaFFT re-calculation of plans every time image size changes slows down things a lot!
            std::stable_sort(sorted_idx.begin(), sorted_idx.end(), compareOpticsGroupsParticles(particles));
        }

        randomised = true;
    }
}

void Experiment::initialiseBodies(int nr_bodies) {
    if (nr_bodies < 2) return;

    this->nr_bodies = nr_bodies;
    MetaDataTable MDbody;
    MDbody.isList = false;
    const bool is_3d = MDimg.containsLabel(EMDL::ORIENT_ORIGIN_Z);
    for (auto i : MDimg) {
        MDbody.addObject();
        const RFLOAT norm = MDimg.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION, i);
        MDbody.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, 0.0, i);
        MDbody.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, 0.0, i);
        MDbody.setValue(EMDL::ORIENT_ROT,   0.0, i);
        MDbody.setValue(EMDL::ORIENT_TILT, 90.0, i);
        MDbody.setValue(EMDL::ORIENT_PSI,   0.0, i);
        MDbody.setValue(EMDL::IMAGE_NORM_CORRECTION, norm, i);
        if (is_3d)
        MDbody.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, 0.0, i);
    }
    // Now just fill all bodies with that MDbody
    MDbodies.resize(nr_bodies, MDbody);
    for (int ibody = 0; ibody < nr_bodies; ibody++) {
        MDbodies[ibody].name = "images_body_" + integerToString(ibody + 1);
    }
}

FileName Experiment::getImageNameOnScratch(long int part_id, int img_id, bool is_ctf_image) {
    const int optics_group = getOpticsGroup(part_id, img_id);
    const long int my_id = particles[part_id].images[img_id].optics_group_id;

    #ifdef DEBUG_SCRATCH
    std::cerr << "part_id = " << part_id << " img_id = " << img_id << " my_id = " << my_id << " nr_parts_on_scratch[" << optics_group << "] = " << nr_parts_on_scratch[optics_group] << std::endl;
    #define RETURN(expression) \
        std::cerr << "getImageNameOnScratch: " << particles[part_id].name << " is cached at " << fn_img << std::endl; \
        return expression;
    #else
    #define RETURN(expression) return expression;
    #endif

    if (fn_scratch.empty() || my_id >= nr_parts_on_scratch[optics_group])
        throw "!";

    if (is_3D)
        RETURN(fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + (is_ctf_image ? "_particle_ctf" : "_particle") + integerToString(my_id + 1) + ".mrc");

    // Write different optics groups into different stacks, as sizes might be different
    const auto fn_tmp = fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particles.mrcs";
    const auto fn_img = FileName::compose(my_id + 1, fn_tmp);
    RETURN(fn_img);
    #undef RETURN
}

void Experiment::setScratchDirectory(const FileName &fn_scratch, bool do_reuse_scratch, int verb) {
    this->fn_scratch = fn_scratch;
    // Make sure fn_scratch ends with a slash
    if (fn_scratch[fn_scratch.length() - 1] != '/')
        this->fn_scratch += '/';
    this->fn_scratch += "relion_volatile/";

    if (do_reuse_scratch) {
        nr_parts_on_scratch.resize(numberOfOpticsGroups(), 0);
        for (int optics_group = 0; optics_group < numberOfOpticsGroups(); optics_group++) {
            if (is_3D) {
                FileName fn_tmp = this->fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particle*.mrc";
                std::vector<FileName> fn_all;
                fn_tmp.globFiles(fn_all, true);
                nr_parts_on_scratch[optics_group] = fn_all.size();

            } else {
                FileName fn_tmp = this->fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particles.mrcs";
                if (exists(fn_tmp)) {
                    Image<RFLOAT> Itmp;
                    Itmp.read(fn_tmp, false);
                    nr_parts_on_scratch[optics_group] = Nsize(Itmp());
                }
                #ifdef DEBUG_SCRATCH
                if (verb > 0)
                    std::cerr << " optics_group= " << (optics_group + 1) << " nr_parts_on_scratch[optics_group]= " << nr_parts_on_scratch[optics_group] << std::endl;
                #endif
            }
        }
    }
}

FileName Experiment::initialiseScratchLock(const FileName &fn_out) {
    // Get a unique lockname for this run
    const int uid = rand() % 100000;
    FileName fn_uniq = fn_out;
    fn_uniq.replaceAllSubstrings("/", "_");
    fn_uniq += "_lock" + integerToString(uid);
    const FileName fn_lock = fn_scratch + fn_uniq;

    if (exists(fn_lock))
        remove(fn_lock.c_str());

    return fn_lock;
}

bool Experiment::prepareScratchDirectory(const FileName &fn_scratch, const FileName &fn_lock) {
    if (!fn_lock.empty() && exists(fn_lock)) {
        // Still measure how much free space there is
        struct statvfs vfs;
        statvfs(fn_scratch.c_str(), &vfs);
        char nodename[64] = "undefined";
        gethostname(nodename, sizeof(nodename));
        free_space_Gb = (RFLOAT) vfs.f_bsize * vfs.f_bfree / (1024 * 1024 * 1024);
        return false;
    } else {
        // Wipe the directory clean and make a new one
        deleteDataOnScratch();

        // Make the scratch directory with write permissions
        const std::string command = "install -d -m 0777 " + this->fn_scratch;
        if (system(command.c_str()))
            REPORT_ERROR("ERROR: cannot execute: " + command);

        // Touch the lock file
        if (fn_lock.empty()) {
            touch(fn_lock);
            const std::string command = "chmod 0777 " + fn_lock;
            if (system(command.c_str()))
                REPORT_ERROR("ERROR: cannot execute: " + command);
        }

        // Measure how much free space there is
        struct statvfs vfs;
        statvfs(fn_scratch.c_str(), &vfs);
        char nodename[64] = "undefined";
        gethostname(nodename, sizeof(nodename));
        free_space_Gb = (RFLOAT) vfs.f_bsize * vfs.f_bfree / (1024 * 1024 * 1024);
        std::cout << " + On host " << std::string(nodename) << ": free scratch space = " << free_space_Gb << " Gb." << std::endl;

        return true;
    }
}

void Experiment::deleteDataOnScratch() {
    // Wipe the scratch directory
    if (!fn_scratch.empty() && exists(fn_scratch)) {
        const std::string command = " rm -rf " + fn_scratch;
        if (system(command.c_str()))
            REPORT_ERROR("ERROR: cannot execute: " + command);
    }
}

void Experiment::copyParticlesToScratch(int verb, bool do_copy, bool also_do_ctf_image, RFLOAT keep_free_scratch_Gb) {
    // This function relies on prepareScratchDirectory() being called before!

    long int nr_part = MDimg.size();
    int barstep;
    if (verb > 0 && do_copy) {
        std::cout << " Copying particles to scratch directory: " << fn_scratch << std::endl;
        init_progress_bar(nr_part);
        barstep = std::max(1, (int) nr_part / 60);
    }

    long int one_part_space, used_space = 0.0;
    long int max_space = (free_space_Gb - keep_free_scratch_Gb) * 1024 * 1024 * 1024; // in bytes
    #ifdef DEBUG_SCRATCH
    std::cerr << " free_space_Gb = " << free_space_Gb << " GB, keep_free_scratch_Gb = " << keep_free_scratch_Gb << " GB.\n";
    std::cerr << " Max space RELION can use = " << max_space << " bytes" << std::endl;
    #endif
    // Loop over all particles and copy them one-by-one
    FileName fn_open_stack = "";
    long int total_nr_parts_on_scratch = 0;
    nr_parts_on_scratch.resize(numberOfOpticsGroups(), 0);

    const int check_abort_frequency = 100;

    FileName prev_img_name = "/Unlikely$filename$?*!";
    int prev_optics_group = -999;
    for (auto i : MDimg) {
        // TODO: think about MPI_Abort here....
        if (i % check_abort_frequency == 0 && pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        const FileName fn_img = MDimg.getValue<std::string>(EMDL::IMAGE_NAME, i);

        const int optics_group = [&] () {
            try {
                return MDimg.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i) - 1;
            } catch (const char *errmsg) {
                return 0;
            };
        }();

        // Get the size of the first particle
        if (nr_parts_on_scratch[optics_group] == 0) {
            const auto tmp = Image<RFLOAT>::from_filename(fn_img, false); // false means: only read the header!
            one_part_space = Xsize(tmp()) * Ysize(tmp()) * Zsize(tmp()) * sizeof(float); // MRC images are stored in floats!
            if (is_3D != (Zsize(tmp()) > 1)) REPORT_ERROR("BUG: inconsistent is_3D values!");
            // add MRC header size for subtomograms, which are stored as 1 MRC file each
            if (is_3D) {
                one_part_space += 1024;
                also_do_ctf_image = MDimg.containsLabel(EMDL::CTF_IMAGE);
                if (also_do_ctf_image) { one_part_space *= 2; }
            }
            #ifdef DEBUG_SCRATCH
            std::cerr << "one_part_space[" << optics_group << "] = " << one_part_space << std::endl;
            #endif
        }

        bool is_duplicate = prev_img_name == fn_img && prev_optics_group == optics_group;
        // Read in the particle image, and write out on scratch
        if (do_copy && !is_duplicate) {
            #ifdef DEBUG_SCRATCH
            std::cerr << "used_space = " << used_space << std::endl;
            #endif
            // Now we have the particle in memory
            // See how much space it occupies
            used_space += one_part_space;
            // If there is no more space, exit the loop over all objects to stop copying files and change filenames in MDimg
            if (used_space > max_space) {
                char nodename[64] = "undefined";
                gethostname(nodename, sizeof(nodename));
                std::cerr << " Warning: scratch space on " << std::string(nodename) << " full. "
                    "Remaining " << nr_part - total_nr_parts_on_scratch << " particles will be read from where they were." << std::endl;
                break;
            }

            if (is_3D) {
                // For subtomograms, write individual .mrc files,possibly also CTF images
                const FileName new_fn_img = fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particle" + integerToString(nr_parts_on_scratch[optics_group] + 1) + ".mrc";
                Image<RFLOAT>::from_filename(fn_img).write(new_fn_img);
                if (also_do_ctf_image) {
                    const FileName fn_ctf = MDimg.getValue<std::string>(EMDL::CTF_IMAGE, MDimg.size() - 1);
                    const FileName new_fn_ctf = fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particle_ctf" + integerToString(nr_parts_on_scratch[optics_group] + 1) + ".mrc";
                    Image<RFLOAT>::from_filename(fn_ctf).write(new_fn_ctf);
                }
            } else {
                // Only open/close new stacks, so check if this is a new stack
                long int imgno;
                FileName fn_stack;
                fn_img.decompose(imgno, fn_stack);
                fImageHandler hFile;
                if (fn_stack != fn_open_stack) {
                    // Manual closing isn't necessary: if still open, then openFile will first close the filehandler
                    // Also closing the last one isn't necessary, as destructor will do this.
                    //if (fn_open_stack != "")
                    //	hFile.closeFile();
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                Image<RFLOAT> img;
                img.readFromOpenFile(fn_img, hFile, -1, false);

                const auto fn_new = FileName::compose(nr_parts_on_scratch[optics_group] + 1, fn_scratch + "opticsgroup" + integerToString(optics_group + 1) + "_particles.mrcs");
                if (nr_parts_on_scratch[optics_group] == 0) {
                    img.write(fn_new, -1, false, WRITE_OVERWRITE);
                } else {
                    img.write(fn_new, -1, true, WRITE_APPEND);
                }

                #ifdef DEBUG_SCRATCH
                std::cerr << "Cached " << fn_img << " to " << fn_new << std::endl;
                #endif
            }
        }

        // Update the counter and progress bar
        if (!is_duplicate)
            nr_parts_on_scratch[optics_group]++;
        total_nr_parts_on_scratch++;

        prev_img_name = fn_img;
        prev_optics_group = optics_group;

        if (verb > 0 && total_nr_parts_on_scratch % barstep == 0)
            progress_bar(total_nr_parts_on_scratch);
    }

    if (verb) {
        progress_bar(nr_part);
        for (int i = 0; i < nr_parts_on_scratch.size(); i++) {
            std::cout << " For optics_group " << i + 1 << ", there are " << nr_parts_on_scratch[i] << " particles on the scratch disk." << std::endl;
        }
    }

    if (do_copy && total_nr_parts_on_scratch > 1) {
        const std::string command = " chmod -R 777 " + fn_scratch + "/";
        if (system(command.c_str()))
            REPORT_ERROR("ERROR in executing: " + command);
    }
}

// Read from file
void Experiment::read(
    FileName fn_exp, bool do_ignore_particle_name, bool do_ignore_group_name, bool do_preread_images,
    bool need_tiltpsipriors_for_helical_refine, int verb
) {

    // #define DEBUG_READ
    #ifdef DEBUG_READ
    std::cerr << "Entering Experiment::read" << std::endl;
    Timer timer;
    int tall   = timer.setNew("ALL");
    int tread  = timer.setNew("read");
    int tsort  = timer.setNew("sort");
    int tfill  = timer.setNew("fill");
    int tgroup = timer.setNew("find group");
    int tdef   = timer.setNew("set defaults");
    int tend   = timer.setNew("ending");
    char c;
    timer.tic(tall);
    timer.tic(tread);
    #endif

    // Only open stacks once and then read multiple images
    fImageHandler hFile;
    long int dump;
    FileName fn_stack, fn_open_stack = "";

    // Initialize by emptying everything
    clear();
    long int group_id = 0, mic_id = 0, part_id = 0;

    if (!fn_exp.isStarFile()) {
        // Read images from stack. Ignore all metadata, just use filenames

        // Add a single Micrograph
        group_id = addGroup("group", 0);
        mic_id = addMicrograph("micrograph");

        // Check that a MRC stack ends in .mrcs, not .mrc (which will be read as a MRC 3D map!)
        if (fn_exp.contains(".mrc") && !fn_exp.contains(".mrcs"))
            REPORT_ERROR("Experiment::read: ERROR: MRC stacks of 2D images should be have extension .mrcs, not .mrc!");

        // Read in header-only information to get the Nsize of the stack
        auto img = Image<RFLOAT>::from_filename(fn_exp, false);  // false means skip data, only read header

        // allocate 1 block of memory
        particles.reserve(Nsize(img()));
        nr_images_per_optics_group.resize(1, 0);

        for (long int n = 0; n < Nsize(img()); n++) {
            const auto fn_img = FileName::compose(n + 1, fn_exp);  // fn_img = integerToString(n) + "@" + fn_exp;
            // Add the particle to my_area = 0
            part_id = addParticle(fn_img, 0);
            // Just add a single image per particle
            addImageToParticle(part_id, fn_img, n, 0, 0, 0, true);

            MDimg.addObject();

            if (do_preread_images) {
                fn_img.decompose(dump, fn_stack);
                if (fn_stack != fn_open_stack) {
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                Image<float> img;
                img.readFromOpenFile(fn_img, hFile, -1, false);
                particles[part_id].images[0].img = img().setXmippOrigin();
            }

            // Set the filename and other metadata parameters
            MDimg.setValue(EMDL::IMAGE_NAME, fn_img, part_id);
            MDimg.setValue(EMDL::IMAGE_OPTICS_GROUP, 1, part_id);
        }
    } else {
        // MDimg and MDopt have to be read at the same time, so that the optics groups can be
        // renamed in case they are non-contiguous or not sorted
        ObservationModel::loadSafely(fn_exp, obsModel, MDimg, "particles", verb);
        nr_images_per_optics_group.resize(obsModel.numberOfOpticsGroups(), 0);

        #ifdef DEBUG_READ
        std::cerr << "Done reading MDimg" << std::endl;
        timer.toc(tread);
        timer.tic(tsort);
        //std::cerr << "Press any key to continue..." << std::endl;
        //std::cin >> c;
        #endif

        // Sort input particles on micrographname
        bool is_mic_a_movie = false, star_contains_micname;
        star_contains_micname = MDimg.containsLabel(EMDL::MICROGRAPH_NAME);
        if (star_contains_micname) {
            // See if the micrograph names contain an "@", i.e. whether they are movies and we are inside polishing or so.
            const FileName fn_mic = MDimg.getValue<std::string>(EMDL::MICROGRAPH_NAME, MDimg.size() - 1);
            if (fn_mic.contains("@")) {
                is_mic_a_movie = true;
                MDimg.newSort<MD::CompareStringsAfterAtAt>(EMDL::MICROGRAPH_NAME);
            } else {
                is_mic_a_movie = false;
                MDimg.newSort<MD::CompareStringsAt>(EMDL::MICROGRAPH_NAME);
            }

            if (do_ignore_group_name)
                group_id = addGroup("group", 0);
        } else {
            // If there is no EMDL::MICROGRAPH_NAME, then just use a single group and micrograph
            group_id = addGroup("group", 0);
            mic_id = addMicrograph("micrograph");
        }
        #ifdef DEBUG_READ
        std::cerr << "Done sorting MDimg" << std::endl;
        std::cerr << " MDimg.size()= " << MDimg.size() << std::endl;
        timer.toc(tsort);
        timer.tic(tfill);
        long nr_read = 0;
        #endif
        // allocate 1 block of memory
        particles.reserve(MDimg.size());

        // Now Loop over all objects in the metadata file and fill the logical tree of the experiment
        long int last_part_id = -1;

        FileName prev_img_name = "/Unlikely$filename$?*!";
        int prev_optics_group = -999;
        //for (auto _ : MDimg)
        // Loop over all objects in MDimg (ori_part_id)
        for (long int ori_img_id = 0; ori_img_id < MDimg.size(); ori_img_id++) {
            // Get the optics group of this particle
            int optics_group = obsModel.getOpticsGroup(MDimg, ori_img_id);

            // Add new micrographs or get mic_id for existing micrograph
            FileName mic_name = ""; // Filename instead of string because will decompose below
            if (star_contains_micname) {
                long int idx = micrographs.size();
                std::string last_mic_name = idx > 0 ? micrographs[idx - 1].name : "";

                mic_name = MDimg.getValue<std::string>(EMDL::MICROGRAPH_NAME, ori_img_id);

                // All frames of a movie belong to the same micrograph
                if (is_mic_a_movie)
                    mic_name = mic_name.substr(mic_name.find("@") + 1);

                mic_id = -1;
                if (last_mic_name == mic_name) {
                    // This particle belongs to the previous micrograph
                    mic_id = micrographs[idx - 1].id;
                } else {
                    // A new micrograph
                    last_part_id = particles.size();
                }

                // Make a new micrograph
                if (mic_id < 0)
                    mic_id = addMicrograph(mic_name);

                #ifdef DEBUG_READ
                timer.tic(tgroup);
                #endif

                // For example in particle_polishing the groups are not needed...
                if (!do_ignore_group_name) {
                    std::string group_name;
                    // Check whether there is a group label, if not use a group for each micrograph
                    if (MDimg.containsLabel(EMDL::MLMODEL_GROUP_NAME)) {
                        group_name = MDimg.getValue<std::string>(EMDL::MLMODEL_GROUP_NAME, ori_img_id);
                    } else {
                        FileName fn_pre, fn_jobnr, fn_post;
                        decomposePipelineFileName(mic_name, fn_pre, fn_jobnr, fn_post);
                        group_name = fn_post;
                    }

                    // If this group did not exist yet, add it to the experiment
                    group_id = -1;
                    for (long int i = groups.size() - 1; i >= 0; i--) {
                    // search backwards to find match faster
                        if (groups[i].name == group_name) {
                            group_id = groups[i].id;
                            break;
                        }
                    }
                    if (group_id < 0) {
                        group_id = addGroup(group_name, optics_group);
                    }
                }

                #ifdef DEBUG_READ
                timer.toc(tgroup);
                #endif

            } else {
                // All images belong to the same micrograph and group
                mic_id = 0;
                group_id = 0;
            }

            // If there is an EMDL::PARTICLE_RANDOM_SUBSET entry in the input STAR-file, then set the random_subset, otherwise use default (0)
            int my_random_subset;
            try {
                my_random_subset = MDimg.getValue<int>(EMDL::PARTICLE_RANDOM_SUBSET, ori_img_id);
            } catch (const char *errmsg) {
                my_random_subset = 0;
            }

            // Add this image to an existing particle, or create a new particle
            std::string part_name = MDimg.getValue<std::string>(
                MDimg.containsLabel(EMDL::PARTICLE_NAME) ? EMDL::PARTICLE_NAME : EMDL::IMAGE_NAME,
                ori_img_id
            );

            long int part_id = -1;
            if (MDimg.containsLabel(EMDL::PARTICLE_NAME) && !do_ignore_particle_name) {
                // Only search ori_particles for the last (original) micrograph
                for (long int i = last_part_id; i < particles.size(); i++) {
                    if (particles[i].name == part_name) {
                        part_id = i;
                        break;
                    }
                }
            }

            // If no particles with this name was found,
            // or if no EMDL::PARTICLE_NAME in the input file, or if do_ignore_original_particle_name
            // then add a new particle
            if (part_id < 0) {
                part_id = addParticle(part_name, my_random_subset);
            }

            // Create a new image in this particle
            FileName img_name = MDimg.getValue<std::string>(EMDL::IMAGE_NAME, ori_img_id);

            bool do_cache = prev_img_name != img_name || prev_optics_group != optics_group;
            #ifdef DEBUG_SCRATCH
            std::cerr << "prev_img_name = " << prev_img_name << " img_name = " << img_name << " prev_optics_group = " << prev_optics_group << " optics_group = " << optics_group << " do_cache = " << do_cache << std::endl;
            #endif
            prev_img_name = img_name;
            prev_optics_group = optics_group;

            int img_id = addImageToParticle(part_id, img_name, ori_img_id, group_id, mic_id, optics_group, do_cache);

            // The group number is only set upon reading: it is not read from the STAR file itself,
            // there the only thing that matters is the order of the micrograph_names
            // Write igroup+1, to start numbering at one instead of at zero
            MDimg.setValue(EMDL::MLMODEL_GROUP_NO, group_id + 1, ori_img_id);

            #ifdef DEBUG_READ
            timer.tic(tori);
            #endif

            if (do_preread_images) {
                img_name.decompose(dump, fn_stack);
                if (fn_stack != fn_open_stack) {
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                Image<float> img;
                img.readFromOpenFile(img_name, hFile, -1, false);
                particles[part_id].images[img_id].img = img().setXmippOrigin();
            }

            #ifdef DEBUG_READ
            timer.toc(tori);
            #endif

            #ifdef DEBUG_READ
            nr_read++;
            #endif
        }

        #ifdef DEBUG_READ
        timer.toc(tfill);
        timer.tic(tdef);
        std::cerr << " nr_read= " << nr_read << " particles.size()= " << particles.size() << " micrographs.size()= " << micrographs.size() << " groups.size()= " << groups.size() << std::endl;
        #endif

        // Check for the presence of multiple bodies (for multi-body refinement)
        nr_bodies = 0;
        while (true) {
            const std::string tablename = "images_body_" + integerToString(nr_bodies + 1);
            MetaDataTable MDimgin;
            if (MDimgin.read(fn_exp, tablename) <= 0) break;
            nr_bodies++;
            MDbodies.push_back(MDimgin);
        }
        // Even if we don't do multi-body refinement, then nr_bodies is still 1
        nr_bodies = std::max(nr_bodies, 1);
    }

    #ifdef DEBUG_READ
    std::cerr << "Done filling MDimg" << std::endl;
    //std::cerr << "Press any key to continue..." << std::endl;
    //std::cin >> c;
    #endif

    // Make sure some things are always set in the MDimg
    bool have_rot  = MDimg.containsLabel(EMDL::ORIENT_ROT);
    bool have_tilt = MDimg.containsLabel(EMDL::ORIENT_TILT);
    bool have_psi  = MDimg.containsLabel(EMDL::ORIENT_PSI);
    bool have_xoff = MDimg.containsLabel(EMDL::ORIENT_ORIGIN_X_ANGSTROM);
    bool have_yoff = MDimg.containsLabel(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);
    bool have_zoff = MDimg.containsLabel(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);
    bool have_zcoord = MDimg.containsLabel(EMDL::IMAGE_COORD_Z);
    bool have_clas = MDimg.containsLabel(EMDL::PARTICLE_CLASS);
    bool have_norm = MDimg.containsLabel(EMDL::IMAGE_NORM_CORRECTION);

    // 20 Jan 2016 - Helical reconstruction
    bool have_tilt_prior = MDimg.containsLabel(EMDL::ORIENT_TILT_PRIOR);
    bool have_psi_prior = MDimg.containsLabel(EMDL::ORIENT_PSI_PRIOR);
    bool have_tiltpsi = have_tilt && have_psi;
    bool have_tiltpsi_prior = have_tilt_prior && have_psi_prior;
    if (need_tiltpsipriors_for_helical_refine && !have_tiltpsi_prior && !have_tiltpsi) {
        REPORT_ERROR("exp_model.cpp: Experiment::read(): Tilt and psi priors of helical segments are missing!");
    }
    for (auto i : MDimg) {
        if (!have_rot)
            MDimg.setValue(EMDL::ORIENT_ROT, 0.0, i);
        if (!have_tilt)
            MDimg.setValue(EMDL::ORIENT_TILT, 0.0, i);
        if (!have_psi)
            MDimg.setValue(EMDL::ORIENT_PSI, 0.0, i);
        if (!have_xoff)
            MDimg.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, 0.0, i);
        if (!have_yoff)
            MDimg.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, 0.0, i);
        if (!have_zoff && have_zcoord)
            MDimg.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, 0.0, i);
        if (!have_clas)
            MDimg.setValue(EMDL::PARTICLE_CLASS, 0, i);
        if (!have_norm)
            MDimg.setValue(EMDL::IMAGE_NORM_CORRECTION, 1.0, i);
        if (need_tiltpsipriors_for_helical_refine && have_tiltpsi_prior) {
            // If doing 3D helical reconstruction and PRIORs exist
            RFLOAT tilt = 0.0, psi = 0.0;
            if (have_tiltpsi)
                tilt = MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT, i);
            // If ANGLEs do not exist or they are all set to 0 (from a Class2D job),
            // copy values of PRIORs to ANGLEs
            if (!have_tiltpsi || abs(tilt) < 0.001) {
                tilt = MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT_PRIOR, i);
                psi  = MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI_PRIOR, i);
                MDimg.setValue(EMDL::ORIENT_TILT, tilt, i);
                MDimg.setValue(EMDL::ORIENT_PSI,  psi,  i);
            }
        }
    }

    // Set is_3D from MDopt
    is_3D = obsModel.opticsMdt.getValue<int>(EMDL::IMAGE_DIMENSIONALITY, 0) == 3;

    #ifdef DEBUG_READ
    timer.toc(tdef);
    std::cerr << "Done setting defaults MDimg" << std::endl;
    timer.toc(tall);
    timer.printTimes(false);
    //std::cerr << "Writing out debug_data.star" << std::endl;
    //write("debug");
    //exit(0);
    #endif
}

// Write to file
void Experiment::write(FileName fn_root) {
    const FileName fn_tmp = fn_root + "_data.star";
    std::ofstream fh ((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR((std::string) "Experiment::write: Cannot write file: " + fn_tmp);

    obsModel.opticsMdt.name = "optics";
    obsModel.opticsMdt.write(fh);

    // Always write MDimg
    MDimg.name = "particles";
    MDimg.write(fh);

    if (nr_bodies > 1) {
        for (int ibody = 0; ibody < nr_bodies; ibody++) {
            MDbodies[ibody].write(fh);
        }
    }

}
