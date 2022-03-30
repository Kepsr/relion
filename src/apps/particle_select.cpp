#include <src/jaz/stack_helper.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/metadata_table.h>

using namespace gravis;


int main(int argc, char *argv[]) {
    IOParser parser;

    parser.setCommandLine(argc, argv);
    int gen_section = parser.addSection("General options");

    std::string sourceFn = parser.getOption("--i", "Input STAR file containing the source particles");
    std::string refFn = parser.getOption("--i_ref", "Input STAR file containing reference particles");

    const bool copyAngles = parser.checkOption("--angles", "Copy particle viewing angles from reference");
    const bool copyOffsets = parser.checkOption("--offsets", "Copy particle offsets from reference");

    std::string outFn = parser.getOption("--o", "Output path", "selected.star");

    if (parser.checkForErrors()) {
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
    }

    MetaDataTable sourceAll, refAll;

    sourceAll.read(sourceFn);
    refAll.read(refFn);

    std::vector<MetaDataTable> sourceByMic = StackHelper::splitByMicrographName(sourceAll);
    std::vector<MetaDataTable> refByMic = StackHelper::splitByMicrographName(refAll);

    std::map<std::string, MetaDataTable*> micToSource;

    for (int m = 0; m < sourceByMic.size(); m++) {
        std::string micName = sourceByMic[m].getValue(EMDL::MICROGRAPH_NAME, 0);
        micToSource[micName] = &sourceByMic[m];
    }

    MetaDataTable out;

    for (int m = 0; m < refByMic.size(); m++) {
        std::string micName = refByMic[m].getValue(EMDL::MICROGRAPH_NAME, 0);

        if (micToSource.find(micName) == micToSource.end()) {
            std::cerr << "Warning: " << micName << " not found.\n";
            continue;
        }

        MetaDataTable* src = micToSource[micName];

        const int pcRef = refByMic[m].numberOfObjects();
        const int pcSrc = src->numberOfObjects();

        std::vector<d2Vector> posSrc(pcSrc);

        for (int p = 0; p < pcSrc; p++) {
            posSrc[p].x = src->getValue(EMDL::IMAGE_COORD_X, p);
            posSrc[p].y = src->getValue(EMDL::IMAGE_COORD_Y, p);
        }

        std::vector<d2Vector> posRef(pcRef);

        for (int p = 0; p < pcRef; p++) {
            posRef[p].x = refByMic[m].getValue(EMDL::IMAGE_COORD_X, p);
            posRef[p].y = refByMic[m].getValue(EMDL::IMAGE_COORD_Y, p);
        }

        int missing = 0, multiple = 0;

        for (int p = 0; p < pcRef; p++) {
            int qBest = -1;

            for (int q = 0; q < pcSrc; q++) {
                double dist = (posRef[p] - posSrc[q]).length();

                if (dist < 1.0) {
                    qBest = qBest == -1 ? q : -2;
                }
            }

            if (qBest >= 0) {
                out.addObject(src->getObject(qBest));
                const int qNew = out.numberOfObjects() - 1;

                int randSubsetSrc = src->       getValue(EMDL::PARTICLE_RANDOM_SUBSET, qBest);
                int randSubsetRef = refByMic[m].getValue(EMDL::PARTICLE_RANDOM_SUBSET, p);

                if (randSubsetSrc != randSubsetRef) {
                    if (copyAngles && copyOffsets) {
                        out.setValue(EMDL::PARTICLE_RANDOM_SUBSET, randSubsetRef, qNew);
                    } else if (copyAngles != copyOffsets) {
                        REPORT_ERROR_STR("Unable to copy only angles or only offsets, since the "
                                         << "particles belong to different random subsets.\n");
                    }
                }

                if (copyAngles) {
                    double rot  = refByMic[m].getValue(EMDL::ORIENT_ROT,  p);
                    double tilt = refByMic[m].getValue(EMDL::ORIENT_TILT, p);
                    double psi  = refByMic[m].getValue(EMDL::ORIENT_PSI,  p);
                    out.setValue(EMDL::ORIENT_ROT,  rot, qNew);
                    out.setValue(EMDL::ORIENT_TILT, tilt, qNew);
                    out.setValue(EMDL::ORIENT_PSI,  psi, qNew);
                }

                if (copyOffsets) {
                    double xoff = refByMic[m].getValue(EMDL::ORIENT_ORIGIN_X, p);
                    double yoff = refByMic[m].getValue(EMDL::ORIENT_ORIGIN_Y, p);
                    out.setValue(EMDL::ORIENT_ORIGIN_X, xoff, qNew);
                    out.setValue(EMDL::ORIENT_ORIGIN_Y, yoff, qNew);
                }
            } else if (qBest == -1) {
                missing++;
            } else {
                // -2
                multiple++;
            }
        }

        if (missing > 0) {
            std::cerr << "    Warning: " << missing << " of " << pcRef
                      << " particles missing from micrograph " << m << "\n";
        }

        if (multiple > 0) {
            std::cerr << "    Warning: " << multiple << " out of " << pcRef
                      << " particles found multiple times in micrograph " << m << "\n"
                      << "    (all will be ignored)\n";
        }
    }

    out.write(outFn);

    return RELION_EXIT_SUCCESS;
}
