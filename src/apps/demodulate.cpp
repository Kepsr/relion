#include "src/jaz/vtk_helper.h"
#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/metadata_table.h"
#include "src/args.h"
#include "src/jaz/parallel_ft.h"

int main(int argc, char *argv[]) {
    std::string particlesFn, outPath;
    MetaDataTable particlesMdt;
    int nr_omp_threads;
    bool r31;

    IOParser parser;

    try {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        particlesFn = parser.getOption("--i", "Input STAR file with a list of particles");
        nr_omp_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
        r31 = parser.checkOption("--r31", "Write output in Relion-3.1 format");
        outPath = parser.getOption("--out", "Output path");

        if (parser.checkForErrors()) return RELION_EXIT_FAILURE;
    } catch (RelionError XE) {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    if (outPath[outPath.length()-1] != '/') {
        outPath += "/";
    }

    std::string command = " mkdir -p " + outPath;
    int res = system(command.c_str());

    ObservationModel obsModel;

    ObservationModel::loadSafely(particlesFn, obsModel, particlesMdt);

    particlesMdt.read(particlesFn);

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    std::vector<MetaDataTable> mdts = StackHelper::splitByStack(particlesMdt);

    const int mc = mdts.size();

    for (MetaDataTable &mdt : mdts) {

        std::vector<Image<Complex>> obs = StackHelper::loadStackFS(mdt, "", nr_omp_threads, false);

        std::string fullName = mdt.getValue<std::string>(EMDL::IMAGE_NAME, 0);
        std::string name = fullName.substr(fullName.find("@") + 1);

        const int pc = mdt.size();
        for (int p = 0; p < pc; p++) {
            int opticsGroup = mdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;
            if (p) {
                obsModel.modulatePhase(opticsGroup, obs[p].data);
            } else {
                obsModel.demodulatePhase(opticsGroup, obs[p].data);
            }
        }

        std::vector<Image<RFLOAT>> demodulated = StackHelper::inverseFourierTransform(obs);
        Image<RFLOAT> out = StackHelper::toSingleImage(demodulated);

        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(name, fn_pre, fn_jobnr, fn_post);

        std::string outFn = outPath + fn_post;

        if (outFn.find_last_of("/") != std::string::npos) {
            system((" mkdir -p " + outFn.substr(0, outFn.find_last_of("/"))).c_str());
        }

        for (int p = 0; p < pc; p++) {
            mdt.setValue(EMDL::IMAGE_NAME, std::to_string(p + 1) + "@" + outFn, p);
        }

        out.write(outFn);
    }

    MetaDataTable mdt1;
    for (const MetaDataTable &mdt: mdts) {
        mdt1.append(mdt);
    }

    if (!r31) {
        const int tpc = mdt1.size();

        std::vector<EMDL::EMDLabel> allOpticsLabels_double(0);

        allOpticsLabels_double.push_back(EMDL::CTF_Q0);
        allOpticsLabels_double.push_back(EMDL::CTF_CS);
        allOpticsLabels_double.push_back(EMDL::CTF_VOLTAGE);
        allOpticsLabels_double.push_back(EMDL::CTF_DETECTOR_PIXEL_SIZE);
        allOpticsLabels_double.push_back(EMDL::CTF_MAGNIFICATION);

        for (int l = 0; l < allOpticsLabels_double.size(); l++) {
            EMDL::EMDLabel lab = allOpticsLabels_double[l];

            mdt1.addLabel(lab);

            for (int p = 0; p < tpc; p++) {
                int opticsGroup = mdt1.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;
                double v = obsModel.opticsMdt.getValue<double>(lab, opticsGroup);
                mdt1.setValue(lab, v, p);
            }
        }

        obsModel.opticsMdt.deactivateLabel(EMDL::IMAGE_OPTICS_GROUP);

        mdt1.version = 30000;
    } else {
        obsModel.opticsMdt.deactivateLabel(EMDL::IMAGE_BEAMTILT_X);
        obsModel.opticsMdt.deactivateLabel(EMDL::IMAGE_BEAMTILT_Y);
        obsModel.opticsMdt.deactivateLabel(EMDL::IMAGE_ODD_ZERNIKE_COEFFS);

        obsModel.opticsMdt.write(outPath + "demodulated_particles_optics.star");
    }

    mdt1.write(outPath + "demodulated_particles.star");

    std::cout << "output written into " << (outPath + "demodulated_particles.star") << "\n";

    return RELION_EXIT_SUCCESS;
}
