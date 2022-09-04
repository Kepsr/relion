/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres" "Jasenko Zivanov"
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

#include <src/args.h>
#include <src/metadata_table.h>
#include <src/micrograph_model.h>
#include <src/jaz/io/star_converter.h>

class star_converter {

    public:

    FileName fn_in, fn_out;
    IOParser parser;
    RFLOAT Cs, Q0;

    void usage() { parser.writeUsage(std::cerr); }

    void read(int argc, char **argv) {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("Options");
        fn_in = parser.getOption("--i", "Input STAR file to be converted", "None");
        fn_out = parser.getOption("--o", "Output STAR file to be written", "None");
        Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration (mm)", "-1"));
        Q0 = textToFloat(parser.getOption("--Q0", "Amplitude contrast", "-1"));

        if (fn_in == "None" || fn_out == "None") {
            usage();
            REPORT_ERROR("Please specify input and output file names");
        }
    }

    void run() {
        MetaDataTable mdt;
        MetaDataTable mdtOut, optOut;
        mdt.read(fn_in);
        const bool isMotionCorrSTAR = mdt.containsLabel(EMDL::MICROGRAPH_METADATA_NAME);
        StarConverter::convert_3p0_particlesTo_3p1(mdt, mdtOut, optOut, "", false); // don't die

        if (mdt.containsLabel(EMDL::IMAGE_NAME)) {
            std::cout << "The input is a particle STAR file" << std::endl;
            mdtOut.name = "particles";
        } else if (isMotionCorrSTAR) {
            std::cout << "The input is a STAR file from a MotionCorr job." << std::endl;
            std::cout << "The (binned) pixel size and the voltage are taken from the first metadata STAR file." << std::endl;
            FileName fn_meta;
            try {
                fn_meta = mdtOut.getValue<std::string>(EMDL::MICROGRAPH_METADATA_NAME, 0);
            } catch (const char *errmsg) {
                REPORT_ERROR("Failed to find the metadata STAR file");
            }

            Micrograph mic(fn_meta);
            const long int i = optOut.size() - 1;

            std::cout << "- voltage: " << mic.voltage << std::endl;
            optOut.setValue(EMDL::CTF_VOLTAGE, mic.voltage, i);

            std::cout << "- unbinned pixel size: " << mic.angpix << std::endl;
            std::cout << "- binning factor: " << mic.getBinningFactor() << std::endl;
            const RFLOAT angpix = mic.angpix * mic.getBinningFactor();
            std::cout << "- binned pixel size: " << angpix << std::endl;
            optOut.setValue(EMDL::MICROGRAPH_PIXEL_SIZE, angpix, i);

            std::cout << "\nThe other microscope parameters must be specified in the command line." << std::endl;
            if (Cs < 0)
                REPORT_ERROR("Please specify the spherical aberration (mm) in the --Cs option.");
            std::cout << "- spherical aberration: " << Cs << std::endl;
            optOut.setValue(EMDL::CTF_CS, Cs, i);
            if (Q0 < 0)
                REPORT_ERROR("Please specify the amplitude contrast in the --Q0 option");
            std::cout << "- amplitude contrast: " << Q0 << std::endl;
            optOut.setValue(EMDL::CTF_Q0, Q0, i);

            std::cout << "\nAll necessary information is ready." << std::endl;

            mdtOut.name = "micrographs";
        } else {
            std::cout << "The input is a micrograph STAR file with CTF information." << std::endl;
            mdtOut.name = "micrographs";
        }
    
        std::ofstream of(fn_out);

        optOut.write(of);
        mdtOut.write(of);
        of.close();

        std::cout << "\nWritten " << fn_out << std::endl;
        std::cout << "Please carefully examine the optics group table at the beginning of the output to make sure the information is correct." << std::endl;
    }
};

int main(int argc, char *argv[]) {
    star_converter app;

    try {
        app.read(argc, argv);
        app.run();
    } catch (RelionError XE) {
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    return RELION_EXIT_SUCCESS;
}
