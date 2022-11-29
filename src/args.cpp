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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/args.h"
#include "src/gcc_version.h"
#include "src/matrix1d.h"
#include "src/colour.h"
#include <algorithm>


int find_option_or_param(int argc, char **argv, const std::string &option) {
    int found_at = -1;
    for (int i = 0; i < argc; i++) {
        // std::cout << i << " " << found_at << " " << argv[i] << " looking for " << option << std::endl;
        if (strcmp(option.c_str(), argv[i]) == 0) {
            if (found_at != -1)
                std::cerr << "WARNING: Command-line option " << option << " was specified more than once. The value specified last will be used." << std::endl;
            found_at = i;
        }
    }
    return found_at;
}

// Get parameters from the command line ====================================
std::string getParameter(
    int argc, char **argv, const std::string &param, const std::string &option
) {
    int found_at = find_option_or_param(argc, argv, param);
    if (0 < found_at && found_at < argc - 1) {
        return argv[found_at + 1];
    } else {
        if (option == "NULL")
            REPORT_ERROR((std::string) "Argument " + param + " not found or invalid argument");

        return option;
    }
}

// Check if a parameter was included the command line =============
bool checkParameter(int argc, char **argv, const std::string &param) {
    for (int i = 0; i < argc; ++i) {
        if (strcmp(param.c_str(), argv[i]) == 0)
            return true;
    }
    return false;
}

IOParser::IOParser() {
    clear();
}

IOParser::IOParser(const IOParser &other):
    options(other.options),
    usages(other.usages),
    optionals(other.optionals),
    defaultvalues(other.defaultvalues),
    argc(other.argc),
    argv(other.argv),
    error_messages(other.error_messages),
    warning_messages(other.warning_messages),
    current_section(other.current_section),
    section_names(other.section_names),
    section_numbers(other.section_numbers) {}

IOParser& IOParser::operator = (IOParser other) {
    swap(other);
    return *this;
}

void IOParser::swap(IOParser &other) {
    std::swap(options,          other.options);
    std::swap(usages,           other.usages);
    std::swap(optionals,        other.optionals);
    std::swap(defaultvalues,    other.defaultvalues);
    std::swap(argc,             other.argc);
    std::swap(argv,             other.argv);
    std::swap(error_messages,   other.error_messages);
    std::swap(warning_messages, other.warning_messages);
    std::swap(current_section,  other.current_section);
    std::swap(section_names,    other.section_names);
    std::swap(section_numbers,  other.section_numbers);
}

void IOParser::clear() {
    argc = 0;
    argv = nullptr;
    options.clear();
    usages.clear();
    optionals.clear();
    defaultvalues.clear();
    error_messages.clear();
    warning_messages.clear();
    section_names.clear();
    section_numbers.clear();
    current_section = 0;
}

void IOParser::setCommandLine(int argc, char **argv) {
    this->argc = argc;
    this->argv = argv;

    // Print version of software and exit
    if (checkParameter(argc, argv, "--version")) {
        PRINT_VERSION_INFO();
        exit(0);
    }
    // Dirty hack to get pipeline control for all programs...
    pipeline_control_outputname = checkParameter(argc, argv, "--pipeline_control") ?
        getParameter(argc, argv, "--pipeline_control") : "";
}

void IOParser::addOption(
    const std::string &option, const std::string &usage, const std::string &defaultvalue,
    bool hidden
) {
    if (hidden) {
        hiddenOptions.push_back(option);
        return;
    }
    if (section_names.empty())
        REPORT_ERROR("IOParser::addOption: ERROR First add a section to the parser, then the options!");
    options.push_back(option);
    usages.push_back(usage);
    section_numbers.push_back(current_section);
    const bool isNULL = defaultvalue == "NULL";
    optionals.push_back(!isNULL);
    defaultvalues.push_back(isNULL ? " " : defaultvalue);
}

int IOParser::addSection(std::string name) {
    current_section = section_names.size();
    section_names.push_back(name);
    return current_section;
}

/** Set the current section to this number */
void IOParser::setSection(int number) {
    current_section = number;
}

bool IOParser::optionExists(std::string option) {
    return std::find(options.begin(), options.end(), option) != options.end() ||
           std::find(hiddenOptions.begin(), hiddenOptions.end(), option) != hiddenOptions.end();
}

std::string IOParser::getOption(
    const std::string &option, const std::string &usage, const std::string &defaultvalue, bool hidden
) {
    // If this option did not exist yet, add it to the list
    if (!optionExists(option))
        addOption(option, usage, defaultvalue, hidden);

    const int found_at = find_option_or_param(argc, argv, option);
    if (0 < found_at && found_at < argc - 1) {
        return argv[found_at + 1];
    } else {
        if (defaultvalue == "NULL") {
            error_messages.push_back((std::string) "ERROR: Argument " + option + " not found or invalid argument");
            return "";
        }
        return defaultvalue;
    }
}

// Checks if a boolean parameter was included the command line =============
bool IOParser::checkOption(
    const std::string &option, const std::string &usage, const std::string &defaultvalue, bool hidden
) {
    // If this option did not exist yet, add it to the list
    if (!optionExists(option))
        addOption(option, usage, defaultvalue, hidden);

    return checkParameter(argc, argv, option);
}

void IOParser::writeCommandLine(std::ostream &out) {
    for (int i = 1; i < argc; i++)
        out << argv[i] << " ";
    out << std::endl;
}

bool IOParser::checkForErrors(int verb) {

    if (checkParameter(argc, argv, "--version")) {
        std::cout << "RELION version " << g_RELION_VERSION << std::endl;
        exit(0);
    }
    if (
        argc == 1
        || argc == 2 && checkParameter(argc, argv, "--continue")
        || checkParameter(argc, argv, "--help")
        || checkParameter(argc, argv, "-h")
    ) {
        writeUsage(std::cout);
        exit(0);
    }

    // First check the command line for unknown arguments
    checkForUnknownArguments();

    // First print warning messages
    if (!warning_messages.empty())
        if (verb > 0) {
            std::cerr << "The following warnings were encountered upon command-line parsing: " << std::endl;
            for (const std::string &wrnmsg : warning_messages)
                std::cerr << wrnmsg << std::endl;
        }

    // Then check for error messages
    if (!error_messages.empty()) {
        if (verb > 0) {
            std::cerr << "The following errors were encountered upon command-line parsing: " << std::endl;
            for (const std::string &errmsg : error_messages)
                std::cerr << errmsg << std::endl;
        }
        return true;
    }
    return false;
}

bool is_ok(IOParser parser, char **argv, int i) {
    // Valid options should start with "--"
    if (strncmp("--", argv[i], 2) == 0)
        return parser.optionExists(std::string(argv[i])) || strncmp("--pipeline_control", argv[i], 18) == 0;
    if (strncmp("-", argv[i], 1) == 0) {
        // If argv[i] starts with one "-", it must be a number and argv[i - 1] must be a valid option.
        float testval;
        return sscanf(argv[i], "%f", &testval) && parser.optionExists(argv[i - 1]);
    }
        return true;
}

void IOParser::checkForUnknownArguments() {
    for (int i = 1; i < argc; i++)
        if (!is_ok(*this, argv, i))
            warning_messages.push_back(
                (std::string) "WARNING: Option " + argv[i] + "\tis not a valid RELION argument"
            );
}

void IOParser::writeUsageOneLine(int i, std::ostream &out) {
    std::string aux = "  " + options[i];
    if (optionals[i])
        aux += " (" + defaultvalues[i] + ")";

    out << std::setw(35) << aux << " : " << usages[i] << std::endl;
}

void IOParser::writeUsageOneSection(int section, std::ostream &out) {
    // First write all compulsory options
    //out << "+++ Compulsory:" << std::endl;
    for (int i = 0; i < options.size(); i++)
        if (!optionals[i] && section_numbers[i] == section)
            writeUsageOneLine(i, out);

    // Then write optional ones
    //out << "+++ Optional (defaults between parentheses):" << std::endl;
    for (int i = 0; i < options.size(); i++)
        if (optionals[i] && section_numbers[i] == section)
            writeUsageOneLine(i, out);
}

void IOParser::writeUsage(std::ostream &out) {
    out << "+++ RELION: command line arguments (with defaults for optional ones between parantheses) +++"<<std::endl;

    for (int i = 0; i < section_names.size(); i++) {
        out << "====== " << section_names[i] << " ===== " << std::endl;
        writeUsageOneSection(i, out);
    }
    out << std::setw(35) << "--version" << " : Print RELION version and exit" << std::endl;
}

ColourScheme IOParser::getColourScheme() {
    return
    checkOption(
        "--colour_fire",
        "Show images in black-grey-white-red colour scheme (highlight high signal)?"
    ) ? black_grey_red :
    checkOption(
        "--colour_ice",
        "Show images in blue-black-grey-white colour scheme (highlight low signal)?"
    ) ? blue_grey_white :
    checkOption(
        "--colour_fire-n-ice",
        "Show images in blue-grey-red colour scheme (highlight high & low signal)?"
    ) ? blue_grey_red :
    checkOption(
        "--colour_rainbow",
        "Show images in cyan-blue-black-red-yellow colour scheme?"
    ) ? rainbow :
    checkOption(
        "--colour_difference",
        "Show images in cyan-blue-black-red-yellow colour scheme (for difference images)?"
    ) ? cyan_black_yellow :
        greyscale;
}

inline void consume(std::string &src, std::vector<std::string> &dest, const std::string &delimiter) {
    size_t pos;
    while ((pos = src.find(delimiter)) != std::string::npos) {
        dest.push_back(src.substr(0, pos));
        src.erase(0, pos + delimiter.length());
    }
    dest.push_back(src);
}

std::vector<std::vector<std::string>> untangleDeviceIDs(std::string &tangled) {

    // Handle GPU (device) assignments for each rank, if specified
    std::vector<std::string> rankIDs;
    consume(tangled, rankIDs, ":");

    std::vector<std::vector<std::string>> untangled (rankIDs.size());
    // Now handle the thread assignments in each rank
    for (int i = 0; i < rankIDs.size(); i++) {
        consume(rankIDs[i], untangled[i], ",");
    }

    #ifdef DEBUG
    std::cout << "untangled.size() == " << untangled.size() << std::endl;
    for (int i = 0; i < untangled.size(); i++) {
        std::cout << "untangled[" << i << "]: ";
        for (const auto &thread : untangled[i])
            std::cout << thread << " ";
        std::cout << std::endl;
    }
    #endif

    return untangled;
}
