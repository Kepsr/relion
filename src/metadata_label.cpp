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
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "src/metadata_label.h"

// This is needed for static memory allocation
std::map<EMDL::EMDLabel, const EMDL::LabelData> EMDL::data;
std::map<std::string, EMDL::EMDLabel> EMDL::labels;
std::map<std::string, std::string> EMDL::definitions;
StaticInitialization EMDL::initialization;  // Just for initialization

struct EMDL::LabelData {

    const std::string name;
    const LabelType type;

    LabelData(): name{}, type{} {}  // Needed for LabelData to be put in a map

    LabelData(const std::string &s, EMDL::LabelType t): name{s}, type{t} {}

};

template <>
EMDL::LabelType EMDL::type2enum<int>() { return EMDL::INT; }

template <>
EMDL::LabelType EMDL::type2enum<bool>() { return EMDL::BOOL; }

template <>
EMDL::LabelType EMDL::type2enum<double>() { return EMDL::DOUBLE; }

template <>
EMDL::LabelType EMDL::type2enum<std::string>() { return EMDL::STRING; }

template <>
EMDL::LabelType EMDL::type2enum<std::vector<double> >() { return EMDL::DOUBLE_VECTOR; }

template <typename T>
void EMDL::addLabel(EMDLabel label, const std::string &name, const std::string &definition) {
    data.insert({label, {name, EMDL::type2enum<T>()}});
    labels[name] = label;
    definitions[name] = definition;
}

void EMDL::addAltLabel(EMDLabel label, std::string name) {
    labels[name] = label;
}

void EMDL::printDefinitions(std::ostream& out) {
    out << "+++ RELION MetaDataLabel (EMDL) definitions: +++" << std::endl;
    std::map<std::string, std::string>::const_iterator strIt;

    for (strIt = definitions.begin(); strIt != definitions.end(); strIt++) {
        out << std::setw(30) << strIt->first;

        const EMDL::EMDLabel label = labels[strIt->first];
        if (EMDL::is<int>(label)) {
            out << " (int)    ";
        } else if (EMDL::is<bool>(label)) {
            out << " (bool)   ";
        } else if (EMDL::is<double>(label)) {
            out << " (double) ";
        } else if (EMDL::is<std::string>(label)) {
            out << " (string) ";
        } else if (EMDL::is<std::vector<double>>(label)) {
            out << " (vector<double>) ";
        } else if (EMDL::is<void>(label)) {
            out << " (string) ";
        } else {
            REPORT_ERROR("EMDL::printDefinitions: unrecognised type");
        }

        out << ": " << strIt->second <<std::endl;
    }
}

EMDL::EMDLabel EMDL::str2Label(const std::string &labelName) {
    if (labels.find(labelName) == labels.end()) return EMDL::UNDEFINED;
    return labels[labelName];
}

std::string EMDL::label2Str(const EMDLabel &label) {
    if (data.find(label) == data.end()) return "";
    return data[label].name;
}

template <typename T>
bool EMDL::is(const EMDL::EMDLabel &label) {
    return data[label].type == EMDL::type2enum<T>();
}

bool EMDL::isValidLabel(const EMDLabel &label) {
    return label > EMDL::UNDEFINED && label < EMDL::LAST_LABEL;
}

bool EMDL::isValidLabel(const std::string &labelName) {
    EMDL::EMDLabel label = EMDL::str2Label(labelName);
    return EMDL::isValidLabel(label);
}
