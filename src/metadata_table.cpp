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
 * Authors:         J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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
 *    All comments concerning this program package may be sent to the
 *    e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/metadata_table.h"
#include "src/metadata_label.h"


inline std::string prependERROR(const std::string &s) {
    return "ERROR: " + s;
}

MetaDataTable::MetaDataTable():
    objects(0),
    label2offset(EMDL::LAST_LABEL, -1),
    current_objectID(0),
    doubleLabels(0),
    intLabels(0),
    boolLabels(0),
    stringLabels(0),
    doubleVectorLabels(0),
    unknownLabels(0),
    isList(false),
    name(""),
    comment(""),
    version(CURRENT_MDT_VERSION),
    activeLabels(0)
{}

MetaDataTable::MetaDataTable(const MetaDataTable &MD):
    objects(MD.objects.size()),
    label2offset(MD.label2offset),
    unknownLabelPosition2Offset(MD.unknownLabelPosition2Offset),
    unknownLabelNames(MD.unknownLabelNames),
    current_objectID(0),
    doubleLabels(MD.doubleLabels),
    intLabels(MD.intLabels),
    boolLabels(MD.boolLabels),
    stringLabels(MD.stringLabels),
    doubleVectorLabels(MD.doubleVectorLabels),
    unknownLabels(MD.unknownLabels),
    isList(MD.isList),
    name(MD.name),
    comment(MD.comment),
    version(MD.version),
    activeLabels(MD.activeLabels)
{
    for (size_t idx = 0; idx < MD.objects.size(); idx++) {
        objects[idx] = new MetaDataContainer(*(MD.objects[idx]));
        objects[idx]->table = this;
    }
}

MetaDataTable& MetaDataTable::operator = (const MetaDataTable &MD) {
    if (this != &MD) {
        clear();

        objects.resize(MD.objects.size());
        label2offset = MD.label2offset;
        unknownLabelPosition2Offset = MD.unknownLabelPosition2Offset;
        unknownLabelNames = MD.unknownLabelNames;
        current_objectID = 0;
        doubleLabels = MD.doubleLabels;
        intLabels = MD.intLabels;
        boolLabels = MD.boolLabels;
        stringLabels = MD.stringLabels;
        doubleVectorLabels = MD.doubleVectorLabels;
        unknownLabels = MD.unknownLabels;

        isList = MD.isList;
        name = MD.name;
        comment = MD.comment;
        version = MD.version;

        activeLabels = MD.activeLabels;

        for (long int idx = 0; idx < MD.objects.size(); idx++) {
            objects[idx] = new MetaDataContainer(this, MD.objects[idx]);
            objects[idx]->table = this;
        }
    }

    return *this;
}

MetaDataTable::~MetaDataTable() {
    for (long i = 0; i < objects.size(); i++) {
        delete objects[i];
    }
}

bool MetaDataTable::isEmpty() const {
    return objects.size() == 0;
}

size_t MetaDataTable::numberOfObjects() const {
    return objects.size();
}

void MetaDataTable::clear() {
    for (long i = 0; i < objects.size(); i++) {
        delete objects[i];
    }
    objects.clear();

    label2offset = std::vector<long>(EMDL::LAST_LABEL, -1);
    current_objectID = 0;
    unknownLabelPosition2Offset.clear();
    unknownLabelNames.clear();

    doubleLabels = 0;
    intLabels = 0;
    boolLabels = 0;
    stringLabels = 0;
    unknownLabels = 0;

    isList = false;
    name = "";
    comment = "";
    version = CURRENT_MDT_VERSION;

    activeLabels.clear();
}

void MetaDataTable::setComment(const std::string newComment) {
    comment = newComment;
}

std::string MetaDataTable::getComment() const {
    return comment;
}

bool MetaDataTable::containsComment() const {
    return comment != std::string("");
}

void MetaDataTable::setName(const std::string newName) {
    name = newName;
}

std::string MetaDataTable::getName() const {
    return name;
}

void MetaDataTable::setVersion(int v) {
    version = v;
}

int MetaDataTable::getVersion() const {
    return version;
}

int MetaDataTable::getCurrentVersion() {
    return CURRENT_MDT_VERSION;
}

std::string MetaDataTable::getUnknownLabelNameAt(int i) const {
    if (activeLabels[i] != EMDL::UNKNOWN_LABEL)
        REPORT_ERROR("MetaDataTable::getUnknownLabelNameAt(): the requested column is not an unknown label.");

    return unknownLabelNames[unknownLabelPosition2Offset[i]];
}

std::string MetaDataTable::getValueToString(EMDL::EMDLabel label, long objectID, bool escape) const {
    // SHWS 18 Jul 2018: this function previously had a stringstream, but it greatly slowed down
    // writing of large STAR files in some strange circumstances
    // (with large data.star and model.star files in refinement)
    // Therefore replaced the strstream with faster snprintf
    //
    // JZ 9 Aug 2018: still using a stringstream for vector<double> fields
    // => Avoid vector-valued columns in particle star-files.
    char buffer[14];

    std::string value;

    if (EMDL::isString(label)) {
        value = getValue<std::string>(label, objectID);
        if (escape) escapeStringForSTAR(value);
        return value;
    } else {
        if (EMDL::isDouble(label)) {
            double v = getValue<double>(label, objectID);

            if (abs(v) > 0.0 && abs(v) < 0.001 || abs(v) > 100000.0) {
                // If the magnitude of v is very small or very large, 
                // use floating-point form (6.02e-23).
                snprintf(buffer, 13, v < 0.0 ? "%12.5e" : "%12.6e", v);
            } else {
                // Otherwise, use fixed-point form (33.3333).
                snprintf(buffer, 13, v < 0.0 ? "%12.5f" : "%12.6f", v);
            }
        } else if (EMDL::isInt(label)) {
            long v = getValue<long>(label, objectID);
            snprintf(buffer, 13, "%12ld", v);
        } else if (EMDL::isBool(label)) {
            bool v = getValue<bool>(label, objectID);
            snprintf(buffer, 13, "%12d", (int) v);
        } else if (EMDL::isDoubleVector(label)) {
            std::vector<double> v = getValue<std::vector<double>>(label, objectID);

            if (v.size() == 0) {
                value = "[]";
            } else {
                std::stringstream sts;

                sts << std::setprecision(12);
                sts << '[';

                for (int i = 0; i < v.size() - 1; i++) {
                    sts << v[i] << ',';
                }
                sts << v[v.size() - 1];

                sts << ']';

                value = sts.str();
            }
            return value;
        }
        return std::string(buffer);
    }
}

bool MetaDataTable::setUnknownValue(int labelPosition, const std::string &value) {
    long offset = unknownLabelPosition2Offset[labelPosition];
    if (offset < 0) REPORT_ERROR("MetaDataTable::setValueFromString BUG: offset should not be negative here....");

    if (offset > -1) {
        objects[current_objectID]->unknowns[offset] = value;
        return true;
    }
    return false;
}

bool MetaDataTable::setValueFromString(
    EMDL::EMDLabel label, const std::string &value, long int objectID
) {
    if (EMDL::isString(label)) {
        return setValue(label, value, objectID);
    } else {
        std::istringstream i(value);

        if (EMDL::isDouble(label)) {
            double v;
            i >> v;
            return setValue(label, v, objectID);
        } else if (EMDL::isInt(label)) {
            long v;
            i >> v;
            return setValue(label, v, objectID);
        } else if (EMDL::isBool(label)) {
            bool v;
            i >> v;
            return setValue(label, v, objectID);
        } else if (EMDL::isDoubleVector(label)) {
            std::vector<double> v;
            v.reserve(32);

            char* temp = new char[value.size()+1];
            strcpy(temp, value.c_str());

            char* token;
            char* rest = temp;

            while ((token = strtok_r(rest, "[,]", &rest)) != 0) {
                double d;
                std::stringstream sts(token);
                sts >> d;

                v.push_back(d);
            }

            delete[] temp;

            return setValue(label, v, objectID);
        }
    }

    REPORT_ERROR("Logic error: should not happen");
    return false;
}

// Comparators used for sorting
namespace MD {

struct DoubleComparator {
    DoubleComparator(long index) : index(index) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->doubles[index] < rh->doubles[index];
    }

    long index;
};

struct IntComparator {

    IntComparator(long index): index(index) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->ints[index] < rh->ints[index];
    }

    long index;
};

struct StringComparator {

    long index;

    StringComparator(long index): index(index) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->strings[index] < rh->strings[index];
    }

};

// Can we get this to inherit from StringComperator?
struct StringAfterAtComparator {

    long index;

    StringAfterAtComparator(long index): index(index) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        std::string slh = lh->strings[index];
        std::string srh = rh->strings[index];
        slh = slh.substr(slh.find("@") + 1);
        srh = srh.substr(srh.find("@") + 1);
        return slh < srh;
    }

};

// Can we get this to inherit from StringComperator?
struct StringBeforeAtComparator {

    long index;

    StringBeforeAtComparator(long index): index(index) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        std::string slh = lh->strings[index];
        std::string srh = rh->strings[index];
        slh = slh.substr(0, slh.find("@"));
        srh = srh.substr(0, srh.find("@"));
        std::stringstream stslh, stsrh;
        stslh << slh;
        stsrh << srh;
        long ilh, irh;
        stslh >> ilh;
        stsrh >> irh;

        return ilh < irh;
    }

};

};

void MetaDataTable::sort(
    EMDL::EMDLabel name, bool do_reverse, bool only_set_index, bool do_random
) {
    if (do_random) {
        srand(time(NULL)); // initialise random seed
    } else if (!EMDL::isNumber(name)) {
        REPORT_ERROR("MetadataTable::sort%% " + prependERROR("can only sorted numbers"));
    }

    std::vector<std::pair<double, long int>> vp;
    vp.reserve(objects.size());
    long int i = 0;

    FOR_ALL_OBJECTS_IN_METADATA_TABLE(*this) {
        double dval;
        if (do_random) {
            dval = (double) rand();
        } else if (EMDL::isInt(name)) {
            dval = (double) getValue<long>(name);
        } else {
            // EMDL::isDouble(name)
            dval = getValue<double>(name);
        }

        vp.push_back(std::make_pair(dval, i));
        i++;
    }

    std::sort(vp.begin(), vp.end());
    if (do_reverse && !do_random)
        std::reverse(vp.begin(), vp.end());

    if (only_set_index) {
        // Add an extra column with the sorted position of each entry
        for (long j = 0; j < vp.size(); j++) {
            (*this).setValue(EMDL::SORTED_IDX, j, vp[j].second);
        }
    } else {
        // Change the actual order in the MetaDataTable
        std::vector<MetaDataContainer*> objs(objects.size());

        for (long j = 0; j < vp.size(); j++) {
            objs[j] = objects[vp[j].second];
        }

        objects = objs;
    }
    // reset pointer to the beginning of the table
    firstObject();
}

void MetaDataTable::newSort(const EMDL::EMDLabel label, bool do_reverse, bool do_sort_after_at, bool do_sort_before_at) {

    if (EMDL::isString(label)) {
        if (do_sort_after_at) {
            std::stable_sort(
                objects.begin(), objects.end(),
                MD::StringAfterAtComparator(label2offset[label])
            );
        } else if (do_sort_before_at) {
            std::stable_sort(
                objects.begin(), objects.end(),
                MD::StringBeforeAtComparator(label2offset[label])
            );
        } else {
            std::stable_sort(
                objects.begin(), objects.end(),
                MD::StringComparator(label2offset[label])
            );
        }
    } else if (EMDL::isDouble(label)) {
        std::stable_sort(
            objects.begin(), objects.end(),
            MD::DoubleComparator(label2offset[label])
        );
    } else if (EMDL::isInt(label)) {
        std::stable_sort(
            objects.begin(), objects.end(),
            MD::IntComparator(label2offset[label])
        );
    } else {
        REPORT_ERROR("Cannot sort this label: " + EMDL::label2Str(label));
    }

    if (do_reverse) {
        std::reverse(objects.begin(), objects.end());
    }
}

// Will be removed in 3.2
bool MetaDataTable::labelExists(EMDL::EMDLabel name) const {
    return containsLabel(name);
}

bool MetaDataTable::containsLabel(const EMDL::EMDLabel label, std::string unknownLabel) const {
    for (int i = 0; i < activeLabels.size(); i++) {
        if (
            activeLabels[i] == label &&
            (label != EMDL::UNKNOWN_LABEL || getUnknownLabelNameAt(i) == unknownLabel)
        ) return true;
    }
    return false;
}

std::vector<EMDL::EMDLabel> MetaDataTable::getActiveLabels() const {
    return activeLabels;
}

void MetaDataTable::deactivateLabel(EMDL::EMDLabel label, std::string unknownLabel) {
    for (int i = 0; i < activeLabels.size(); i++) {
        if (
            activeLabels[i] == label &&
            (label != EMDL::UNKNOWN_LABEL || unknownLabelNames[unknownLabelPosition2Offset[i]] == unknownLabel)
        ) {
            activeLabels.erase(activeLabels.begin() + i);
            unknownLabelPosition2Offset.erase(unknownLabelPosition2Offset.begin() + i);

            if (label != EMDL::UNKNOWN_LABEL)
                label2offset[label] = -1;
        }
    }
}

void MetaDataTable::addLabel(EMDL::EMDLabel label, std::string unknownLabel) {
    if (label >= EMDL::LAST_LABEL)
        REPORT_ERROR(std::string(
            "MetaDataTable::addLabel: unrecognised label: "
        ) + EMDL::label2Str(label));
    if (label == EMDL::UNKNOWN_LABEL && unknownLabel.empty())
        REPORT_ERROR("MetaDataTable::addLabel: unknownLabel is empty");

    if (label2offset[label] < 0 || label == EMDL::UNKNOWN_LABEL) {
        // keep pushing the same unknown label...
        long id;

        if (EMDL::isDouble(label)) {
            id = doubleLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->doubles.push_back(0);
            }

            doubleLabels++;
        } else if (EMDL::isInt(label)) {
            id = intLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->ints.push_back(0);
            }

            intLabels++;
        } else if (EMDL::isBool(label)) {
            id = boolLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->bools.push_back(false);
            }

            boolLabels++;
        } else if (EMDL::isString(label)) {
            id = stringLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->strings.push_back("empty");
            }

            stringLabels++;
        } else if (EMDL::isDoubleVector(label)) {
            id = doubleVectorLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->doubleVectors.push_back(std::vector<double>());
            }

            doubleVectorLabels++;
        } else if (EMDL::isUnknown(label)) {
            id = unknownLabels;

            for (long i = 0; i < objects.size(); i++) {
                objects[i]->unknowns.push_back("empty");
            }

            unknownLabelNames.push_back(unknownLabel);
            unknownLabels++;
        }

        activeLabels.push_back(label);
        unknownLabelPosition2Offset.push_back(EMDL::isUnknown(label) ? id : -1);

        label2offset[label] = id;
    }
}

void MetaDataTable::addMissingLabels(const MetaDataTable* mdt) {
    for (long i = 0; i < mdt->activeLabels.size(); i++) {
        EMDL::EMDLabel l = mdt->activeLabels[i];

        if (l == EMDL::UNKNOWN_LABEL) {
            std::string unknownLabel = mdt->getUnknownLabelNameAt(i);
            if (!containsLabel(l, unknownLabel)) {
                addLabel(l, unknownLabel);
            }
        } else if (label2offset[l] < 0) {
            addLabel(l);
        }
    }
}

void MetaDataTable::append(const MetaDataTable& mdt) {
    if (activeLabels.size() == 0) {
        // If the current one is empty, add missing labels and append the new one:
        addMissingLabels(&mdt);
    } else {
        // If the current one is not-empty, check all labels are the same before appending. Otherwise, raise error
        if (!compareLabels(*this, mdt))
            REPORT_ERROR("ERROR in appending metadata tables with not the same columns!");
    }

    // Now append
    objects.reserve(objects.size() + mdt.numberOfObjects());
    for (long i = 0; i < mdt.objects.size(); i++) {
        objects.push_back(new MetaDataContainer(
            this,
            doubleLabels, intLabels, boolLabels, stringLabels,
            doubleVectorLabels, unknownLabels
        ));

        setObjectUnsafe(mdt.getObject(i), objects.size() - 1);
    }

    // reset pointer to the beginning of the table
    firstObject();
}


MetaDataContainer* MetaDataTable::getObject(long objectID) const {
    if (objectID < 0) {
        objectID = current_objectID;
    }

    checkObjectID(objectID, "MetaDataTable::getObject");

    return objects[objectID];
}

void MetaDataTable::setObject(MetaDataContainer* data, long objectID) {
    if (objectID < 0) {
        objectID = current_objectID;
    }

    checkObjectID(objectID, "MetaDataTable::setObject");
    addMissingLabels(data->table);

    setObjectUnsafe(data, objectID);
}

void MetaDataTable::setValuesOfDefinedLabels(MetaDataContainer* data, long objectID) {
    if (objectID < 0) {
        objectID = current_objectID;
    }

    checkObjectID(objectID, "MetaDataTable::setValuesOfDefinedLabels");

    setObjectUnsafe(data, objectID);
}

void MetaDataTable::reserve(size_t capacity) {
    objects.reserve(capacity);
}

void MetaDataTable::setObjectUnsafe(MetaDataContainer* data, long objectID) {
    MetaDataContainer* obj = objects[objectID];

    for (long i = 0; i < data->table->activeLabels.size(); i++) {
        EMDL::EMDLabel label = data->table->activeLabels[i];

        if (label != EMDL::UNKNOWN_LABEL) {
            long myOff = label2offset[label];
            long srcOff = data->table->label2offset[label];

            if (myOff < 0)
                continue;

            if (EMDL::isDouble(label)) {
                obj->doubles[myOff] = data->doubles[srcOff];
            } else if (EMDL::isInt(label)) {
                obj->ints[myOff] = data->ints[srcOff];
            } else if (EMDL::isBool(label)) {
                obj->bools[myOff] = data->bools[srcOff];
            } else if (EMDL::isString(label)) {
                obj->strings[myOff] = data->strings[srcOff];
            } else if (EMDL::isDoubleVector(label)) {
                obj->doubleVectors[myOff] = data->doubleVectors[srcOff];
            }
        } else {
            std::string unknownLabel = data->table->getUnknownLabelNameAt(i);
            long srcOff = data->table->unknownLabelPosition2Offset[i];
            long myOff = -1;

            for (int j = 0; j < unknownLabelNames.size(); j++) {
                if (unknownLabelNames[j] == unknownLabel) {
                    myOff = j;
                    break;
                }
            }

            if (myOff < 0)
                REPORT_ERROR("MetaDataTable::setObjectUnsafe: logic error. cannot find srcOff.");

            obj->unknowns[myOff] = data->unknowns[srcOff];
        }
    }
}

void MetaDataTable::addObject() {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels, stringLabels,
        doubleVectorLabels, unknownLabels
    ));

    current_objectID = objects.size()-1;
}

void MetaDataTable::addObject(MetaDataContainer* data) {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels, stringLabels,
        doubleVectorLabels, unknownLabels
    ));

    setObject(data, objects.size() - 1);
    current_objectID = objects.size() - 1;
}

void MetaDataTable::addValuesOfDefinedLabels(MetaDataContainer* data) {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels,
        stringLabels, doubleVectorLabels, unknownLabels
    ));

    setValuesOfDefinedLabels(data, objects.size()-1);
    current_objectID = objects.size()-1;
}

void MetaDataTable::removeObject(long objectID) {
    long i = objectID < 0 ? current_objectID : objectID;

    checkObjectID(i, "MetaDataTable::removeObject");

    delete objects[i];
    objects.erase(objects.begin() + i);

    current_objectID = objects.size() - 1;
}

long int MetaDataTable::firstObject() {
    current_objectID = 0;
    return 0;
}

long int MetaDataTable::nextObject() {
    current_objectID++;

    if (current_objectID >= objects.size()) {
        return NO_MORE_OBJECTS;
    } else {
        return current_objectID;
    }
}

long int MetaDataTable::goToObject(long int objectID) {
    checkObjectID(objectID, "MetaDataTable::goToObject");

    current_objectID = objectID;
    return current_objectID;
}

long int MetaDataTable::readStarLoop(std::ifstream &in, bool do_only_count) {
    isList = false;

    //Read column labels
    int labelPosition = 0;
    std::string line, token;

    // First read all the column labels
    while (getline(in, line, '\n')) {
        line = simplify(line);
        // TODO: handle comments...
        if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
            continue;

        if (line[0] == '_') {
            // label definition line
            // Take string from "_" until "#"
            size_t start = line.find("_");
            size_t end = line.find("#");

            token = line.substr(start + 1, end - start - 2);

            EMDL::EMDLabel label = EMDL::str2Label(token);

            if (label == EMDL::UNDEFINED) {
                std::cerr << " + WARNING: will ignore (but maintain) values for the unknown label: " << token << std::endl;
                label = EMDL::UNKNOWN_LABEL;
            }

            addLabel(label, token);

            labelPosition++;
        } else {
            // found first data line
            break;
        }
    }

    // Then fill the table (dont read another line until the one from above has been handled)
    bool is_first = true;
    long int nr_objects = 0;
    const int num_labels = activeLabels.size();

    while (is_first || getline(in, line, '\n')) {
        is_first = false;

        line = simplify(line);
        // Stop at empty line
        if (line[0] == '\0')
            break;

        nr_objects++;
        if (!do_only_count) {
            // Add a new line to the table
            addObject();

            // Parse data values
            int pos = 0;
            std::string value;
            labelPosition = 0;
            while (nextTokenInSTAR(line, pos, value)) {
                if (labelPosition >= num_labels) {
                    std::cerr << "Error in line: " << line << std::endl;
                    REPORT_ERROR("A line in the STAR file contains more columns than the number of labels.");
                }
                // Check whether this is an unknown label
                if (activeLabels[labelPosition] == EMDL::UNKNOWN_LABEL) {
                    setUnknownValue(labelPosition, value);
                } else {
                    setValueFromString(activeLabels[labelPosition], value);
                }
                labelPosition++;
            }
            if (labelPosition < num_labels && num_labels > 2) {
                // For backward-compatibility for cases like "fn_mtf <empty>", don't die if num_labels == 2.
                std::cerr << "Error in line: " << line << std::endl;
                REPORT_ERROR("A line in the STAR file contains fewer columns than the number of labels. Expected = " + integerToString(num_labels) + " Found = " +  integerToString(labelPosition));
            }
        }
    }

    return nr_objects;
}

bool MetaDataTable::readStarList(std::ifstream &in) {
    isList = true;
    addObject();
    long int objectID = objects.size() - 1;

    std::string line, firstword, value;

    bool also_has_loop = false;

    // Read data and fill structures accordingly
    int labelPosition = 0;
    while (getline(in, line, '\n')) {
        int pos = 0;
        // Ignore empty lines
        if (!nextTokenInSTAR(line, pos, firstword))
            continue;

        // Get label-value pairs
        if (firstword[0] == '_') {
            std::string token = firstword.substr(1); // get rid of leading underscore
            EMDL::EMDLabel label = EMDL::str2Label(token);
            if (!nextTokenInSTAR(line, pos, value))
                REPORT_ERROR("MetaDataTable::readStarList: did not encounter a single word after "+firstword);

            if (label == EMDL::UNDEFINED) {
                label = EMDL::UNKNOWN_LABEL;
                addLabel(label, token);
                setUnknownValue(labelPosition, value);
                std::cerr << " + WARNING: will ignore (but maintain) values for the unknown label: " << token << std::endl;
            } else {
                addLabel(label);
                setValueFromString(label, value, objectID);
            }
            labelPosition++;
        } else if (firstword[0] == '#' || firstword[0] == ';') {
            // Check whether there is a comment or an empty line
            // TODO: handle comments?
            continue;
        } else if (firstword.find("loop_") == 0) {
            // Check whether a loop structure comes after this list
            also_has_loop = true;
            return also_has_loop;
        } else if (firstword.find("data_") == 0) {
            // Check whether this data blocks ends (because a next one is there)
            // Should I reverse the pointer one line?
            return also_has_loop;
        }
    }
    // Reached the end of the file
    return also_has_loop;
}

long int MetaDataTable::readStar(
    std::ifstream &in, const std::string &name, bool do_only_count
) {
    std::string line, token, value;
    clear();
    bool also_has_loop;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Set the version to 30000 by default, in case there is no version tag
    // (version tags were introduced in version 31000)
    version = 30000;

    // Proceed until the next data_ or _loop statement
    // The loop statement may be necessary for data blocks that have a list AND a table inside them
    while (getline(in, line, '\n')) {
        trim(line);
        if (line.find("# version ") != std::string::npos) {
            token = line.substr(line.find("# version ") + std::string("# version ").length());

            std::istringstream sts(token);
            sts >> version;
        }

        // Find data_ lines
        if (line.find("data_") != std::string::npos) {
            token = line.substr(line.find("data_") + 5);
            // If a name has been given, only read data_thatname
            // Otherwise, just read the first data_ block
            if (name.empty() || name == token) {
                setName(token);
                // Get the next item that starts with "_somelabel" or with "loop_"
                int current_pos = in.tellg();
                while (getline(in, line, '\n')) {
                    if (line.find("loop_") != std::string::npos) {
                        return readStarLoop(in, do_only_count);
                    } else if (line[0] == '_') {
                        // go back one line in the ifstream
                        in.seekg(current_pos);
                        also_has_loop = readStarList(in);
                        return !also_has_loop;
                    }
                }
            }
        }
    }

    // Clear the eofbit so we can perform more actions on the stream.
    in.clear();

    return 0;
}

long int MetaDataTable::read(
    const FileName &filename, const std::string &name, bool do_only_count
) {

    // Clear current table
    clear();

    // Check for an :star extension
    FileName fn_read = filename.removeFileFormat();

    std::ifstream in(fn_read.data(), std::ios_base::in);

    if (in.fail()) {
        REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_read + " does not exist" );
    }

    return readStar(in, name, do_only_count);

    in.close();

    // Go to the first object
    firstObject();
}

void MetaDataTable::write(std::ostream& out) {
    // Only write tables that have something in them
    if (isEmpty()) {
        return;
    }

    if (version >= 30000) {
        out << "\n";
        out << "# version " << getCurrentVersion() <<"\n";
    }

    out << "\n";
    out << "data_" << getName() <<"\n";

    if (containsComment()) {
        out << "# "<< comment << "\n";
    }

    out << "\n";

    if (!isList) {
        // Write loop header structure
        out << "loop_ \n";

        for (long i = 0, n_printed = 1; i < activeLabels.size(); i++) {
            EMDL::EMDLabel l = activeLabels[i];
            if (l == EMDL::UNKNOWN_LABEL) {
                const long offset = unknownLabelPosition2Offset[i];
                out << "_" << unknownLabelNames[offset]<< " #" << (n_printed++) << " \n";
            } else if (l != EMDL::COMMENT && l != EMDL::SORTED_IDX) {
                // EMDL::SORTED_IDX is only for internal use, never write it out!
                out << "_" << EMDL::label2Str(l) << " #" << (n_printed++) << " \n";
            }
        }

        // Write actual data block
        for (long int idx = 0; idx < objects.size(); idx++) {
            std::string entryComment = "";

            for (long i = 0; i < activeLabels.size(); i++) {
                EMDL::EMDLabel l = activeLabels[i];

                if (l == EMDL::UNKNOWN_LABEL) {
                    out.width(10);
                    std::string token, val;
                    long offset = unknownLabelPosition2Offset[i];
                    val = objects[idx]->unknowns[offset];
                    escapeStringForSTAR(val);
                    out << val << " ";
                } else if (l != EMDL::COMMENT && l != EMDL::SORTED_IDX) {
                    out.width(10);
                    std::string val = getValueToString(l, idx, true); // escape=true
                    out << val << " ";
                }
                if (l == EMDL::COMMENT) {
                    entryComment = getValue<std::string>(EMDL::COMMENT, idx);
                }
            }
            if (entryComment != std::string("")) {
                out << "# " << entryComment;
            }
            out << "\n";
        }
        // Finish table with a white-line
        out << " \n";

    } else {
        // isList
        // Get first object. In this case (row format) there is a single object
        std::string entryComment = "";
        int maxWidth = 10;

        for (long i = 0; i < activeLabels.size(); i++) {
            EMDL::EMDLabel l = activeLabels[i];

            if (l != EMDL::COMMENT && l != EMDL::UNKNOWN_LABEL) {
                int w = EMDL::label2Str(l).length();
                if (w > maxWidth) maxWidth = w;
            } else if (l == EMDL::UNKNOWN_LABEL) {
                long offset = unknownLabelPosition2Offset[i];
                int w = unknownLabelNames[offset].length();
                if (w > maxWidth) maxWidth = w;
            } else {
                entryComment = getValue<std::string>(EMDL::COMMENT, 0);
            }
        }

        for (long i = 0; i < activeLabels.size(); i++) {
            EMDL::EMDLabel l = activeLabels[i];

            if (l == EMDL::UNKNOWN_LABEL) {
                long offset = unknownLabelPosition2Offset[i];
                int w = unknownLabelNames[offset].length();
                out << "_" << unknownLabelNames[offset] << std::setw(12 + maxWidth - w) << " " << objects[0]->unknowns[offset] << "\n";
            } else if (l != EMDL::COMMENT) {
                int w = EMDL::label2Str(l).length();
                out << "_" << EMDL::label2Str(l) << std::setw(12 + maxWidth - w) << " ";

                std::string val = getValueToString(l, 0, true); // escape=true
                out << val << "\n";
            }
        }
        if (entryComment != std::string("")) {
            out << "# " << entryComment << "\n";
        }

        // End a data block with a white line
        out << " \n";
    }
}

void MetaDataTable::write(const FileName &fn_out) {
    std::ofstream fh;
    FileName fn_tmp = fn_out + ".tmp";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MetaDataTable::write: cannot write to file: " + fn_out);
    //    fh << "# RELION; version " << g_RELION_VERSION << std::endl;
    write(fh);
    fh.close();
    // Rename to prevent errors with programs in pipeliner reading in incomplete STAR files
    std::rename(fn_tmp.c_str(), fn_out.c_str());

}

void MetaDataTable::columnHistogram(
    EMDL::EMDLabel label, std::vector<RFLOAT> &histX, std::vector<RFLOAT> &histY,
    int verb, CPlot2D *plot2D,
    long int nr_bin, RFLOAT hist_min, RFLOAT hist_max,
    bool do_fractional_instead, bool do_cumulative_instead
) {
    if (!containsLabel(label))
        REPORT_ERROR(prependERROR("The column specified is not present in the MetaDataTable."));

    std::vector<RFLOAT> values;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(*this) {
        RFLOAT val;
        if (EMDL::isDouble(label)) {
            val = getValue<double>(label);
        } else if (EMDL::isInt(label)) {
            val = getValue<long>(label);
        } else if (EMDL::isBool(label)) {
            val = getValue<bool>(label);
        } else {
            REPORT_ERROR("Cannot use --stat_column for this type of column");
        }
        values.push_back(val);
    }

    std::string title = EMDL::label2Str(label);
    histogram(values, histX, histY, verb, title, plot2D, nr_bin, hist_min, hist_max, do_fractional_instead, do_cumulative_instead);
}

void MetaDataTable::histogram(
    std::vector<RFLOAT> &values, std::vector<RFLOAT> &histX, std::vector<RFLOAT> &histY,
    int verb, std::string title, CPlot2D *plot2D,
    long int nr_bin, RFLOAT hist_min, RFLOAT hist_max,
    bool do_fractional_instead, bool do_cumulative_instead
) {
    double sum = 0, sumsq = 0;
    for (size_t i = 0, ilim = values.size(); i < ilim; i++) {
        RFLOAT value = values[i];
        sum += value;
        sumsq += value * value;
    }

    long long n_row = values.size();
    std::sort(values.begin(), values.end());
    sum /= n_row; sumsq /= n_row;

    if (verb > 0) {
        std::cout << "Number of items: " << n_row << std::endl;
        std::cout << "Min: " << values[0] << " Q1: " << values[n_row / 4];
        std::cout << " Median: " << values[n_row / 2] << " Q3: " << values[n_row * 3 / 4] << " Max: " << values[n_row - 1] << std::endl;
        std::cout << "Mean: " << sum << " Std: " << std::sqrt(sumsq - sum * sum) << std::endl;
    }

    RFLOAT iqr = values[n_row * 3 / 4] - values[n_row / 2];
    RFLOAT bin_width = 1;
    unsigned int bin_size = 1;

    // change bin parameters only when there are many values
    if (iqr != 0) {
        if (nr_bin <= 0) {
            hist_min = values[0];
            hist_max = values[n_row - 1];
            bin_width = 2 * iqr / std::pow(n_row, 1.0 / 3); // Freedman-Diaconis rule
            bin_size = (unsigned int)(std::ceil((hist_max - hist_min) / bin_width));
            if (bin_size > 5000) bin_size = 5000; // FIXME: Ad hoc upper limit to avoid using too much memory
        } else {
            if (!std::isfinite(hist_min) || hist_min == -LARGE_NUMBER) { hist_min = values[0]; }
            if (!std::isfinite(hist_max) || hist_max == +LARGE_NUMBER) { hist_max = values[n_row - 1]; }
            bin_size = nr_bin;
        }
        bin_width = (hist_max - hist_min) / bin_size;
    } else {
        if (!std::isfinite(hist_min) || hist_min == -LARGE_NUMBER) { hist_min = values[0]; }
        if (!std::isfinite(hist_max) || hist_max == +LARGE_NUMBER) { hist_max = values[n_row - 1]; }
    }

    bin_size += 2; // for -inf and +inf
    if (verb > 0) std::cout << "Bin size: " << bin_size << " width: " << bin_width << std::endl;

    std::vector<long> hist(bin_size);
    histY.resize(4 * bin_size, 0.0);
    histX.resize(4 * bin_size, 0.0);
    for (int i = 0; i < n_row; i++) {
        int ibin = (values[i] - hist_min) / bin_width + 1;
        if (ibin < 0) { ibin = 0; }
        if (ibin >= bin_size) { ibin = bin_size - 1; }
        hist[ibin]++;
    }

    long cum = 0;
    for (int i = 0; i < bin_size; i++) {
        if (i == 0) {
            if (verb > 0) std::cout << "[-INF, " << hist_min << "): ";
            histX[4 * i]     = hist_min - bin_width;
            histX[4 * i + 1] = hist_min - bin_width;
            histX[4 * i + 2] = hist_min;
            histX[4 * i + 3] = hist_min;
        } else if (i == bin_size - 1) {
            if (verb > 0) std::cout << "[" << hist_max << ", +INF]: ";
            histX[4 * i]     = hist_max;
            histX[4 * i + 1] = hist_max;
            histX[4 * i + 2] = hist_max + bin_width;
            histX[4 * i + 3] = hist_max + bin_width;
        } else {
            if (verb > 0) std::cout << "[" << (hist_min + bin_width * (i - 1)) << ", " << (hist_min + bin_width * i) << "): ";
            histX[4 * i]     = hist_min + bin_width * (i - 1);
            histX[4 * i + 1] = hist_min + bin_width * (i - 1);
            histX[4 * i + 2] = hist_min + bin_width * i;
            histX[4 * i + 3] = hist_min + bin_width * i;
        }

        cum += hist[i];
        if (do_fractional_instead) {
            hist[i] = 100.0 * hist[i] / (float) n_row;
        } else if (do_cumulative_instead) {
            hist[i] = 100 * cum / (float) n_row;
        }

        if (verb > 0) { std::cout  << hist[i] << std::endl; }

        histY[4 * i + 1] = histY[4 * i + 2] = hist[i];
        histY[4 * i] = histY[4 * i + 3] = 0.0;

    }
    histX[histX.size() - 1] = histX[histX.size() - 2];

    if (plot2D != NULL) {
        plot2D->SetTitle(" Histogram of " + title);
        plot2D->SetDrawLegend(false);
        plot2D->AddDataSet(histX, histY);
        plot2D->SetXAxisTitle(title);
        plot2D->SetYAxisTitle("# entries");
    }
}

void MetaDataTable::addToCPlot2D(
    CPlot2D *plot2D, EMDL::EMDLabel xaxis, EMDL::EMDLabel yaxis,
    double red, double green, double blue, double linewidth, std::string marker
) {
    CDataSet dataSet;
    if (marker.empty()) {
        dataSet.SetDrawMarker(false);
    } else {
        dataSet.SetDrawMarker(true);
        dataSet.SetMarkerSymbol(marker);
    }
    dataSet.SetLineWidth(linewidth);
    dataSet.SetDatasetColor(red, green, blue);
    dataSet.SetDatasetTitle(EMDL::label2Str(yaxis));

    double mydbl;
    long int myint;
    double xval, yval;
    for (long int idx = 0; idx < objects.size(); idx++) {
        const long offx = label2offset[xaxis];
        if (offx < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find x-axis label");

        if (xaxis == EMDL::UNDEFINED) {
            xval = idx + 1;
        } else if (EMDL::isDouble(xaxis)) {
            objects[idx]->getValue(offx, mydbl);
            xval = mydbl;
        } else if (EMDL::isInt(xaxis)) {
            objects[idx]->getValue(offx, myint);
            xval = myint;
        } else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot x-axis double, int or long int");

        const long offy = label2offset[yaxis];
        if (offy < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find y-axis label");

        if (EMDL::isDouble(yaxis)) {
            objects[idx]->getValue(offy, mydbl);
            yval = mydbl;
        } else if (EMDL::isInt(yaxis)) {
            objects[idx]->getValue(offy, myint);
            yval = myint;
        } else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot y-axis double, int or long int");

        CDataPoint point(xval, yval);
        dataSet.AddDataPoint(point);

    }

    plot2D->AddDataSet(dataSet);

    if (xaxis != EMDL::UNDEFINED) {
        plot2D->SetXAxisTitle(EMDL::label2Str(xaxis));
    }
    plot2D->SetYAxisTitle(EMDL::label2Str(yaxis));

}

void MetaDataTable::printLabels(std::ostream &ost) {
    for (int i = 0; i < activeLabels.size(); i++) {
        ost << EMDL::label2Str(activeLabels[i]) << "\n";
    }
}

void MetaDataTable::randomiseOrder() {
    std::random_shuffle(objects.begin(), objects.end());
}

void MetaDataTable::checkObjectID(long id, std::string caller) const {
    if (id >= objects.size() || id < 0) {
        std::stringstream sts0, sts1;
        sts0 << id;
        sts1 << objects.size();
        REPORT_ERROR(
            caller + ": object " + sts0.str() + " out of bounds! (" + sts1.str() + " objects present)"
        );
    }
}

//FIXME: does not support unknownLabels but this function is only used by relion_star_handler
//       so I will leave this for future...
void compareMetaDataTable(MetaDataTable &MD1, MetaDataTable &MD2,
    MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
    EMDL::EMDLabel label1, double eps, EMDL::EMDLabel label2, EMDL::EMDLabel label3
) {
    if (!MD1.containsLabel(label1))
        REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label1.");
    if (!MD2.containsLabel(label1))
        REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label1.");

    if (label2 != EMDL::UNDEFINED) {
        if (!EMDL::isDouble(label1) || !EMDL::isDouble(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 2D or 3D distances are only allowed for doubles.");
        if (!MD1.containsLabel(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label2.");
        if (!MD2.containsLabel(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label2.");
    }

    if (label3 != EMDL::UNDEFINED) {
        if (!EMDL::isDouble(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 3D distances are only allowed for doubles.");
        if (!MD1.containsLabel(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label3.");
        if (!MD2.containsLabel(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label3.");
    }

    MDboth.clear();
    MDonly1.clear();
    MDonly2.clear();

    std::string mystr1, mystr2;
    long int myint1, myint2;
    double myd1, myd2, mydy1 = 0.0, mydy2 = 0.0, mydz1 = 0.0, mydz2 = 0.0;

    // loop over MD1
    std::vector<long int> to_remove_from_only2;
    for (
        long int current_object1 = MD1.firstObject();
        current_object1 != MetaDataTable::NO_MORE_OBJECTS && current_object1 != MetaDataTable::NO_OBJECTS_STORED;
        current_object1 = MD1.nextObject()
    ) {
        if (EMDL::isString(label1)) {
            mystr1 = MD1.getValue<std::string>(label1);
        } else if (EMDL::isInt(label1)) {
            myint1 = MD1.getValue<int>(label1);
        } else if (EMDL::isDouble(label1)) {
            myd1 = MD1.getValue<double>(label1);
            if (label2 != EMDL::UNDEFINED)
            mydy1 = MD1.getValue<double>(label2);
            if (label3 != EMDL::UNDEFINED)
            mydz1 = MD1.getValue<double>(label3);
        } else {
            REPORT_ERROR("compareMetaDataTableEqualLabel ERROR: only implemented for strings, integers or doubles");
        }
        // loop over MD2
        bool have_in_2 = false;
        for (
            long int current_object2 = MD2.firstObject();
            current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object2 != MetaDataTable::NO_OBJECTS_STORED;
            current_object2 = MD2.nextObject()
        ) {
            if (EMDL::isString(label1)) {
                mystr2 = MD2.getValue<std::string>(label1);
                if (strcmp(mystr1.c_str(), mystr2.c_str()) == 0) {
                    have_in_2 = true;
                    to_remove_from_only2.push_back(current_object2);
                    MDboth.addObject(MD1.getObject());
                    break;
                }
            } else if (EMDL::isInt(label1)) {
                myint2 = MD2.getValue<int>(label1);
                if (abs(myint2 - myint1) <= round(eps)) {
                    have_in_2 = true;
                    to_remove_from_only2.push_back(current_object2);
                    MDboth.addObject(MD1.getObject());
                    break;
                }
            } else if (EMDL::isDouble(label1)) {
                myd2 = MD2.getValue<double>(label1);
                if (label2 != EMDL::UNDEFINED)
                mydy2 = MD2.getValue<double>(label2);
                if (label3 != EMDL::UNDEFINED)
                mydz2 = MD2.getValue<double>(label3);

                double dist = sqrt(
                    (myd1  - myd2)  * (myd1  - myd2)  +
                    (mydy1 - mydy2) * (mydy1 - mydy2) +
                    (mydz1 - mydz2) * (mydz1 - mydz2)
                );
                if (abs(dist) <= eps) {
                    have_in_2 = true;
                    to_remove_from_only2.push_back(current_object2);
                    //std::cerr << " current_object1= " << current_object1 << std::endl;
                    //std::cerr << " myd1= " << myd1 << " myd2= " << myd2 << " mydy1= " << mydy1 << " mydy2= " << mydy2 << " dist= "<<dist<<std::endl;
                    //std::cerr << " to be removed current_object2= " << current_object2 << std::endl;
                    MDboth.addObject(MD1.getObject());
                    break;
                }
            }
        }

        if (!have_in_2) {
            MDonly1.addObject(MD1.getObject());
        }
    }

    for (
        long int current_object2 = MD2.firstObject();
        current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object2 != MetaDataTable::NO_OBJECTS_STORED;
        current_object2 = MD2.nextObject()
    ) {

        bool to_be_removed = false;
        for (long int i = 0; i < to_remove_from_only2.size(); i++) {
            if (to_remove_from_only2[i] == current_object2) {
                to_be_removed = true;
                break;
            }
        }
        if (!to_be_removed) {
            //std::cerr << " doNOT remove current_object2= " << current_object2 << std::endl;
            MDonly2.addObject(MD2.getObject(current_object2));
        }
    }
}

MetaDataTable MetaDataTable::combineMetaDataTables(std::vector<MetaDataTable> &MDin) {
    MetaDataTable MDc;

    if (MDin.size() == 0) {
        REPORT_ERROR("combineMetaDataTables ERROR: No input STAR files selected!");
    } else if (MDin.size() == 1) {
        MDc = MDin[0];
    } else {
        // Find which taTable combineMetaDataTables
        MetaDataTable commonLabels;

        // Loop over all labels in first
        // activeLabels is private but accessible from other instances of the same class in C++.
        for (size_t i = 0; i < MDin[0].activeLabels.size(); i++) {
            // Check their presence in each of the input files
            bool is_present = true;

            EMDL::EMDLabel thisLabel = MDin[0].activeLabels[i];
            std::string unknownLabel = "";
            if (thisLabel == EMDL::UNKNOWN_LABEL)
                unknownLabel = MDin[0].getUnknownLabelNameAt(i);

            for (size_t j = 1; j < MDin.size(); j++) {
                is_present = MDin[j].containsLabel(thisLabel, unknownLabel);

                if (!is_present) {
                    std::cerr << " + WARNING: ignoring label " << (unknownLabel.empty() ? EMDL::label2Str(thisLabel): unknownLabel) << " in " << j+1 << "th STAR file because it is not present in all STAR files to be combined." << std::endl;
                    break;
                }
            }

            if (is_present) {
                commonLabels.addLabel(thisLabel, unknownLabel);
            }
        }

        // Also disable any labels that do not occur in all input tables
        for (int i = 0; i < MDin.size(); i++) {
            for (int j = 0; j < MDin[i].activeLabels.size(); j++) {
                EMDL::EMDLabel thisLabel = MDin[i].activeLabels[j];
                std::string unknownLabel = "";
                if (thisLabel == EMDL::UNKNOWN_LABEL)
                    unknownLabel = MDin[i].getUnknownLabelNameAt(j);

                if (!commonLabels.containsLabel(thisLabel, unknownLabel)) {
                    MDin[i].deactivateLabel(thisLabel, unknownLabel);
                    std::cerr << " + WARNING: ignoring label " << (unknownLabel.empty() ? EMDL::label2Str(thisLabel) : unknownLabel) << " in " << i+1 << "th STAR file because it is not present in all STAR files to be combined." << std::endl;
                }
            }
        }

        // Then we can just append entire tables
        for (const MetaDataTable &mdt: MDin) { MDc.append(mdt); }
    }

    return MDc;
}

bool MetaDataTable::compareLabels(
    const MetaDataTable &MD1, const MetaDataTable &MD2
) {
    if (MD1.activeLabels.size() != MD2.activeLabels.size())
        return false;

    // Since we have the same number of labels, it suffices to check
    // all labels in MD1 are present in MD2.
    for (size_t id = 0; id < MD1.activeLabels.size(); id++) {
        EMDL::EMDLabel l = MD1.activeLabels[id];
        std::string unknownLabel = "";
        if (l == EMDL::UNKNOWN_LABEL)
            unknownLabel = MD1.getUnknownLabelNameAt(id);

        if (!MD2.containsLabel(l, unknownLabel))
            return false;
    }

    return true;
}

MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label, RFLOAT min_value, RFLOAT max_value
) {
    if (!EMDL::isInt(label) && !EMDL::isDouble(label))
        REPORT_ERROR("subsetMetadataTable ERROR: can only make a subset selection based on numbers");

    if (!MDin.containsLabel(label))
        REPORT_ERROR("subsetMetadataTable ERROR: input MetaDataTable does not contain label: " +  EMDL::label2Str(label));

    MetaDataTable MDout;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin) {
        RFLOAT val = EMDL::isInt(label) ? MDin.getValue<long>(label) : MDin.getValue<RFLOAT>(label);

        if (val <= max_value && val >= min_value) {
            MDout.addObject(MDin.getObject(index));
        }
    }

    return MDout;

}

MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label, std::string search_str, bool exclude
) {

    if (!EMDL::isString(label))
        REPORT_ERROR("subsetMetadataTable ERROR: can only make a subset selection based on strings");

    if (!MDin.containsLabel(label))
        REPORT_ERROR("subsetMetadataTable ERROR: input MetaDataTable does not contain label: " +  EMDL::label2Str(label));

    MetaDataTable MDout;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin) {
        if ((
            MDin.getValue<std::string>(label).find(search_str) != std::string::npos
        ) != exclude) {
            MDout.addObject(MDin.getObject(index));
        }
    }

    return MDout;
}

MetaDataTable removeDuplicatedParticles(
    MetaDataTable &MDin, EMDL::EMDLabel mic_label, RFLOAT threshold,
    RFLOAT origin_scale, FileName fn_removed, bool verb
) {
    // Sanity check
    if (!MDin.containsLabel(EMDL::ORIENT_ORIGIN_X_ANGSTROM) || !MDin.containsLabel(EMDL::ORIENT_ORIGIN_Y_ANGSTROM))
        REPORT_ERROR("You need rlnOriginXAngst and rlnOriginYAngst to remove duplicated particles");

    if (!MDin.containsLabel(EMDL::IMAGE_COORD_X) && !MDin.containsLabel(EMDL::IMAGE_COORD_Y))
        REPORT_ERROR("You need rlnCoordinateX, rlnCoordinateY to remove duplicated particles");

    if (!MDin.containsLabel(mic_label))
        REPORT_ERROR("STAR file does not contain " + EMDL::label2Str(mic_label));

    std::vector<bool> valid(MDin.numberOfObjects(), true);
    std::vector<RFLOAT> xs(MDin.numberOfObjects(), 0.0);
    std::vector<RFLOAT> ys(MDin.numberOfObjects(), 0.0);
    std::vector<RFLOAT> zs;

    bool dataIs3D = false;
    if (MDin.containsLabel(EMDL::IMAGE_COORD_Z)) {
        if (!MDin.containsLabel(EMDL::ORIENT_ORIGIN_Z_ANGSTROM))
            REPORT_ERROR("You need rlnOriginZAngst to remove duplicated 3D particles");
        dataIs3D = true;
        zs.resize(MDin.numberOfObjects(), 0.0);
    }

    // group by micrograph
    std::map<std::string, std::vector<long>> grouped;
    // Populate xs and ys
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin) {

        std::string mic_name = MDin.getValue<std::string>(mic_label);

        RFLOAT origin     = MDin.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM);
        RFLOAT coordinate = MDin.getValue<RFLOAT>(EMDL::IMAGE_COORD_X);
        xs[index] = coordinate - origin * origin_scale;

        origin     = MDin.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);
        coordinate = MDin.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y);
        ys[index] = coordinate - origin * origin_scale;

        if (dataIs3D) {
        origin     = MDin.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);
        coordinate = MDin.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z);
        zs[index] = coordinate - origin * origin_scale;
        }

        grouped[mic_name].push_back(index);
    }

    // The minimal permitted distance between any two particles
    // (technically the maximal forbidden distance)
    RFLOAT threshold_sq = threshold * threshold;

    // Remove duplicates
    // For each particle group
    for (
        std::map<std::string, std::vector<long>>::iterator it = grouped.begin();
        it != grouped.end(); ++it
    ) {

        // For every ordered pair of non-identical particles
        long n_particles = it->second.size();
        for (long i = 0; i < n_particles; i++) {
            long part_id1 = it->second[i];
            for (long j = i + 1; j < n_particles; j++) {
                long part_id2 = it->second[j];

                RFLOAT dx = xs[part_id1] - xs[part_id2];
                RFLOAT dy = ys[part_id1] - ys[part_id2];
                // The squared distance between the two particles
                RFLOAT dist_sq = dx * dx + dy * dy;
                if (dataIs3D) {
                    RFLOAT dz = zs[part_id1] - zs[part_id2];
                    dist_sq += dz * dz;
                }

                // If the particles are too close, invalidate one.
                if (dist_sq <= threshold_sq) {
                    // std::cout << it->first << " " << part_id1 << " " << part_id2 << " " << dist_sq << std::endl;
                    valid[part_id1] = false;
                    break;
                }
            }
        }
    }

    MetaDataTable MDout, MDremoved;
    long n_removed = 0;
    // (Bookkeeping) Make a note of which particles were kept and which removed.
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin) {
        if (valid[index]) {
            MDout.addObject(MDin.getObject(index));
        } else {
            MDremoved.addObject(MDin.getObject(index));
            n_removed++;
        }
    }

    if (fn_removed != "") MDremoved.write(fn_removed);

    std::cout << "Removed " << n_removed << " duplicated objects from " << MDin.numberOfObjects() << " objects." << std::endl;

    return MDout;
}
