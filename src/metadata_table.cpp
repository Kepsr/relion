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
    label_indices(EMDL::LAST_LABEL, -1),
    current_object(nullptr),
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
    label_indices(MD.label_indices),
    unknown_label_indices(MD.unknown_label_indices),
    unknownLabelNames(MD.unknownLabelNames),
    current_object(nullptr),
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
    if (this == &MD) return *this;

    clear();

    objects.resize(MD.objects.size());
    label_indices = MD.label_indices;
    unknown_label_indices = MD.unknown_label_indices;
    unknownLabelNames = MD.unknownLabelNames;
    current_object = &*objects.begin();
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
    }

    return *this;
}

MetaDataTable::~MetaDataTable() {
    for (const auto &object : objects) delete object;
}

bool MetaDataTable::isEmpty() const {
    return objects.empty();
}

size_t MetaDataTable::numberOfObjects() const {
    return objects.size();
}

void MetaDataTable::clear() {
    for (const auto &object : objects) delete object;
    objects.clear();

    label_indices = std::vector<long>(EMDL::LAST_LABEL, -1);
    current_object = &*objects.begin();
    unknown_label_indices.clear();
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

void MetaDataTable::setComment(const std::string &comment) {
    this->comment = comment;
}

std::string MetaDataTable::getComment() const {
    return comment;
}

bool MetaDataTable::containsComment() const {
    return !comment.empty();
}

void MetaDataTable::setName(const std::string &name) {
    this->name = name;
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

std::map<EMDL::EMDLabel, std::string> labels_and_unknowns(const MetaDataTable& mdt) {
    std::map<EMDL::EMDLabel, std::string> map;
    for (int i = 0; i < mdt.getActiveLabels().size(); i++)
        map.insert(mdt.label_and_unknown(i));
    return map;
}

std::string MetaDataTable::getUnknownLabelNameAt(int i) const {
    if (activeLabels[i] != EMDL::UNKNOWN_LABEL)
        REPORT_ERROR("MetaDataTable::getUnknownLabelNameAt(): the requested column is not an unknown label.");

    return unknownLabelNames[unknown_label_indices[i]];
}

std::string MetaDataTable::getValueToString(EMDL::EMDLabel label, long objectID) const {
    // SHWS 18 Jul 2018: this function previously had a stringstream, but it greatly slowed down
    // writing of large STAR files in some strange circumstances
    // (with large data.star and model.star files in refinement)
    // Therefore replaced the strstream with faster snprintf
    //
    // JZ 9 Aug 2018: still using a stringstream for vector<double> fields
    // => Avoid vector-valued columns in particle star-files.

    if (EMDL::is<std::string>(label)) return getValue<std::string>(label, objectID);

    if (EMDL::is<std::vector<double>>(label)) {
        const std::vector<double> v = getValue<std::vector<double>>(label, objectID);

        if (v.empty()) return "[]";

        std::stringstream sts;
        sts << std::setprecision(12);
        sts << '[';
        // It would be nice to use join for this
        // (but it doesn't currently work with string streams)
        for (int i = 0; i < v.size() - 1; i++) {
            sts << v[i] << ',';
        }
        sts << v[v.size() - 1];
        sts << ']';
        return sts.str();
    }

    char buffer[14];
    if (EMDL::is<double>(label)) {
        const double v = getValue<double>(label, objectID);
        auto mag = abs(v);
        if (mag > 0.0 && mag < 0.001 || mag > 100000.0) {
            // If the magnitude of v is very small or very large,
            // use floating-point form (6.02e-23).
            snprintf(buffer, 13, v < 0.0 ? "%12.5e" : "%12.6e", v);
        } else {
            // Otherwise, use fixed-point form (33.3333).
            snprintf(buffer, 13, v < 0.0 ? "%12.5f" : "%12.6f", v);
        }
    } else if (EMDL::is<int>(label)) {
        const long v = getValue<long>(label, objectID);
        snprintf(buffer, 13, "%12ld", v);
    } else if (EMDL::is<bool>(label)) {
        const bool v = getValue<bool>(label, objectID);
        snprintf(buffer, 13, "%12d", (int) v);
    }
    return std::string(buffer);
}

void MetaDataTable::setUnknownValue(int i, const std::string &value) {
    long j = unknown_label_indices[i];
    if (j < 0) REPORT_ERROR("MetaDataTable::setUnknownValue BUG: j should not be negative here....");
    (*current_object)->unknowns[j] = value;
}

void MetaDataTable::setValueFromString(
    EMDL::EMDLabel label, const std::string &value, long int objectID
) {
    if (EMDL::is<std::string>(label)) {
        setValue(label, value, objectID);
        return;
    } else {
        std::istringstream i(value);

        if (EMDL::is<double>(label)) {
            double v;
            i >> v;
            setValue(label, v, objectID);
            return;
        } else if (EMDL::is<int>(label)) {
            long v;
            i >> v;
            setValue(label, v, objectID);
            return;
        } else if (EMDL::is<bool>(label)) {
            bool v;
            i >> v;
            setValue(label, v, objectID);
            return;
        } else if (EMDL::is<std::vector<double>>(label)) {
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

            setValue(label, v, objectID);
            return;
        }
    }

    REPORT_ERROR("Logic error: should not happen");
}

// Comparators used for sorting (essentially lambdas)
namespace MD {

struct CompareAt {

    long i;

    CompareAt(long i): i(i) {}

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return true;  // Dummy implementation
    }

};

struct CompareDoublesAt: public CompareAt {

    using CompareAt::CompareAt;

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->doubles[i] < rh->doubles[i];
    }

};

struct CompareIntsAt: public CompareAt {

    using CompareAt::CompareAt;

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->ints[i] < rh->ints[i];
    }

};

struct CompareStringsAt: public CompareAt {

    using CompareAt::CompareAt;

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        return lh->strings[i] < rh->strings[i];
    }

};

struct CompareStringsAfterAtAt: public CompareAt {

    using CompareAt::CompareAt;

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        std::string slh = lh->strings[i];
        std::string srh = rh->strings[i];
        return slh.substr(slh.find("@") + 1) < srh.substr(srh.find("@") + 1);
    }

};

struct CompareStringsBeforeAtAt: public CompareAt {

    using CompareAt::CompareAt;

    bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const {
        std::string slh = lh->strings[i];
        std::string srh = rh->strings[i];
        std::stringstream stslh, stsrh;
        stslh << slh.substr(0, slh.find("@"));
        stsrh << srh.substr(0, srh.find("@"));
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
        srand(time(nullptr)); // initialise random seed
    } else if (!EMDL::is<int>(name) && !EMDL::is<double>(name)) {
        REPORT_ERROR("MetadataTable::sort%% " + prependERROR("can only sort numbers"));
    }

    std::vector<std::pair<double, long int>> vp;
    const auto N = objects.size();
    vp.reserve(N);
    for (long int i = 0; i < N; ++i) {
        if (do_random) {
            vp.emplace_back(rand(), i);
        } else if (EMDL::is<int>(name)) {
            vp.emplace_back(getValue<long>(name), i);
        } else {
            // EMDL::is<double>(name)
            vp.emplace_back(getValue<double>(name), i);
        }
    }

    std::sort(vp.begin(), vp.end());
    if (do_reverse && !do_random)
        std::reverse(vp.begin(), vp.end());

    if (only_set_index) {
        // Add an extra column with the sorted position of each entry
        for (long j = 0; j < vp.size(); j++) {
            setValue(EMDL::SORTED_IDX, j, vp[j].second);
        }
    } else {
        // Change the actual order in the MetaDataTable
        std::vector<MetaDataContainer*> objs;
        objs.reserve(N);

        for (const auto &x : vp) {
            objs.push_back(objects[x.second]);
        }

        objects = objs;
    }
    // reset pointer to the beginning of the table
    firstObject();
}

void MetaDataTable::newSort(const EMDL::EMDLabel label, bool do_sort_after_at, bool do_sort_before_at) {

    MD::CompareAt comp (0);  // Ideally, we wouldn't have to initialise the base class.

    if (EMDL::is<std::string>(label)) {
        if (do_sort_after_at) {
            comp = MD::CompareStringsAfterAtAt(label_indices[label]);
        } else if (do_sort_before_at) {
            comp = MD::CompareStringsBeforeAtAt(label_indices[label]);
        } else {
            comp = MD::CompareStringsAt(label_indices[label]);
        }
    } else if (EMDL::is<double>(label)) {
        comp = MD::CompareDoublesAt(label_indices[label]);
    } else if (EMDL::is<int>(label)) {
        comp = MD::CompareIntsAt(label_indices[label]);
    } else {
        REPORT_ERROR("Cannot sort this label: " + EMDL::label2Str(label));
    }

    std::stable_sort(objects.begin(), objects.end(), comp);

}

bool MetaDataTable::containsLabel(const EMDL::EMDLabel label, const std::string &unknownLabel) const {
    const bool is_known = label != EMDL::UNKNOWN_LABEL;
    const auto map = labels_and_unknowns(*this);
    return std::any_of(map.begin(), map.end(),
        [&] (const std::pair<EMDL::EMDLabel, std::string> &lu) {
        return label == lu.first && (is_known || unknownLabel == lu.second);
    });
}

std::vector<EMDL::EMDLabel> MetaDataTable::getActiveLabels() const {
    return activeLabels;
}

void MetaDataTable::deactivateLabel(EMDL::EMDLabel label, const std::string &unknownLabel) {
    const bool is_known = label != EMDL::UNKNOWN_LABEL;
    /// TODO: Get these erases out of this loop!
    for (int i = 0; i < activeLabels.size(); i++) {
        const auto lu = label_and_unknown(i);
        if (label == lu.first && (is_known || unknownLabel == lu.second)) {
            activeLabels.erase(activeLabels.begin() + i);
            unknown_label_indices.erase(unknown_label_indices.begin() + i);

            if (is_known)
                label_indices[label] = -1;  // This has to be a bug.
        }
    }
}

void MetaDataTable::addLabel(const EMDL::EMDLabel label, const std::string &unknownLabel) {
    if (label >= EMDL::LAST_LABEL)
        REPORT_ERROR(std::string(
            "MetaDataTable::addLabel: unrecognised label: "
        ) + EMDL::label2Str(label));
    const bool is_known = label != EMDL::UNKNOWN_LABEL;
    if (!is_known && unknownLabel.empty())
        REPORT_ERROR("MetaDataTable::addLabel: unknownLabel is empty");

    if (label_indices[label] >= 0 && is_known) return;
    // keep pushing the same unknown label...
    long i;

    if (EMDL::is<double>(label)) {
        for (const auto &object : objects) {
            object->doubles.push_back(0);
        }
        i = doubleLabels++;
    } else if (EMDL::is<int>(label)) {
        for (const auto &object : objects) {
            object->ints.push_back(0);
        }
        i = intLabels++;
    } else if (EMDL::is<bool>(label)) {
        for (const auto &object : objects) {
            object->bools.push_back(false);
        }
        i = boolLabels++;
    } else if (EMDL::is<std::string>(label)) {
        for (const auto &object : objects) {
            object->strings.emplace_back("empty");
        }
        i = stringLabels++;
    } else if (EMDL::is<std::vector<double>>(label)) {
        for (const auto &object : objects) {
            object->doubleVectors.emplace_back();
        }
        i = doubleVectorLabels++;
    } else if (EMDL::is<void>(label)) {
        for (const auto &object : objects) {
            object->unknowns.emplace_back("empty");
        }
        unknownLabelNames.push_back(unknownLabel);
        i = unknownLabels++;
    }

    activeLabels.push_back(label);
    unknown_label_indices.push_back(EMDL::is<void>(label) ? i : -1);

    label_indices[label] = i;
}

void MetaDataTable::addMissingLabels(const MetaDataTable &mdt) {
    for (const auto &lu : labels_and_unknowns(mdt)) {
        if (lu.first == EMDL::UNKNOWN_LABEL && !containsLabel(lu.first, lu.second) ||
            label_indices[lu.first] < 0) {
            addLabel(lu.first, lu.second);
        }
    }
}

void MetaDataTable::append(const MetaDataTable &mdt) {
    if (activeLabels.empty()) {
        // If the current one is empty, add missing labels and append the new one:
        addMissingLabels(mdt);
    } else {
        // If the current one is not-empty, check all labels are the same before appending. Otherwise, raise error
        if (!compareLabels(*this, mdt))
            REPORT_ERROR("ERROR in appending metadata tables with not the same columns!");
    }

    // Now append
    objects.reserve(objects.size() + mdt.objects.size());
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
    try {
        if (objectID < 0) { objectID = current_object - &*objects.begin(); }
        checkObjectID(objectID);
        return objects[objectID];
    } catch (const std::string &errmsg) {
        REPORT_ERROR((std::string) __func__ + ": " + errmsg);
    }
}

void MetaDataTable::setObject(MetaDataContainer* data, long objectID) {
    try { 
        if (objectID < 0) { objectID = current_object - &*objects.begin(); }
        checkObjectID(objectID); 
        addMissingLabels(*data->table);
        setObjectUnsafe(data, objectID);
    } catch (const std::string &errmsg) {
        REPORT_ERROR((std::string) __func__ + ": " + errmsg);
    }
}

void MetaDataTable::setValuesOfDefinedLabels(MetaDataContainer* data, long objectID) {
    try {
        if (objectID < 0) { objectID = current_object - &*objects.begin(); }
        checkObjectID(objectID);
        setObjectUnsafe(data, objectID);
    } catch (const std::string &errmsg) {
        REPORT_ERROR((std::string) __func__ + ": " + errmsg);
    }
}

void MetaDataTable::reserve(size_t capacity) {
    objects.reserve(capacity);
}

void MetaDataTable::setObjectUnsafe(MetaDataContainer* data, long objectID) {
    MetaDataContainer* obj = objects[objectID];

    for (long i = 0; i < data->table->activeLabels.size(); i++) {
        EMDL::EMDLabel label = data->table->activeLabels[i];

        if (label != EMDL::UNKNOWN_LABEL) {
            const long this_off =              label_indices[label];
            const long that_off = data->table->label_indices[label];

            if (this_off < 0) continue;

            if (EMDL::is<double>(label)) {
                obj->doubles[this_off] = data->doubles[that_off];
            } else if (EMDL::is<int>(label)) {
                obj->ints[this_off] = data->ints[that_off];
            } else if (EMDL::is<bool>(label)) {
                obj->bools[this_off] = data->bools[that_off];
            } else if (EMDL::is<std::string>(label)) {
                obj->strings[this_off] = data->strings[that_off];
            } else if (EMDL::is<std::vector<double>>(label)) {
                obj->doubleVectors[this_off] = data->doubleVectors[that_off];
            }
        } else {
            const long that_off            = data->table->unknown_label_indices[i];
            const std::string unknownLabel = data->table->unknownLabelNames[that_off];
            const auto search = std::find(unknownLabelNames.begin(), unknownLabelNames.end(), unknownLabel);

            if (search == unknownLabelNames.end())
                REPORT_ERROR("MetaDataTable::setObjectUnsafe: logic error."
                             "unknownLabel was not found.");

            obj->unknowns[search - unknownLabelNames.begin()] = data->unknowns[that_off];
        }
    }
}

void MetaDataTable::addObject() {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels, stringLabels,
        doubleVectorLabels, unknownLabels
    ));

    current_object = &*objects.end() - 1;
}

void MetaDataTable::addObject(MetaDataContainer* data) {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels, stringLabels,
        doubleVectorLabels, unknownLabels
    ));

    setObject(data, objects.size() - 1);
    current_object = &*objects.end() - 1;
}

void MetaDataTable::addValuesOfDefinedLabels(MetaDataContainer* data) {
    objects.push_back(new MetaDataContainer(
        this,
        doubleLabels, intLabels, boolLabels,
        stringLabels, doubleVectorLabels, unknownLabels
    ));

    setValuesOfDefinedLabels(data, objects.size() - 1);
    current_object = &*objects.end() - 1;
}

void MetaDataTable::removeObject(long objectID) {
    try {
        if (objectID < 0) { objectID = current_object - &*objects.begin(); }
        checkObjectID(objectID);
        delete objects[objectID];
        objects.erase(objects.begin() + objectID);
        current_object = &*objects.end() - 1;
    } catch (const std::string &errmsg) {
        REPORT_ERROR((std::string) __func__ + ": " + errmsg);
    }
}

MetaDataContainer** MetaDataTable::firstObject() {
    return current_object = &*objects.begin();
}

MetaDataContainer** MetaDataTable::nextObject() {
    return ++current_object;
}

MetaDataContainer** MetaDataTable::goToObject(long int objectID) {
    try {
        checkObjectID(objectID);
        return current_object = &*objects.begin() + objectID;
    } catch (const std::string &errmsg) {
        REPORT_ERROR((std::string) __func__ + ": " + errmsg);
    }
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
    for (auto _ : *this) {
        // As the internal state of MDin gets incremented...
        RFLOAT val;
        if (EMDL::is<double>(label)) {
            val = getValue<double>(label);
        } else if (EMDL::is<int>(label)) {
            val = getValue<long>(label);
        } else if (EMDL::is<bool>(label)) {
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
    for (RFLOAT value : values) {
        sum   += value;
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
            bin_size = (unsigned int) (std::ceil((hist_max - hist_min) / bin_width));
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
            histX[4 * i + 0] = hist_min - bin_width;
            histX[4 * i + 1] = hist_min - bin_width;
            histX[4 * i + 2] = hist_min;
            histX[4 * i + 3] = hist_min;
        } else if (i == bin_size - 1) {
            if (verb > 0) std::cout << "[" << hist_max << ", +INF]: ";
            histX[4 * i + 0] = hist_max;
            histX[4 * i + 1] = hist_max;
            histX[4 * i + 2] = hist_max + bin_width;
            histX[4 * i + 3] = hist_max + bin_width;
        } else {
            if (verb > 0) std::cout << "[" << (hist_min + bin_width * (i - 1)) << ", " << (hist_min + bin_width * i) << "): ";
            histX[4 * i + 0] = hist_min + bin_width * (i - 1);
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
        histY[4 * i + 0] = histY[4 * i + 3] = 0.0;

    }
    histX[histX.size() - 1] = histX[histX.size() - 2];

    if (plot2D) {
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
        const long offx = label_indices[xaxis];
        if (offx < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find x-axis label");

        if (xaxis == EMDL::UNDEFINED) {
            xval = idx + 1;
        } else if (EMDL::is<double>(xaxis)) {
            xval = objects[idx]->getValue<double>(offx);
        } else if (EMDL::is<int>(xaxis)) {
            xval = objects[idx]->getValue<int>(offx);
        } else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot x-axis double, int or long int");

        const long offy = label_indices[yaxis];
        if (offy < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find y-axis label");

        if (EMDL::is<double>(yaxis)) {
            yval = objects[idx]->getValue<double>(offy);
        } else if (EMDL::is<int>(yaxis)) {
            yval = objects[idx]->getValue<int>(offy);
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
    for (EMDL::EMDLabel label : activeLabels)
        ost << EMDL::label2Str(label) << "\n";
}

void MetaDataTable::randomiseOrder() {
    std::random_shuffle(objects.begin(), objects.end());
}

// Bounds checking
void MetaDataTable::checkObjectID(long i) const throw (std::string) {
    if (i < 0 || i >= objects.size())
        throw "object " + std::to_string(i) + " out of bounds! (" + std::to_string(objects.size()) + " objects present)";
}

//FIXME: does not support unknownLabels but this function is only used by relion_star_handler
//       so I will leave this for future...
// Decompose two tables (A, B) into (A - B, A & B, B - A)
void compareMetaDataTable(
    MetaDataTable &MD1, MetaDataTable &MD2,
    MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
    EMDL::EMDLabel label1, double eps, EMDL::EMDLabel label2, EMDL::EMDLabel label3
) {
    if (!MD1.containsLabel(label1))
        REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label1.");
    if (!MD2.containsLabel(label1))
        REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label1.");

    if (label2 != EMDL::UNDEFINED) {
        if (!EMDL::is<double>(label1) || !EMDL::is<double>(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 2D or 3D distances are only allowed for doubles.");
        if (!MD1.containsLabel(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label2.");
        if (!MD2.containsLabel(label2))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label2.");
    }

    if (label3 != EMDL::UNDEFINED) {
        if (!EMDL::is<double>(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 3D distances are only allowed for doubles.");
        if (!MD1.containsLabel(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label3.");
        if (!MD2.containsLabel(label3))
            REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label3.");
    }

    MDboth.clear();
    MDonly1.clear();
    MDonly2.clear();

    std::string mystr1;
    long int myint1;
    double myd1, mydy1 = 0.0, mydz1 = 0.0;

    // loop over MD1
    std::vector<long int> to_remove_from_only2;
    for (long int current_object1 : MD1) {
        if (EMDL::is<std::string>(label1)) {
            mystr1 = MD1.getValue<std::string>(label1);
        } else if (EMDL::is<int>(label1)) {
            myint1 = MD1.getValue<int>(label1);
        } else if (EMDL::is<double>(label1)) {
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
        for (long int current_object2 : MD2) {
            if (EMDL::is<std::string>(label1)) {
                std::string mystr2 = MD2.getValue<std::string>(label1);
                if (mystr1 == mystr2) {
                    have_in_2 = true;
                    to_remove_from_only2.push_back(current_object2);
                    MDboth.addObject(MD1.getObject());
                    break;
                }
            } else if (EMDL::is<int>(label1)) {
                long int myint2 = MD2.getValue<int>(label1);
                if (abs(myint2 - myint1) <= round(eps)) {
                    have_in_2 = true;
                    to_remove_from_only2.push_back(current_object2);
                    MDboth.addObject(MD1.getObject());
                    break;
                }
            } else if (EMDL::is<double>(label1)) {
                double myd2 = MD2.getValue<double>(label1);
                double mydy2 = label2 == EMDL::UNDEFINED ? 0.0 :
                    MD2.getValue<double>(label2);
                double mydz2 = label3 == EMDL::UNDEFINED ? 0.0 :
                    MD2.getValue<double>(label3);

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

    for (long int current_object2 : MD2) {
        // If there is no current_object2 in to_remove_from_only2
        if (std::find(
            to_remove_from_only2.begin(), to_remove_from_only2.end(), current_object2
        ) == to_remove_from_only2.end()) {
            // std::cerr << " doNOT remove current_object2= " << current_object2 << std::endl;
            MDonly2.addObject(MD2.getObject(current_object2));
        }
    }
}

MetaDataTable MetaDataTable::getCommonLabels(const std::vector<MetaDataTable> &mdts) {
    MetaDataTable common_labels;
    // Assume mdts.size() > 1
    const MetaDataTable &first_mdt = mdts[0];
    // For each label in the first table
    for (const auto &lu : labels_and_unknowns(first_mdt)) {
        // Is this label present in all of the other MetaDataTables?
        for (size_t j = 1; j < mdts.size(); j++) {
            if (!mdts[j].containsLabel(lu.first, lu.second)) {
                std::cerr << " + WARNING: ignoring label " << (lu.second.empty() ? EMDL::label2Str(lu.first) : lu.second)
                << " in " << j + 1 << "th STAR file because it is not present in all STAR files to be combined." << std::endl;
                goto next_label;
            }
        }
        common_labels.addLabel(lu.first, lu.second);
        next_label: {}
    }
    return common_labels;
}

MetaDataTable MetaDataTable::combineMetaDataTables(std::vector<MetaDataTable> &MDin) {

    if (MDin.empty())
        REPORT_ERROR("combineMetaDataTables ERROR: No input STAR files selected!");
    if (MDin.size() == 1)
        return MDin[0];

    MetaDataTable MDc;
    // Find which taTable combineMetaDataTables

    MetaDataTable common_labels = getCommonLabels(MDin);

    /// TODO: Can the loop used to construct common_labels also be used to disable labels?

    // And disable those labels that do not occur in all input tables
    for (int i = 0; i < MDin.size(); i++) {
        MetaDataTable &mdt = MDin[i];
        for (const auto &lu : labels_and_unknowns(mdt)) {
            if (!common_labels.containsLabel(lu.first, lu.second)) {
                mdt.deactivateLabel(lu.first, lu.second);
                std::cerr << " + WARNING: ignoring label " << (lu.second.empty() ? EMDL::label2Str(lu.first) : lu.second) << " in " << i + 1 << "th STAR file"
                "because it is not present in all STAR files to be combined." << std::endl;
            }
        }
    }

    // Then we can just append entire tables
    for (const MetaDataTable &mdt: MDin) { MDc.append(mdt); }
}

bool MetaDataTable::compareLabels(
    const MetaDataTable &MD1, const MetaDataTable &MD2
) {
    if (MD1.activeLabels.size() != MD2.activeLabels.size())
        return false;

    // Since we have the same number of labels,
    // it suffices to check that all labels in MD1 are present in MD2.
    const auto map = labels_and_unknowns(MD1);
    return std::all_of(map.begin(), map.end(),
        [&] (const std::pair<EMDL::EMDLabel, std::string> &lu) {
        return MD2.containsLabel(lu.first, lu.second);
    });
}

MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label, RFLOAT min_value, RFLOAT max_value
) {
    if (!EMDL::is<int>(label) && !EMDL::is<double>(label))
        REPORT_ERROR("subsetMetadataTable ERROR: can only make a subset selection based on numbers");

    if (!MDin.containsLabel(label))
        REPORT_ERROR("subsetMetadataTable ERROR: input MetaDataTable does not contain label: " + EMDL::label2Str(label));

    MetaDataTable MDout;
    for (long int i : MDin) {
        RFLOAT x = EMDL::is<int>(label) ? MDin.getValue<long>(label) :
            MDin.getValue<RFLOAT>(label);

        if (x <= max_value && x >= min_value)
            MDout.addObject(MDin.getObject(i));
    }
    return MDout;
}

MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label, const std::string &search_str, bool exclude
    // exclude determines whether to perform A * +B or A * -B
) {

    if (!EMDL::is<std::string>(label))
        REPORT_ERROR("subsetMetadataTable ERROR: can only make a subset selection based on strings");

    if (!MDin.containsLabel(label))
        REPORT_ERROR("subsetMetadataTable ERROR: input MetaDataTable does not contain label: " + EMDL::label2Str(label));

    MetaDataTable MDout;
    for (long int i : MDin) {
        if (exclude == (
            MDin.getValue<std::string>(label).find(search_str) == std::string::npos
        )) MDout.addObject(MDin.getObject(i));
    }
    return MDout;
}

// Map each micrograph name to a collection of object indices
// Also, populate xs, ys, zs
static std::map<std::string, std::vector<long>> group_particles_by_micrograph(
    MetaDataTable mdt, EMDL::EMDLabel mic_label, RFLOAT origin_scale,
    std::vector<RFLOAT> &xs, std::vector<RFLOAT> &ys, std::vector<RFLOAT> &zs,
    bool dataIs3D
) {
    std::map<std::string, std::vector<long>> grouped;
    for (long int i : mdt) {

        RFLOAT origin, coordinate;

        origin     = mdt.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM);
        coordinate = mdt.getValue<RFLOAT>(EMDL::IMAGE_COORD_X);
        xs[i] = coordinate - origin * origin_scale;

        origin     = mdt.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);
        coordinate = mdt.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y);
        ys[i] = coordinate - origin * origin_scale;

        if (dataIs3D) {
        origin     = mdt.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);
        coordinate = mdt.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z);
        zs[i] = coordinate - origin * origin_scale;
        }

        std::string mic_name = mdt.getValue<std::string>(mic_label);
        grouped[mic_name].push_back(i);  // Will insert key-value pair if key not found

    }
    return grouped;
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

    std::vector<RFLOAT> xs (MDin.numberOfObjects(), 0.0);
    std::vector<RFLOAT> ys (MDin.numberOfObjects(), 0.0);
    std::vector<RFLOAT> zs;
    bool dataIs3D = false;
    if (MDin.containsLabel(EMDL::IMAGE_COORD_Z)) {
        if (!MDin.containsLabel(EMDL::ORIENT_ORIGIN_Z_ANGSTROM))
            REPORT_ERROR("You need rlnOriginZAngst to remove duplicated 3D particles");
        dataIs3D = true;
        zs.resize(MDin.numberOfObjects(), 0.0);
    }

    std::map<std::string, std::vector<long>> grouped = group_particles_by_micrograph(
        MDin, mic_label, origin_scale, xs, ys, zs, dataIs3D
    );

    // The minimal permitted distance between any two particles
    // (technically the maximal forbidden distance)
    RFLOAT threshold_sq = threshold * threshold;

    // For each particle group, remove duplicates
    std::vector<bool> valid (MDin.numberOfObjects(), true);
    for (auto mic_name_and_object_indices : grouped) {

        // For every ordered pair of non-identical particles
        long n_particles = mic_name_and_object_indices.second.size();
        for (long i = 0; i < n_particles; i++) {
        long part_id1 = mic_name_and_object_indices.second[i];
        for (long j = i + 1; j < n_particles; j++) {
        long part_id2 = mic_name_and_object_indices.second[j];

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
                // std::cout << mic_name_and_object_indices.first << " " << part_id1 << " " << part_id2 << " " << dist_sq << std::endl;
                valid[part_id1] = false;
                break;
            }
        }
        }
    }

    MetaDataTable MDout, MDremoved;
    // (Bookkeeping) Make a note of which particles were kept and which removed.
    for (long int i : MDin) {
        (valid[i] ? MDout : MDremoved).addObject(MDin.getObject(i));
    }

    if (!fn_removed.empty()) MDremoved.write(fn_removed);

    std::cout << "Removed " << MDremoved.numberOfObjects() << " duplicated objects from " << MDin.numberOfObjects() << " objects." << std::endl;

    return MDout;
}

/// STAR file i/o

// Reading

long int MetaDataTable::read(
    const FileName &filename, const std::string &name, bool do_only_count
) {

    clear();  // Clear current table

    FileName fn_read = filename.removeFileFormat();  // Check for a :star extension

    std::ifstream in (fn_read.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR((std::string) "MetaDataTable::read: File " + fn_read + " does not exist");

    return readStar(in, name, do_only_count);

    firstObject();  // Go to the first object
}

long int MetaDataTable::readStarLoop(std::ifstream &in, bool do_only_count) {
    isList = false;

    // Read column labels
    int labelPosition = 0;
    std::string line, token;

    // First read all the column labels
    while (getline(in, line, '\n')) {
        line = simplify(line);
        // TODO: handle comments...
        if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
            continue;

        if (line[0] != '_') break;  // found first data line

        // label definition line
        // Take string from "_" until "#"
        size_t start = line.find("_");
        size_t end   = line.find("#");

        token = line.substr(start + 1, end - start - 2);

        EMDL::EMDLabel label = EMDL::str2Label(token);

        if (label == EMDL::UNDEFINED) {
            std::cerr << " + WARNING: will ignore (but maintain) values for the unknown label: " << token << std::endl;
            label = EMDL::UNKNOWN_LABEL;
        }

        addLabel(label, token);

        labelPosition++;
    }

    // Then fill the table (dont read another line until the one from above has been handled)
    long int nr_objects = 0;
    const int num_labels = activeLabels.size();
    do {

        line = simplify(line);
        // Stop at empty line
        if (line[0] == '\0') break;

        nr_objects++;
        if (!do_only_count) {

            addObject();  // Add a new line to the table

            // Parse data values
            int pos = 0;
            std::string value;
            labelPosition = 0;
            for (; nextTokenInSTAR(line, pos, value); labelPosition++) {
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
            }
            if (labelPosition < num_labels && num_labels > 2) {
                // For backward-compatibility for cases like "fn_mtf <empty>", don't die if num_labels == 2.
                std::cerr << "Error in line: " << line << std::endl;
                REPORT_ERROR("A line in the STAR file contains fewer columns than the number of labels. Expected = " + integerToString(num_labels) + " Found = " +  integerToString(labelPosition));
            }
        }
    } while (getline(in, line, '\n'));

    return nr_objects;
}

bool MetaDataTable::readStarList(std::ifstream &in) {
    isList = true;
    addObject();
    long int objectID = objects.size() - 1;

    std::string line, firstword, value;

    // Read data and fill structures accordingly
    int labelPosition = 0;
    while (getline(in, line, '\n')) {
        int pos = 0;
        // Ignore empty lines
        if (!nextTokenInSTAR(line, pos, firstword)) continue;

        // Get label-value pairs
        if (firstword[0] == '_') {
            std::string token = firstword.substr(1); // get rid of leading underscore
            EMDL::EMDLabel label = EMDL::str2Label(token);
            if (!nextTokenInSTAR(line, pos, value))
                REPORT_ERROR("MetaDataTable::readStarList: did not encounter a single word after " + firstword);

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
            return true;
        } else if (firstword.find("data_") == 0) {
            // Check whether this data blocks ends (because a next one is there)
            // Should I reverse the pointer one line?
            return false;
        }
    }
    // Reached the end of the file
    return false;
}

long int MetaDataTable::readStar(
    std::ifstream &in, const std::string &name, bool do_only_count
) {
    std::string line, token, value;
    clear();

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
                        return !readStarList(in);
                    }
                }
            }
        }
    }

    // Clear the eofbit so we can perform more actions on the stream.
    in.clear();

    return 0;
}

// Writing

void MetaDataTable::write(std::ostream& out) {
    // Only write tables that have something in them
    if (isEmpty()) return;

    if (version >= 30000) {
        out << "\n"
            << "# version " << getCurrentVersion() << "\n";
    }

    out << "\n"
        << "data_" << getName() << "\n";

    if (containsComment())
    out << "# " << comment << "\n";

    out << "\n";

    if (!isList) {

        // Write loop header
        out << "loop_\n";

        long int i = 0;
        for (const auto &lu : labels_and_unknowns(*this)) {
            if (lu.first != EMDL::COMMENT && lu.first != EMDL::SORTED_IDX) {
                // EMDL::SORTED_IDX is only for internal use, never write it out!
                out << "_" << (lu.first == EMDL::UNKNOWN_LABEL ? lu.second :
                    EMDL::label2Str(lu.first)) << " #" << ++i << "\n";
            }
        }

        // Write data block
        for (long int j = 0; j < objects.size(); j++) {
            std::string entryComment = "";

            for (long int i = 0; i < activeLabels.size(); i++) {
                EMDL::EMDLabel l = activeLabels[i];
                if (l == EMDL::SORTED_IDX) continue;
                if (l == EMDL::COMMENT) {
                    entryComment = "# " + getValue<std::string>(EMDL::COMMENT, j);
                } else {
                    out.width(10);
                    out << escapeStringForSTAR(l == EMDL::UNKNOWN_LABEL ?
                        objects[j]->unknowns[unknown_label_indices[i]] :
                        getValueToString(l, j))
                        << " ";
                }
            }
            out << entryComment << "\n";
        }

    } else {
        // isList
        // Get first object. In this case (row format) there is a single object
        std::string entryComment = "";

        // Determine column width
        int maxWidth = 10;
        for (const auto &lu : labels_and_unknowns(*this)) {
            if (lu.first == EMDL::COMMENT) {
                entryComment = "# " + getValue<std::string>(EMDL::COMMENT, 0);
            } else {
                int w = (lu.first == EMDL::UNKNOWN_LABEL ? lu.second : 
                    EMDL::label2Str(lu.first)).length();
                if (w > maxWidth) maxWidth = w;
            }
        }

        for (long i = 0; i < activeLabels.size(); i++) {
            const auto lu = label_and_unknown(i);
            if (lu.first == EMDL::COMMENT) continue;
            const bool is_known = lu.first == EMDL::UNKNOWN_LABEL;
            const std::string key = is_known ? EMDL::label2Str(lu.first) : lu.second;
            const std::string val = escapeStringForSTAR(is_known ?
                getValueToString(lu.first, 0) :
                objects[0]->unknowns[unknown_label_indices[i]]);
            out << "_" << key << std::setw(12 + maxWidth - key.length()) << " " << val << "\n";
        }
        out << entryComment << "\n";
    }
    // Finish table/data block with an empty line
    out << "\n";
}

void MetaDataTable::write(const FileName &fn_out) {
    FileName fn_tmp = fn_out + ".tmp";
    std::ofstream fh (fn_tmp.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR((std::string) "MetaDataTable::write: cannot write to file: " + fn_out);
        // fh << "# RELION; version " << g_RELION_VERSION << std::endl;
    write(fh);
    // Rename to prevent errors with programs in pipeliner reading in incomplete STAR files
    std::rename(fn_tmp.c_str(), fn_out.c_str());
}
