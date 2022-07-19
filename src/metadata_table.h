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
 * Redesigned by: J. Zivanov in June 2017
 * MRC Laboratory of Molecular Biology
 *
 * Original author: J.R. Bilbao-Castro (jrbcast@ace.ual.es)
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
 ***************************************************************************/

#ifndef METADATA_TABLE_H
#define METADATA_TABLE_H

#include <map>
#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <sstream>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include "src/funcs.h"
#include "src/args.h"
#include "src/CPlot2D.h"
#include "src/metadata_container.h"
#include "src/metadata_label.h"

const int CURRENT_MDT_VERSION = 30001;

/*	class MetaDataTable:
 *
 *	- stores a table of values for an arbitrary subset of predefined EMDLabels
 *	- each column corresponds to a label
 *	- each row represents a data point
 *	- the rows are stored in per-type contiguous blocks of memory
 *
 *	2020/Nov/12:
 *	  This class is organized as an array (`objects`) of structures (`MetaDataContainer`).
 *
 *        `activeLabels` contains all valid labels.
 *        Even when a label is `deactivateLabel`-ed, the values remain in `MetaDataContainer`s.
 *        The label is only removed from `activeLabels`.
 *
 *        Each data type (int, double, etc) contains its own storage array inside `MetaDataContainer`.
 *        Thus, values in `label2offsets` are NOT unique. Accessing columns via a wrong type is
 *        very DANGEROUS. Use `cmake -DMDT_TYPE_CHECK=ON` to enable runtime checks.
 *
 *        Handling of labels unknown to RELION needs care.
 *        They all share the same label, EMD_UNKNOWN_LABEL. Thus, `addLabel`, `containsLabel`,
 *        `compareLabels` etc must check not only EMDLabel in `activeLabels` but also the
 *        real labels stored in `unknownLabelNames`. This should be done via `getUnknownLabelNameAt`.
 *	  Note that two STAR files might contain the same set of unknown labels, but in different orders.
 *
 *        Whenever `activeLabels` is modified, `unknownLabelPosition2Offset` MUST be updated accordingly.
 *        When the label for a column is EMD_UNKNOWN_LABEL, the corresponding element in
 *        `unknownLabelPosition2Offset` must store the offset in `unknownLabelNames` and
 *        `MetaDataContainer->unknowns`. Otherwise, the value does not matter.
 */
class MetaDataTable {

    // Effectively stores all metadata
    std::vector<MetaDataContainer*> objects;

    // Maps labels to corresponding indices in the vectors in MetaDataContainer.
    // The length of label2offset is always equal to the number of defined labels (~320)
    // e.g.:
    // the value of "defocus-U" for row r is stored in:
    //	 objects[r]->doubles[label2offset[EMDL::CTF_DEFOCUSU]]
    // the value of "image name" is stored in:
    //	 objects[r]->strings[label2offset[EMDL::IMAGE_NAME]]
    std::vector<long> label2offset;

    /** What labels have been read from a docfile/metadata file
     *  and/or will be stored on a new metadata file when "save" is
     *  called
     **/
    std::vector<EMDL::EMDLabel> activeLabels;

    std::vector<std::string> unknownLabelNames;
    std::vector<long> unknownLabelPosition2Offset;

    // Current object id
    long current_objectID;

    // Number of labels of each type
    long doubleLabels, intLabels, boolLabels, stringLabels, doubleVectorLabels, unknownLabels;

    // Name of the metadata table
    std::string name;

    // A comment for the metadata table
    std::string comment;

    // The version number of the file format (multiplied by 10,000)
    int version;

    public:

    MetaDataTable();

    // Copy constructor
    // Fill the new table with *copies* of all objects
    MetaDataTable(const MetaDataTable &c);

    // Assignment
    // Fill the new table with *copies* of all objects
    MetaDataTable& operator = (const MetaDataTable &MD);

    ~MetaDataTable();

    // Is this a 1D list (as opposed to a 2D table)?
    bool isList;

    bool isEmpty() const;
    size_t numberOfObjects() const;
    void clear();

    void setComment(const std::string Comment);
    std::string getComment() const;
    bool containsComment() const;

    void setName(const std::string Name);
    std::string getName() const;

    void setVersion(int v);
    int getVersion() const;
    static int getCurrentVersion();

    template<class T>
    T getValue(EMDL::EMDLabel label, long objectID = -1) const;

    std::string getValueToString(EMDL::EMDLabel label, long int objectID = -1, bool escape=false) const;

    std::string getUnknownLabelNameAt(int i) const;

    // Set the value of label for a specified object.
    // If no objectID is given, the internal iterator 'current_objectID' is used
    // objectID is 0-indexed.
    template<class T>
    bool setValue(EMDL::EMDLabel name, const T &value, long int objectID = -1);

    void setUnknownValue(int labelPosition, const std::string &value);
    bool setValueFromString(EMDL::EMDLabel label, const std::string &value, long int objectID = -1);

    // Sort the order of the elements based on the values in the input label
    // (only numbers, no strings/bools)
    void sort(EMDL::EMDLabel name, bool do_reverse = false, bool only_set_index = false, bool do_random = false);

    void newSort(const EMDL::EMDLabel name, bool do_reverse = false, bool do_sort_after_at = false, bool do_sort_before_at = false);

    // Check whether a label is defined in the table.
    // This is redundant and will be removed in 3.2.
    bool labelExists(EMDL::EMDLabel name) const;

    // Does 'activeLabels' contain 'label'?
    bool containsLabel(const EMDL::EMDLabel label, const std::string unknownLabel="") const;

    std::vector<EMDL::EMDLabel> getActiveLabels() const;

    // Deactivate a column from a table, so that it is no longer written out
    void deactivateLabel(EMDL::EMDLabel label, std::string unknownLabel="");

    // Add a new label and update all objects
    void addLabel(EMDL::EMDLabel label, std::string unknownLabel="");

    // Add missing labels that are present in 'app'.
    void addMissingLabels(const MetaDataTable *app);

    // Append all rows from 'app' to the end of the table and insert all missing labels.
    void append(const MetaDataTable &app);

    // Get metadatacontainer for objectID (current_objectID if objectID < 0)
    MetaDataContainer* getObject(long objectID = -1) const;

    /* setObject(data, objectID)
     *  copies values from 'data' to object 'objectID'.
     *  The target object is assumed to exist.
     *  If objectID < 0, then current_objectID is set.
     *  Undefined labels are inserted.
     *
     *  Use addObject() to set an object that does not yet exist */
    void setObject(MetaDataContainer* data, long objectID = -1);

    /* setValuesOfDefinedLabels(data, objectID)
     * copies values from 'data' to object 'objectID'.
     * The target object is assumed to exist.
     * If objectID < 0, then current_objectID is set.
     * Only already defined labels are considered.
     *
     * Use addValuesOfDefinedLabels() to add an object that does not yet exist */
    void setValuesOfDefinedLabels(MetaDataContainer* data, long objectID = -1);

    // reserve memory for this many lines
    void reserve(size_t capacity);

    /* addObject()
     *  Adds a new object and initializes the defined labels with default values.
     *  Afterwards, 'current_objectID' points to the newly added object.*/
    void addObject();

    /* addObject(data)
     *  Adds a new object and sets its values to those from 'data'.
     *  The set of labels for the table is extended as necessary.
     *  Afterwards, 'current_objectID' points to the newly added object.*/
    void addObject(MetaDataContainer* data);

    /* addValuesOfDefinedLabels(data)
     *  Adds a new object and sets the already defined values to those from 'data'.
     *  Labels from 'data' that are not already defined are ignored.
     *  Afterwards, 'current_objectID' points to the newly added object.*/
    void addValuesOfDefinedLabels(MetaDataContainer* data);

    /* removeObject(objectID)
     *  If objectID is not given, 'current_objectID' will be removed.
     *  'current_objectID' is set to the last object in the list. */
    void removeObject(long objectID = -1);

    long firstObject();
    long nextObject();

    /** @TODO: remove nextObject() after removing calls in:
     * - "helix.cpp"
     */
    
    /** MetaDataTable::iterator
     *    
     * This struct lets us iterate over the object indices in a MetaDataTable:
     * @code
     * for (long int i : mdt) {
     *     foo(i);
     * }
     * @endcode
     *
     * This is not thread-safe, because current_objectID is updated.
     */
    struct iterator {

        bool isDone;
        MetaDataTable *const mdt;
        long int i;

        iterator(MetaDataTable *mdt, bool isDone = false):
        isDone(isDone), mdt(mdt), i(isDone ? 0 : mdt->firstObject()) {}

        long int operator *() const {
            return i;
        }

        iterator &operator ++() {
            // assert !isDone
            if (!(isDone = i + 1 >= mdt->numberOfObjects())) { i = mdt->nextObject(); }
            return *this;
        }

        bool operator != (const iterator &other) const {
            return isDone != other.isDone || mdt != other.mdt || i != other.i;
        }

    };

    iterator begin() {
        return iterator(this, false);
    }

    iterator end() {
        return iterator(this, true);
    }

    long goToObject(long objectID);

    // Read a STAR loop structure
    long int readStarLoop(std::ifstream &in, bool do_only_count = false);

    /* Read a STAR list
     * The function returns true if the list is followed by a loop, false otherwise */
    bool readStarList(std::ifstream &in);

    /* Read a MetaDataTable from a STAR-format data block
     *
     * If the data block contains a list and a table, the function will return 2,
     * the first time it is called and the list is read into the MetaDataTable
     * in that case the function needs to be called another time. The second time
     * it will read the _loop structure into the MetaDataTable and 1 will be returned
     *
     * If the data block contains only a list or a table, it is read in the MetaDataTable and the function will return 1
     *
     * If no data block is found the function will return 0 and the MetaDataTable remains empty
     */
    long int readStar(std::ifstream &in, const std::string &name = "", bool do_only_count = false);

    // Read a MetaDataTable (get file format from extension)
    long int read(const FileName &filename, const std::string &name = "", bool do_only_count = false);

    // Write a MetaDataTable in STAR format
    void write(std::ostream &out = std::cout);

    // Write to a single file
    void write(const FileName &fn_out);

    // Make a histogram of a column
    void columnHistogram(
        EMDL::EMDLabel label, std::vector<RFLOAT> &histX, std::vector<RFLOAT> &histY,
        int verb = 0, CPlot2D *plot2D = NULL, long int nr_bin = -1,
        RFLOAT hist_min = -LARGE_NUMBER, RFLOAT hist_max = LARGE_NUMBER,
        bool do_fractional_instead = false, bool do_cumulative_instead = false
    );

    static void histogram(
        std::vector<RFLOAT> &values, std::vector<RFLOAT> &histX, std::vector<RFLOAT> &histY,
        int verb = 0, std::string title="Histogram", CPlot2D *plot2D = NULL, long int nr_bin = -1,
        RFLOAT hist_min = -LARGE_NUMBER, RFLOAT hist_max = LARGE_NUMBER,
        bool do_fractional_instead = false, bool do_cumulative_instead = false
    );

    void addToCPlot2D(
        CPlot2D *plot2D, EMDL::EMDLabel xaxis, EMDL::EMDLabel yaxis,
        double red=0., double green=0., double blue=0., double linewidth = 1.0,
        std::string marker=""
    );

    void printLabels(std::ostream &ost);

    // Randomise the order inside the STAR file
    void randomiseOrder();

    // 14 Feb 2017 - Shaoda, Check whether the two MetaDataTables contain the same set of activeLabels
    static bool compareLabels(const MetaDataTable &MD1, const MetaDataTable &MD2);

    // Join 2 metadata tables. Only include labels that are present in both of them.
    static MetaDataTable combineMetaDataTables(std::vector<MetaDataTable> &MDin);

    template<class T>
    bool isTypeCompatible(EMDL::EMDLabel label, T &value) const;

    private:

    // Check if 'id' corresponds to an actual object.
    // Throw if it does not.
    void checkObjectID(long id) const throw (std::string);

    /* setObjectUnsafe(data)
     *  Same as setObject, but assumes that all labels are present. */
    void setObjectUnsafe(MetaDataContainer *data, long objId);

};

void compareMetaDataTable(
    MetaDataTable &MD1, MetaDataTable &MD2,
    MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
    EMDL::EMDLabel label1, double eps = 0.0,
    EMDL::EMDLabel label2 = EMDL::UNDEFINED,
    EMDL::EMDLabel label3 = EMDL::UNDEFINED
);

// find a subset of the input metadata table that has corresponding entries between the specified min and max values
MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label,
    RFLOAT min_value, RFLOAT max_value
);

// find a subset of the input metadata table that has corresponding entries with or without a given substring
MetaDataTable subsetMetaDataTable(
    MetaDataTable &MDin, EMDL::EMDLabel label,
    std::string search_str, bool exclude=false
);

// remove duplicated particles that are in the same micrograph (mic_label) and within a given threshold [px]
// OriginX/Y are multiplied by origin_scale before added to CoordinateX/Y to compensate for down-sampling
MetaDataTable removeDuplicatedParticles(
    MetaDataTable &MDin, EMDL::EMDLabel mic_label,
    RFLOAT threshold, RFLOAT origin_scale=1.0,
    FileName fn_removed="", bool verb=true
);

#ifdef METADATA_TABLE_TYPE_CHECK
//#pragma message("typecheck enabled")
template<class T>
bool MetaDataTable::isTypeCompatible(EMDL::EMDLabel label, T& value) const {
    // remove const appended by setValue()
    typedef typename std::remove_const<T>::type U;

    // In C++11, this repeat can be avoided by using "if constexpr(...) else static_assert"
    static_assert(
        std::is_same<bool, U>::value ||
        std::is_same<FileName, U>::value || std::is_same<std::string, U>::value ||
        std::is_same<double, U>::value || std::is_same<float, U>::value ||
        std::is_same<int, U>::value || std::is_same<long, U>::value ||
        std::is_same<std::vector<double>, U>::value || std::is_same<std::vector<float>, U>::value,
        "Compile error: wrong type given to MetaDataTable::getValur or setValue"
    );

    if (std::is_same<bool, U>::value) {
        return EMDL::isBool(label);
    } else if (std::is_same<FileName, U>::value || std::is_same<std::string, U>::value) {
        return EMDL::isString(label);
    } else if (std::is_same<double, U>::value || std::is_same<float, U>::value) {
        return EMDL::isDouble(label);
    } else if (std::is_same<int, U>::value || std::is_same<long, U>::value) {
        return EMDL::isInt(label);
    } else if (std::is_same<std::vector<double>, U>::value || std::is_same<std::vector<float>, U>::value) {
        return EMDL::isVector(label);
    } else {
        return false;
    }
}
#endif

/// Return value or, if the label does not exist, raise error.
// objectID is 0-indexed.
template<typename T>
T MetaDataTable::getValue(EMDL::EMDLabel label, long objectID) const {

    T value;
    if (label < 0 || label >= EMDL::LAST_LABEL) throw "Label not recognised";

    if (label == EMDL::UNKNOWN_LABEL) REPORT_ERROR("MetaDataTable::setValue does not support unknown label.");

    #ifdef METADATA_TABLE_TYPE_CHECK
    if (!isTypeCompatible(label, value))
        REPORT_ERROR("Runtime error: wrong type given to MetaDataTable::getValue for label " + EMDL::label2Str(label));
    #endif

    const long off = label2offset[label];
    if (off > -1) {
        if (objectID < 0) {
            objectID = current_objectID;
        } else {
            try { checkObjectID(objectID); }
            catch (const std::string &errmsg) {
                REPORT_ERROR((std::string) __func__ + ": " + errmsg);
            }
        }

        value = objects[objectID]->getValue<T>(off);
    } else {
        throw "Negative offset";
    }
    return value;
}


/// TODO: value rather than &value
template<class T>
bool MetaDataTable::setValue(EMDL::EMDLabel label, const T &value, long int objectID) {

    if (label < 0 || label >= EMDL::LAST_LABEL) return false;

    if (label == EMDL::UNKNOWN_LABEL) REPORT_ERROR("MetaDataTable::setValue does not support unknown label.");

    #ifdef METADATA_TABLE_TYPE_CHECK
    if (!isTypeCompatible(label, value)) REPORT_ERROR("Runtime error: wrong type given to MetaDataTable::setValue for label " + EMDL::label2Str(label));
    #endif

    long off = label2offset[label];

    if (off < 0) {
        addLabel(label);
        off = label2offset[label];
    }

    if (objectID < 0) {
        objectID = current_objectID;
    } else {
        try { checkObjectID(objectID); }
        catch (const std::string &errmsg) {
            REPORT_ERROR((std::string) __func__ + ": " + errmsg);
        }
    }

    if (off <= -1) return false;

    objects[objectID]->setValue(off, value);
    return true;
}

#endif
