/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef METADATA_CONTAINER_H
#define METADATA_CONTAINER_H

#include "src/funcs.h"
#include "src/metadata_label.h"

class MetaDataTable;

struct MetaDataContainer {

    MetaDataTable* table;

    std::vector<double> doubles;  // Extended precision
    std::vector<long> ints;       // Extended precision
    std::vector<bool> bools;
    std::vector<std::string> strings;
    std::vector<std::vector<double>> doubleVectors;
    std::vector<std::string> unknowns;

    MetaDataContainer();
    MetaDataContainer(
        MetaDataTable* table, long doubleCount, long intCount, long boolCount,
        long stringCount, long doubleVectorCount, long unknownCount
    );
    MetaDataContainer(MetaDataTable* table, MetaDataContainer* mdc);

    template <typename T>
    T getValue(long i) const;

    void setValue(long i, double src);
    void setValue(long i, float src);
    void setValue(long i, int src);
    void setValue(long i, long src);
    void setValue(long i, bool src);
    void setValue(long i, const std::string &src);
    void setValue(long i, const std::vector<double> &src);
    void setValue(long i, const std::vector<float> &src);
};

#endif
