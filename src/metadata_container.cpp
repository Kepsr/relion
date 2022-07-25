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
#include "src/metadata_container.h"

MetaDataContainer::MetaDataContainer():
doubles(0), ints(0), bools(0), strings(0), doubleVectors(0), unknowns(0) {}


MetaDataContainer::MetaDataContainer(
    MetaDataTable *table, long doubleCount, long intCount,
    long boolCount, long stringCount, long doubleVectorCount, long unknownCount
):
table(table),
doubles(doubleCount, 0),
ints(intCount, 0),
bools(boolCount, false),
strings(stringCount, ""),
doubleVectors(doubleVectorCount),
unknowns(unknownCount) {}

MetaDataContainer::MetaDataContainer(MetaDataTable *table, MetaDataContainer* mdc):
table(table),
doubles(mdc->doubles),
ints(mdc->ints),
bools(mdc->bools),
strings(mdc->strings),
doubleVectors(mdc->doubleVectors),
unknowns(mdc->unknowns) {}

template <>
double MetaDataContainer::getValue(long i) const {
    return doubles[i];
}

template <>
float MetaDataContainer::getValue(long i) const {
    return (float) doubles[i];
}

template <>
int MetaDataContainer::getValue(long i) const {
    return (int) ints[i];
}

template <>
long MetaDataContainer::getValue(long i) const {
    return ints[i];
}

template <>
bool MetaDataContainer::getValue(long i) const {
    return bools[i];
}

template <>
std::vector<double> MetaDataContainer::getValue(long i) const {
    return doubleVectors[i];
}

template <>
std::vector<float> MetaDataContainer::getValue(long i) const {
    const auto &v = doubleVectors[i];
    std::vector<float> result (v.size());
    std::copy(v.begin(), v.end(), result.begin());
    return result;
}

template <>
std::string MetaDataContainer::getValue(long i) const {
    const auto &s = strings[i];
    return s == "\"\"" ? "" : s;
}

void MetaDataContainer::setValue(long i, double src) {
    doubles[i] = src;
}

void MetaDataContainer::setValue(long i, float src) {
    doubles[i] = src;
}

void MetaDataContainer::setValue(long i, int src) {
    ints[i] = src;
}

void MetaDataContainer::setValue(long i, long src) {
    ints[i] = src;
}

void MetaDataContainer::setValue(long i, bool src) {
    bools[i] = src;
}

void MetaDataContainer::setValue(long i, const std::string &src) {
    strings[i] = src.empty() ? "\"\"" : src;
}

void MetaDataContainer::setValue(long i, const std::vector<double> &src) {
    doubleVectors[i] = src;
}

void MetaDataContainer::setValue(long i, const std::vector<float> &src) {
    auto &v = doubleVectors[i];
    v.resize(src.size());
    std::copy(src.begin(), src.end(), v.begin());
}

