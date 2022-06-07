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

MetaDataContainer::MetaDataContainer()
    :   doubles(0), ints(0), bools(0), strings(0), doubleVectors(0), unknowns(0) {}


MetaDataContainer::MetaDataContainer(
        MetaDataTable *table, long doubleCount, long intCount,
        long boolCount, long stringCount, long doubleVectorCount, long unknownCount)
: table(table),
  doubles(doubleCount, 0),
  ints(intCount, 0),
  bools(boolCount, false),
  strings(stringCount, ""),
  doubleVectors(doubleVectorCount),
  unknowns(unknownCount) {}

MetaDataContainer::MetaDataContainer(
        MetaDataTable *table, MetaDataContainer* mdc)
: table(table),
  doubles(mdc->doubles),
  ints(mdc->ints),
  bools(mdc->bools),
  strings(mdc->strings),
  doubleVectors(mdc->doubleVectors),
  unknowns(mdc->unknowns) {}

template <>
double MetaDataContainer::getValue(long offset) const {
    return doubles[offset];
}

template <>
float MetaDataContainer::getValue(long offset) const {
    return (float) doubles[offset];
}

template <>
int MetaDataContainer::getValue(long offset) const {
    return (int) ints[offset];
}

template <>
long MetaDataContainer::getValue(long offset) const {
    return ints[offset];
}

template <>
bool MetaDataContainer::getValue(long offset) const {
    return bools[offset];
}

template <>
std::vector<double> MetaDataContainer::getValue(long offset) const {
	return doubleVectors[offset];
}

template <>
std::vector<float> MetaDataContainer::getValue(long offset) const {
    std::vector<float> result (doubleVectors[offset].size());
	std::copy(doubleVectors[offset].begin(), doubleVectors[offset].end(), result.begin());
    return result;
}

template <>
std::string MetaDataContainer::getValue(long offset) const {
	return strings[offset] == "\"\"" ? "" : strings[offset];
}

void MetaDataContainer::setValue(long offset, double src) {
    doubles[offset] = src;
}

void MetaDataContainer::setValue(long offset, float src) {
    doubles[offset] = src;
}

void MetaDataContainer::setValue(long offset, int src) {
    ints[offset] = src;
}

void MetaDataContainer::setValue(long offset, long src) {
    ints[offset] = src;
}

void MetaDataContainer::setValue(long offset, bool src) {
    bools[offset] = src;
}

void MetaDataContainer::setValue(long offset, const std::string &src) {
	strings[offset] = src.empty() ? "\"\"" : src;
}

void MetaDataContainer::setValue(long offset, const std::vector<double> &src) {
	doubleVectors[offset] = src;
}

void MetaDataContainer::setValue(long offset, const std::vector<float> &src) {
	doubleVectors[offset].resize(src.size());
	std::copy(src.begin(), src.end(), doubleVectors[offset].begin());
}

