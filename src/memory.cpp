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
#include "src/memory.h"

char* askMemory(size_t memsize) throw (RelionError) {
    char *ptr = (char *) calloc(memsize, sizeof(char));
    if (!ptr) {
        std::cerr << "Failed to allocate " <<  memsize << " bytes." << std::endl;
        REPORT_ERROR("Error in askMemory");
    }
    return ptr;
}

int freeMemory(void *ptr, size_t memsize) {

    if (memsize < 1) return -1;

    free(ptr);
    ptr = nullptr;
    return 0;
}
