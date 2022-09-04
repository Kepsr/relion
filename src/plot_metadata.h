#pragma once
#include "src/metadata_table.h"
#include "src/CPlot2D.h"

struct PlotMetaData {

    // Make a histogram of a column
    static void columnHistogram(
        const MetaDataTable &mdt,
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

    static void addToCPlot2D(
        const MetaDataTable &mdt,
        CPlot2D *plot2D, EMDL::EMDLabel xaxis, EMDL::EMDLabel yaxis,
        double red=0.0, double green=0.0, double blue=0.0, double linewidth = 1.0,
        std::string marker=""
    );

};
