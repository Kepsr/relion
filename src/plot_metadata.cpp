#include "src/plot_metadata.h"

void PlotMetaData::columnHistogram(
    const MetaDataTable &mdt,
    EMDL::EMDLabel label, std::vector<RFLOAT> &histX, std::vector<RFLOAT> &histY,
    int verb, CPlot2D *plot2D,
    long int nr_bin, RFLOAT hist_min, RFLOAT hist_max,
    bool do_fractional_instead, bool do_cumulative_instead
) {
    if (!mdt.containsLabel(label))
        REPORT_ERROR("ERROR: The column specified is not present in the MetaDataTable.");

    std::vector<RFLOAT> values;
    values.reserve(mdt.numberOfObjects());
    for (auto i : mdt) {
        if (EMDL::is<double>(label)) {
            values.push_back(mdt.getValue<double>(label, i));
        } else if (EMDL::is<int>(label)) {
            values.push_back(mdt.getValue<long>(label, i));
        } else if (EMDL::is<bool>(label)) {
            values.push_back(mdt.getValue<bool>(label, i));
        } else {
            REPORT_ERROR("Cannot use --stat_column for this type of column");
        }
    }

    const std::string title = EMDL::label2Str(label);
    histogram(values, histX, histY, verb, title, plot2D, nr_bin, hist_min, hist_max, do_fractional_instead, do_cumulative_instead);
}

void PlotMetaData::histogram(
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

void PlotMetaData::addToCPlot2D(
    const MetaDataTable &mdt,
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
    for (long int idx = 0; idx < mdt.objects.size(); idx++) {
        const long offx = mdt.label_indices[xaxis];
        if (offx < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find x-axis label");

        if (xaxis == EMDL::UNDEFINED) {
            xval = idx + 1;
        } else if (EMDL::is<double>(xaxis)) {
            xval = mdt.objects[idx]->getValue<double>(offx);
        } else if (EMDL::is<int>(xaxis)) {
            xval = mdt.objects[idx]->getValue<int>(offx);
        } else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot x-axis double, int or long int");

        const long offy = mdt.label_indices[yaxis];
        if (offy < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find y-axis label");

        if (EMDL::is<double>(yaxis)) {
            yval = mdt.objects[idx]->getValue<double>(offy);
        } else if (EMDL::is<int>(yaxis)) {
            yval = mdt.objects[idx]->getValue<int>(offy);
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
