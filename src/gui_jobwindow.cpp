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
#include "src/gui_jobwindow.h"
#include "src/gui_mainwindow.h"  // For GroupContext
#include "src/pipeliner.h"
#include <unordered_map>

JobWindow::JobWindow(int x, int y, int w, int h, const char* title):
Fl_Box(this->x, this->y, this->w, this->h, title) {
    clear();
    this->x = x; this->y = y; this->w = w; this->h = h;
}

void JobWindow::clear() {
    tabs = nullptr;
    for (int i = 0; i < 7; ++i) { tabcarr[i] = nullptr; }
    runtab = nullptr;
    group1 = group2 = group3 = group4 = group5 = group6 = group7 = queue_group = nullptr;
    current_y = start_y = 0;
    is_continue = false;
    guientries.clear();
}

void JobWindow::setupTabs(int nr_tabs) {
    current_y = y;  // top of the GUI

    char *my_allow_change_dedicated = getenv("RELION_ALLOW_CHANGE_MINIMUM_DEDICATED");
    if (!my_allow_change_dedicated) {
        do_allow_change_minimum_dedicated = DEFAULT::MINIMUMDEDICATED;
    } else {
        int check_allow = textToInteger(my_allow_change_dedicated);
        do_allow_change_minimum_dedicated = check_allow != 0;
    }

    // Set up tabs
    if (nr_tabs >= 1) {
        // there is always the running tab, which is not counted on the input nr_tabs!
        tabs = new Fl_Tabs(x, current_y, w, h - MENUHEIGHT);
        current_y += TABHEIGHT;
        tabs->begin();
        for (int i = 0; i < nr_tabs && i < 7; ++i) {
            GroupContext (tabcarr[i] = new Fl_Group(x, current_y, w, h - MENUHEIGHT, ""));
            tabcarr[i]->color(GUI_BACKGROUND_COLOR);
            tabcarr[i]->selection_color(GUI_BACKGROUND_COLOR2);
        }
        if (nr_tabs >= 8) {
            std::cerr << "ERROR: only 7 job-specific tabs implemented..." << std::endl;
            exit(1);
        }
        current_y += 15;
        start_y = current_y;

        {
        GroupContext context (runtab = new Fl_Group(x, current_y, w, h - MENUHEIGHT, ""));
        runtab->label("Running");
        // Fill this in later, when we have the joboptions
        }
        setupRunTab();
        runtab->color(GUI_BACKGROUND_COLOR);
        runtab->selection_color(GUI_BACKGROUND_COLOR2);

        tabs->end();
    }
}

void JobWindow::setupRunTab() {
    GroupContext context (runtab);

    resetHeight();

    bool has_parallel = false;

    if (myjob.joboptions.find("nr_mpi") != myjob.joboptions.end()) {
        place("nr_mpi", TOGGLE_LEAVE_ACTIVE);
        has_parallel = true;
    }

    if (myjob.joboptions.find("nr_threads") != myjob.joboptions.end()) {
        place("nr_threads", TOGGLE_LEAVE_ACTIVE);
        has_parallel = true;
    }

    // Add a little spacer
    if (has_parallel) { current_y += STEPY / 4; }

    // Set up queue groups for running tab
    GroupContext (queue_group = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_queue", TOGGLE_LEAVE_ACTIVE, queue_group);
    {
    GroupContext context (queue_group);

    place("queuename");

    place("qsub");

    char *extra_count_text = getenv("RELION_QSUB_EXTRA_COUNT");
    const char extra_count_val = extra_count_text ? atoi(extra_count_text) : 2;
    for (int i = 1; i <= extra_count_val; i++) {
        const std::string fn = "qsub_extra" + std::to_string(i);
        if (myjob.joboptions.find(fn) != myjob.joboptions.end())
            place(fn);
    }

    place("qsubscript");

    place("min_dedicated");
    guientries["min_dedicated"].deactivate(!do_allow_change_minimum_dedicated);

    }
    guientries["do_queue"].cb_menu_i();  // This is to make the default effective

    // Add a little spacer
    current_y += STEPY / 4;

    place("other_args");

}

void JobWindow::place(std::string key, int deactivate_option, Fl_Group *deactivate_this_group, bool actually_activate) {
    if (myjob.joboptions.find(key) == myjob.joboptions.end())
        std::cerr << "WARNING: cannot find " << key << " in the defined joboptions of jobtype= " << myjob.type << std::endl;

    guientries[key].place(myjob.joboptions[key], current_y, deactivate_option, deactivate_this_group, actually_activate);
}

void JobWindow::place2(std::string key1, std::string key2, std::string label, int deactivate_option) {

    if (myjob.joboptions.find(key1) == myjob.joboptions.end())
    std::cerr << "WARNING: cannot find " << key1 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
    if (myjob.joboptions.find(key2) == myjob.joboptions.end())
    std::cerr << "WARNING: cannot find " << key2 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;

    myjob.joboptions[key1].label_gui = label;
    myjob.joboptions[key2].label_gui = "";
    int old_y = current_y;
    guientries[key1].place(
        myjob.joboptions[key1], current_y, deactivate_option, nullptr, false,
        XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2
    );
    current_y = old_y;
    guientries[key2].place(
        myjob.joboptions[key2], current_y, deactivate_option, nullptr, false,
        XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2
    );
}

void JobWindow::place3(
    std::string key1, std::string key2, std::string key3, std::string label,
    int deactivate_option
) {

    if (myjob.joboptions.find(key1) == myjob.joboptions.end())
    std::cerr << "WARNING: cannot find " << key1 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
    if (myjob.joboptions.find(key2) == myjob.joboptions.end())
    std::cerr << "WARNING: cannot find " << key2 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
    if (myjob.joboptions.find(key3) == myjob.joboptions.end())
    std::cerr << "WARNING: cannot find " << key3 << " in the defined joboptions of jobtype= "  << myjob.type<< std::endl;

    myjob.joboptions[key1].label_gui = label;
    myjob.joboptions[key2].label_gui = "";
    myjob.joboptions[key3].label_gui = "";
    int old_y = current_y;
    guientries[key1].place(
        myjob.joboptions[key1], current_y, deactivate_option, nullptr, false,
        XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3
    );
    current_y = old_y;
    guientries[key2].place(
        myjob.joboptions[key2], current_y, deactivate_option, nullptr, false,
        XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3
    );
    current_y = old_y;
    guientries[key3].place(
        myjob.joboptions[key3], current_y, deactivate_option, nullptr, false,
        XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3
    );
}

void JobWindow::toggle_new_continue(bool is_continue) {
    myjob.is_continue = this->is_continue = is_continue;

    for (auto &guientry : guientries) {
        switch (guientry.second.deactivate_option) {

            case TOGGLE_DEACTIVATE:
            guientry.second.deactivate(is_continue);
            break;

            case TOGGLE_REACTIVATE:
            guientry.second.deactivate(!is_continue);
            break;

            case TOGGLE_ALWAYS_DEACTIVATE:
            guientry.second.deactivate(true);
            break;

            case TOGGLE_LEAVE_ACTIVE:
            // do nothing
            break;

            default:
            REPORT_ERROR("ERROR: unrecognised deactivate-option for GUI entry " + guientry.first);

        }
    }
}

void JobWindow::resetHeight() { current_y = start_y; }

// Update all values in the Fl_Input entries from the corresponding job_options
void JobWindow::updateMyGui() {
    for (auto &guientry : guientries) {
        if (myjob.joboptions.find(guientry.first) == myjob.joboptions.end())
            std::cerr << "WARNING: cannot find " << guientry.first << " in the defined joboptions!" <<std::endl;

        guientry.second.setValue(myjob.joboptions[guientry.first].value);
    }
}

// Update all values in the Fl_Input entries into the corresponding job_options
void JobWindow::updateMyJob() {
    for (auto &option : myjob.joboptions) {
        if (guientries.find(option.first) == guientries.end()) {
            std::cerr << "ERROR: cannot find " << option.first << " in the defined joboptions!" << std::endl;
            REPORT_ERROR("Stopping now...");
        }

        option.second.value = std::string(guientries[option.first].inp->value());
    }
}

void JobWindow::initialise(int job_type) {

    // A mapping from ints (job types) to member function pointers
    const std::unordered_map<int, void (JobWindow::*)()> i2mfp {
        {Process::IMPORT,       &JobWindow::initialiseImportWindow},
        {Process::MOTIONCORR,   &JobWindow::initialiseMotioncorrWindow},
        {Process::CTFFIND,      &JobWindow::initialiseCtffindWindow},
        {Process::MANUALPICK,   &JobWindow::initialiseManualpickWindow},
        {Process::AUTOPICK,     &JobWindow::initialiseAutopickWindow},
        {Process::EXTRACT,      &JobWindow::initialiseExtractWindow},
        {Process::CLASSSELECT,  &JobWindow::initialiseSelectWindow},
        {Process::CLASS2D,      &JobWindow::initialiseClass2DWindow},
        {Process::INIMODEL,     &JobWindow::initialiseInimodelWindow},
        {Process::CLASS3D,      &JobWindow::initialiseClass3DWindow},
        {Process::AUTO3D,       &JobWindow::initialiseAutorefineWindow},
        {Process::MULTIBODY,    &JobWindow::initialiseMultiBodyWindow},
        {Process::MASKCREATE,   &JobWindow::initialiseMaskcreateWindow},
        {Process::JOINSTAR,     &JobWindow::initialiseJoinstarWindow},
        {Process::SUBTRACT,     &JobWindow::initialiseSubtractWindow},
        {Process::POST,         &JobWindow::initialisePostprocessWindow},
        {Process::RESMAP,       &JobWindow::initialiseLocresWindow},
        {Process::MOTIONREFINE, &JobWindow::initialiseMotionrefineWindow},
        {Process::CTFREFINE,    &JobWindow::initialiseCtfrefineWindow},
        {Process::EXTERNAL,     &JobWindow::initialiseExternalWindow}
    };

    auto it = i2mfp.find(job_type);
    if (it == i2mfp.end())
        REPORT_ERROR("ERROR: unrecognised job-type to add to the GUI");

    myjob.initialise(job_type);
    (this->*(it->second))();  // Call the member function through its pointer

    // read settings if hidden file exists
    myjob.read("", is_continue);

    // update the window
    updateMyGui();
}

void JobWindow::initialiseImportWindow() {

    setupTabs(2);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("Movies/mics");
    resetHeight();

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_raw", TOGGLE_DEACTIVATE, group1, false);
    current_y += STEPY / 2;  // Add a little spacer
    {
    GroupContext context (group1);
    place("fn_in_raw");
    place("is_multiframe");

    current_y += STEPY / 2;  // Add a little spacer

    place("optics_group_name");
    place("fn_mtf");
    place("angpix");
    place("kV");
    place("Cs");
    place("Q0");
    place("beamtilt_x");
    place("beamtilt_y");
    }
    guientries["do_raw"].cb_menu_i();  // make default active

    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Others");
    resetHeight();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_other", TOGGLE_DEACTIVATE, group2, false);
    {
    GroupContext context (group2);
    current_y += STEPY / 2;  // Add a little spacer
    place("fn_in_other");
    place("node_type");
    current_y += STEPY / 2;  // Add a little spacer
    place("optics_group_particles");
    }
    guientries["do_other"].cb_menu_i();  // make default active

    }

    // Always deactivate the queue option
    guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
    myjob.joboptions["do_queue"].setString("No");
}

void JobWindow::initialiseMotioncorrWindow() {

    setupTabs(2);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("input_star_mics", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("first_frame_sum", TOGGLE_DEACTIVATE);
    place("last_frame_sum", TOGGLE_DEACTIVATE);
    place("dose_per_frame", TOGGLE_DEACTIVATE);
    place("pre_exposure", TOGGLE_DEACTIVATE);
    place("eer_grouping", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_dose_weighting", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("do_save_noDW", TOGGLE_DEACTIVATE);
    }
    guientries["do_dose_weighting"].cb_menu_i();  // make default active

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_save_ps", TOGGLE_DEACTIVATE, group2);

    {
    GroupContext context (group2);
    place("group_for_ps", TOGGLE_DEACTIVATE);
    }

    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Motion");
    resetHeight();

    place("bfactor", TOGGLE_DEACTIVATE);
    place2("patch_x", "patch_y", "Number of patches X, Y", TOGGLE_DEACTIVATE);
    place("group_frames", TOGGLE_DEACTIVATE);
    place("bin_factor", TOGGLE_DEACTIVATE);
    place("fn_gain_ref", TOGGLE_DEACTIVATE);
    place("gain_rot", TOGGLE_DEACTIVATE);
    place("gain_flip", TOGGLE_DEACTIVATE);
    place("fn_defect", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;
    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_own_motioncor", TOGGLE_DEACTIVATE, group4, true);
    {
    GroupContext context (group4);
    place("fn_motioncor2_exe", TOGGLE_DEACTIVATE);
    place("gpu_ids");
    place("other_motioncor2_args", TOGGLE_DEACTIVATE);
    }
    guientries["do_own_motioncor"].cb_menu_i();  // make default active

    }
}

void JobWindow::initialiseCtffindWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("input_star_mics", TOGGLE_DEACTIVATE);
    place("use_noDW", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_phaseshift", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place3("phase_min", "phase_max", "phase_step", "Phase shift - Min, Max, Step (deg)", TOGGLE_DEACTIVATE);
    }
    guientries["do_phaseshift"].cb_menu_i();  // make default active

    current_y += STEPY / 2;  // Add a little spacer

    place("dast", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("CTFFIND-4.1");
    resetHeight();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_ctffind4", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place("fn_ctffind_exe", TOGGLE_DEACTIVATE);
    place("use_given_ps", TOGGLE_DEACTIVATE);
    place("slow_search", TOGGLE_DEACTIVATE);
    place("ctf_win", TOGGLE_DEACTIVATE);
    }
    guientries["use_ctffind4"].cb_menu_i();  // make default active

    current_y += STEPY / 2;  // Add a little spacer

    place("box", TOGGLE_DEACTIVATE);
    place("resmin", TOGGLE_DEACTIVATE);
    place("resmax", TOGGLE_DEACTIVATE);
    place("dfmin", TOGGLE_DEACTIVATE);
    place("dfmax", TOGGLE_DEACTIVATE);
    place("dfstep", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Gctf");
    resetHeight();

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_gctf", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);
    place("fn_gctf_exe", TOGGLE_DEACTIVATE);
    place("do_ignore_ctffind_params", TOGGLE_DEACTIVATE);
    place("do_EPA", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer
    place("other_gctf_args", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer
    place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["use_gctf"].cb_menu_i();  // make default active

    }
}

void JobWindow::initialiseManualpickWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_in", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Display");
    resetHeight();

    place("diameter");
    place("micscale");
    place("sigma_contrast");
    place("white_val");
    place("black_val");

    current_y += STEPY / 2;
    place("lowpass");
    place("highpass");
    place("angpix");

    current_y += STEPY / 2;
    place ("do_startend");

    current_y += STEPY / 2;
    place("ctfscale");
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Colors");

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("do_color", TOGGLE_LEAVE_ACTIVE, group1);
    {
    GroupContext context (group1);
    place("color_label");
    place("fn_color");
    place("blue_value");
    place("red_value");
    }
    guientries["do_color"].cb_menu_i();  // make default active

    }
    // Always deactivate the queue option
    guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
    myjob.joboptions["do_queue"].setString("No");
}

void JobWindow::initialiseAutopickWindow() {

    setupTabs(5);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_input_autopick", TOGGLE_DEACTIVATE);
    place("angpix", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;

    place("fn_refs_autopick", TOGGLE_DEACTIVATE);

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_ref3d", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("fn_ref3d_autopick", TOGGLE_DEACTIVATE);
    place("ref3d_symmetry", TOGGLE_DEACTIVATE);
    place("ref3d_sampling", TOGGLE_DEACTIVATE);
    }
    guientries["do_ref3d"].cb_menu_i();

    place("do_log", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Laplacian");
    resetHeight();

    place("log_diam_min", TOGGLE_DEACTIVATE);
    place("log_diam_max", TOGGLE_DEACTIVATE);
    place("log_invert", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer
    place("log_maxres", TOGGLE_DEACTIVATE);
    place("log_adjust_thr");
    place("log_upper_thr");
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("References");
    resetHeight();

    //set up group
    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("lowpass", TOGGLE_DEACTIVATE);
    place("highpass", TOGGLE_DEACTIVATE);
    place("angpix_ref", TOGGLE_DEACTIVATE);
    place("particle_diameter", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("psi_sampling_autopick", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("do_invert_refs", TOGGLE_DEACTIVATE);

    place("do_ctf_autopick", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place("do_ignore_first_ctfpeak_autopick", TOGGLE_DEACTIVATE);  //(current_y, "Ignore CTFs until first peak?", false, "Set this to Yes, only if this option was also used to generate the references.");
    }
    guientries["do_ctf_autopick"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Autopicking");
    resetHeight();

    place("threshold_autopick");
    place("mindist_autopick");
    place("maxstddevnoise_autopick");
    place("minavgnoise_autopick");

    current_y += STEPY / 2;

    place("do_write_fom_maps");
    place("do_read_fom_maps");

    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    place("shrink", TOGGLE_DEACTIVATE);

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group3);
    {
    GroupContext context (group3);
    place("gpu_ids");
    }
    guientries["use_gpu"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[4]);
    tabcarr[4]->label("Helix");
    resetHeight();

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_pick_helical_segments", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);

    place("do_amyloid");
    place("helical_tube_outer_diameter");

    current_y += STEPY / 2;

    place("helical_nr_asu");
    place("helical_rise");

    current_y += STEPY / 2;

    place("helical_tube_kappa_max");
    place("helical_tube_length_min");
    }
    guientries["do_pick_helical_segments"].cb_menu_i();

    }
}

void JobWindow::initialiseExtractWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("star_mics", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;
    place("coords_suffix", TOGGLE_DEACTIVATE);
    current_y += STEPY / 2;

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_reextract", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("fndata_reextract", TOGGLE_DEACTIVATE);
    place("do_reset_offsets", TOGGLE_DEACTIVATE);

    GroupContext (group7 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_recenter", TOGGLE_DEACTIVATE, group7);
    {
    GroupContext context (group7);
    place3("recenter_x", "recenter_y", "recenter_z", "Recenter on - X, Y, Z (pix):", TOGGLE_DEACTIVATE);
    }
    guientries["do_recenter"].cb_menu_i();

    }
    guientries["do_reextract"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Extract");
    resetHeight();

    place("extract_size", TOGGLE_DEACTIVATE);  //(current_y, "Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
    place("do_invert", TOGGLE_DEACTIVATE);  //(current_y, "Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted.");

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_norm", TOGGLE_DEACTIVATE, group3);
    {
    GroupContext context (group3);

    place("bg_diameter", TOGGLE_DEACTIVATE);  //(current_y, "Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");

    place("white_dust", TOGGLE_DEACTIVATE);  //(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

    place("black_dust", TOGGLE_DEACTIVATE);  //(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
    }
    guientries["do_norm"].cb_menu_i();

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_rescale", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);
    place("rescale", TOGGLE_DEACTIVATE);
    }
    guientries["do_rescale"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Helix");
    resetHeight();

    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_extract_helix", TOGGLE_DEACTIVATE, group5);
    {
    GroupContext context (group5);

    place("helical_tube_outer_diameter", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;

    place("helical_bimodal_angular_priors", TOGGLE_DEACTIVATE);

    GroupContext (group6 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));
    current_y += STEPY / 2;

    place("do_extract_helical_tubes", TOGGLE_DEACTIVATE, group6);
    {
    GroupContext context (group6);
    place("do_cut_into_segments", TOGGLE_DEACTIVATE);
    place("helical_nr_asu", TOGGLE_DEACTIVATE);
    place("helical_rise", TOGGLE_DEACTIVATE);
    }
    guientries["do_extract_helical_tubes"].cb_menu_i();

    }
    guientries["do_extract_helix"].cb_menu_i();

    }
}

void JobWindow::initialiseSelectWindow() {

    setupTabs(4);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_model", TOGGLE_DEACTIVATE);
    place("fn_mic", TOGGLE_DEACTIVATE);
    place("fn_data", TOGGLE_DEACTIVATE);
    place("fn_coords", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Class options");
    resetHeight();

    place("do_recenter", TOGGLE_DEACTIVATE);
    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_regroup", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("nr_groups", TOGGLE_DEACTIVATE);
    }
    guientries["do_regroup"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Subsets");
    resetHeight();

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_select_values", TOGGLE_DEACTIVATE, group3);
    {
    GroupContext context (group3);
    place("select_label", TOGGLE_DEACTIVATE);
    place("select_minval", TOGGLE_DEACTIVATE);
    place("select_maxval", TOGGLE_DEACTIVATE);
    }
    guientries["do_select_values"].cb_menu_i();

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    current_y += STEPY / 2;  // Add a little spacer

    place("do_discard", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);
    place("discard_label", TOGGLE_DEACTIVATE);
    place("discard_sigma", TOGGLE_DEACTIVATE);
    }
    guientries["do_discard"].cb_menu_i();

    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    current_y += STEPY / 2;  // Add a little spacer

    place("do_split", TOGGLE_DEACTIVATE, group5);
    {
    GroupContext context (group5);
    place("do_random", TOGGLE_DEACTIVATE);
    place("split_size", TOGGLE_DEACTIVATE);
    place("nr_split", TOGGLE_DEACTIVATE);
    }
    guientries["do_split"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Duplicates");
    resetHeight();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_remove_duplicates", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place("duplicate_threshold", TOGGLE_DEACTIVATE);
    place("image_angpix", TOGGLE_DEACTIVATE);
    }
    guientries["do_remove_duplicates"].cb_menu_i();

    }
    // Always deactivate the queue option
    guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
    myjob.joboptions["do_queue"].setString("No");
}

void JobWindow::initialiseClass2DWindow() {

    setupTabs(6);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_img", TOGGLE_DEACTIVATE);
    place("fn_cont", TOGGLE_REACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("CTF");

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
    }
    guientries["do_ctf_correction"].cb_menu_i();  // To make default effective

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Optimisation");
    resetHeight();

    //set up groups
    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("nr_classes", TOGGLE_DEACTIVATE);
    place("tau_fudge");

    current_y += STEPY / 2;  // Add a little spacer

    place("nr_iter");
    place("do_fast_subsets", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("particle_diameter");
    place("do_zero_mask", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("highres_limit", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer
    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Sampling");

    //set up groups
    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("dont_skip_align", TOGGLE_LEAVE_ACTIVE, group3);
    {
    GroupContext context (group3);
    place("psi_sampling");
    place("offset_range");
    place("offset_step");

    current_y += STEPY / 2;
    place("allow_coarser");

    }
    guientries["dont_skip_align"].cb_menu_i();  // to make default effective

    }

    {
    GroupContext context (tabcarr[4]);
    tabcarr[4]->label("Helix");
    resetHeight();

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_helix", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);

    place("helical_tube_outer_diameter");
    place("do_bimodal_psi");
    place("range_psi");

    GroupContext (group7 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_restrict_xoff", TOGGLE_LEAVE_ACTIVE, group7);
    {
    GroupContext context (group7);
    place("helical_rise", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["do_restrict_xoff"].cb_menu_i();

    }
    guientries["do_helix"].cb_menu_i();  // to make default effective

    }

    {
    GroupContext context (tabcarr[5]);
    tabcarr[5]->label("Compute");
    resetHeight();

    place("do_parallel_discio");
    place("nr_pool");

    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group5, true);
    {
    GroupContext context (group5);
    place("scratch_dir");
    }

    place("do_combine_thru_disc");

    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    GroupContext (group6 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group6);
    {
    GroupContext context (group6);
    place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["use_gpu"].cb_menu_i();

    }
}

void JobWindow::initialiseInimodelWindow() {

    setupTabs(5);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_img", TOGGLE_DEACTIVATE);
    place("fn_cont", TOGGLE_REACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("CTF");

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    #ifdef ALLOW_CTF_IN_SGD
    place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);

    {
    GroupContext context (group1);
    place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
    }

    guientries["do_ctf_correction"].cb_menu_i();  // To make default effective
    #else
    Fl_Text_Buffer *textbuff1 = new Fl_Text_Buffer();
    textbuff1->text("CTF-modulation, as mentioned in claim 1 of patent US10,282,513B2, is disabled\nYou can enable it by rebuilding, using -DALLOW_CTF_IN_SGD=ON in cmake.");
    Fl_Text_Display* textdisp1 = new Fl_Text_Display(XCOL1, current_y, WCOL1+WCOL2+WCOL3+10, STEPY*1.8);
    textdisp1->textsize(11);
    textdisp1->color(GUI_BACKGROUND_COLOR);
    textdisp1->buffer(textbuff1);

    current_y += STEPY * 2.5;

    place("do_ctf_correction", TOGGLE_ALWAYS_DEACTIVATE);

    {
    GroupContext context (group1);
    place("ctf_phase_flipped", TOGGLE_ALWAYS_DEACTIVATE);
    place("ctf_intact_first_peak", TOGGLE_ALWAYS_DEACTIVATE);
    }
    #endif
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Optimisation");
    resetHeight();

    place("nr_classes", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("particle_diameter");
    place("do_solvent", TOGGLE_DEACTIVATE);
    place("sym_name", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("sampling");
    place("offset_range");
    place("offset_step");
    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("SGD");

    resetHeight();

    place("sgd_ini_iter");
    place("sgd_inbetween_iter");
    place("sgd_fin_iter");
    place("sgd_write_iter");

    current_y += STEPY / 2;  // Add a little spacer

    place("sgd_ini_resol");
    place("sgd_fin_resol");

    current_y += STEPY / 2;  // Add a little spacer

    place("sgd_ini_subset_size");
    place("sgd_fin_subset_size");

    current_y += STEPY / 2;  // Add a little spacer

    place("sgd_sigma2fudge_halflife", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[4]);
    tabcarr[4]->label("Compute");
    resetHeight();

    place("do_parallel_discio");
    place("nr_pool");
    place("do_pad1");
    place("skip_gridding");
    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group5, true);
    {
    GroupContext context (group5);
    place("scratch_dir");
    }

    place("do_combine_thru_disc");

    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    GroupContext (group6 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group6);
    {
    GroupContext context (group6);
    place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["use_gpu"].cb_menu_i();

    }
}

void JobWindow::initialiseClass3DWindow() {

    setupTabs(7);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_img", TOGGLE_DEACTIVATE);
    place("fn_cont", TOGGLE_REACTIVATE);
    place("fn_ref", TOGGLE_DEACTIVATE);
    place("fn_mask");
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Reference");
    resetHeight();

    place("ref_correct_greyscale", TOGGLE_DEACTIVATE);
    place("ini_high", TOGGLE_DEACTIVATE);
    current_y += STEPY / 2;  // Add a little spacer

    place("sym_name", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("CTF");

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);

    place("ctf_corrected_ref", TOGGLE_DEACTIVATE);
    place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);

    }
    guientries["do_ctf_correction"].cb_menu_i();  // To make default effective

    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Optimisation");
    resetHeight();

    //set up groups
    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("nr_classes", TOGGLE_DEACTIVATE);
    place("tau_fudge");

    current_y += STEPY / 2;  // Add a little spacer

    place("nr_iter");
    place("do_fast_subsets", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("particle_diameter");
    place("do_zero_mask", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("highres_limit", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[4]);
    tabcarr[4]->label("Sampling");

    //set up groups
    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("dont_skip_align", TOGGLE_LEAVE_ACTIVE, group3);
    {
    GroupContext context (group3);

    place("sampling");
    place("offset_range");
    place("offset_step");

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_local_ang_searches", TOGGLE_LEAVE_ACTIVE, group4);
    {
    GroupContext context (group4);
    place("sigma_angles");
    place("relax_sym");
    }
    guientries["do_local_ang_searches"].cb_menu_i();  // to make default effective

    current_y += STEPY / 2;
    place("allow_coarser");
    }
    guientries["dont_skip_align"].cb_menu_i();  // to make default effective

    }

    {
    GroupContext context (tabcarr[5]);
    tabcarr[5]->label("Helix");
    resetHeight();
    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    // helix_text", TOGGLE_DEACTIVATE);  //(current_y, "Nov 21, 2015");

    place("do_helix", TOGGLE_DEACTIVATE, group5);
    {
    GroupContext context (group5);
    place2("helical_tube_inner_diameter", "helical_tube_outer_diameter", "Tube diameter - inner, outer (A):", TOGGLE_DEACTIVATE);
    place3("range_rot", "range_tilt", "range_psi", "Angular search range - rot, tilt, psi (deg):", TOGGLE_DEACTIVATE);
    place("helical_range_distance", TOGGLE_DEACTIVATE);
    place("keep_tilt_prior_fixed", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group8 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_apply_helical_symmetry", TOGGLE_DEACTIVATE, group8);
    {
    GroupContext context (group8);
    place("helical_nr_asu", TOGGLE_DEACTIVATE);
    place2("helical_twist_initial", "helical_rise_initial", "Initial twist (deg), rise (A):", TOGGLE_DEACTIVATE);
    place("helical_z_percentage", TOGGLE_DEACTIVATE);
    }
    guientries["do_apply_helical_symmetry"].cb_menu_i();  // to make default effective

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group6 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_local_search_helical_symmetry", TOGGLE_DEACTIVATE, group6);
    {
    GroupContext context (group6);
    place3("helical_twist_min", "helical_twist_max", "helical_twist_inistep", "Twist search - Min, Max, Step (deg):", TOGGLE_DEACTIVATE);
    place3("helical_rise_min", "helical_rise_max", "helical_rise_inistep", "Rise search - Min, Max, Step (A):", TOGGLE_DEACTIVATE);
    }
    guientries["do_local_search_helical_symmetry"].cb_menu_i();  // to make default effective

    }
    guientries["do_helix"].cb_menu_i();  // to make default effective

    }

    {
    GroupContext context (tabcarr[6]);
    tabcarr[6]->label("Compute");
    resetHeight();

    place("do_parallel_discio");
    place("nr_pool");
    place("do_pad1");
    place("skip_gridding");
    GroupContext (group7 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group7, true);
    {
    GroupContext context (group7);
    place("scratch_dir");
    }

    place("do_combine_thru_disc");
    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    GroupContext (group8 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group8);
    {
    GroupContext context (group8);
    place("gpu_ids");
    }
    guientries["use_gpu"].cb_menu_i();  // This is to make the default effective

    }
}

void JobWindow::initialiseAutorefineWindow() {

    setupTabs(7);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_img", TOGGLE_DEACTIVATE);
    place("fn_cont", TOGGLE_REACTIVATE);
    place("fn_ref", TOGGLE_DEACTIVATE);
    place("fn_mask");
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Reference");
    resetHeight();

    place("ref_correct_greyscale", TOGGLE_DEACTIVATE);
    place("ini_high", TOGGLE_DEACTIVATE);
    current_y += STEPY / 2;  // Add a little spacer
    place("sym_name", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("CTF");

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    resetHeight();

    place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("ctf_corrected_ref", TOGGLE_DEACTIVATE);
    place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
    }
    guientries["do_ctf_correction"].cb_menu_i();  // To make default effective

    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Optimisation");
    resetHeight();

    place("particle_diameter");
    place("do_zero_mask", TOGGLE_DEACTIVATE);
    current_y += STEPY / 2;  // Add a little spacer

    place("do_solvent_fsc");
    }

    {
    GroupContext context (tabcarr[4]);
    tabcarr[4]->label("Auto-sampling");
    resetHeight();

    place("sampling", TOGGLE_DEACTIVATE);
    place("offset_range", TOGGLE_DEACTIVATE);
    place("offset_step", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;
    place("auto_local_sampling", TOGGLE_DEACTIVATE);
    place("relax_sym");
    current_y += STEPY / 2;
    place("auto_faster");
    }

    {
    GroupContext context (tabcarr[5]);
    tabcarr[5]->label("Helix");
    resetHeight();
    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_helix", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place2("helical_tube_inner_diameter", "helical_tube_outer_diameter", "Tube diameter - inner, outer (A):",TOGGLE_DEACTIVATE);
    place3("range_rot", "range_tilt", "range_psi", "Angular search range - rot, tilt, psi (deg):", TOGGLE_DEACTIVATE);
    place("helical_range_distance", TOGGLE_DEACTIVATE);
    place("keep_tilt_prior_fixed", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_apply_helical_symmetry", TOGGLE_DEACTIVATE, group5);
    {
    GroupContext context (group5);
    place("helical_nr_asu", TOGGLE_DEACTIVATE);
    place2("helical_twist_initial", "helical_rise_initial", "Initial twist (deg), rise (A):",TOGGLE_DEACTIVATE);
    place("helical_z_percentage", TOGGLE_DEACTIVATE);
    }
    guientries["do_apply_helical_symmetry"].cb_menu_i();  // to make default effective

    current_y += STEPY / 2;  // Add a little spacer

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_local_search_helical_symmetry", TOGGLE_DEACTIVATE, group3);
    {
    GroupContext context (group3);
    place3("helical_twist_min", "helical_twist_max", "helical_twist_inistep", "Twist search - Min, Max, Step (deg):", TOGGLE_DEACTIVATE);
    place3("helical_rise_min", "helical_rise_max", "helical_rise_inistep", "Rise search - Min, Max, Step (A):", TOGGLE_DEACTIVATE);
    }
    guientries["do_local_search_helical_symmetry"].cb_menu_i();  // to make default effective

    }
    guientries["do_helix"].cb_menu_i();  // to make default effective

    }

    {
    GroupContext context (tabcarr[6]);
    tabcarr[6]->label("Compute");
    resetHeight();

    place("do_parallel_discio");
    place("nr_pool");
    place("do_pad1");
    place("skip_gridding");

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group4, true);
    {
    GroupContext context (group4);
    place("scratch_dir");
    }

    place("do_combine_thru_disc");

    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));
    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group5);
    {
    GroupContext context (group5);
    place("gpu_ids");
    }
    guientries["use_gpu"].cb_menu_i();  // This is to make the default effective
    }
}

void JobWindow::initialiseMultiBodyWindow() {

    setupTabs(4);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_in", TOGGLE_DEACTIVATE);
    place("fn_cont", TOGGLE_REACTIVATE);
    place("fn_bodies", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;  // Add a little spacer

    place("do_subtracted_bodies", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Auto-sampling");
    resetHeight();

    place("sampling", TOGGLE_DEACTIVATE);
    place("offset_range", TOGGLE_DEACTIVATE);
    place("offset_step", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Analyse");
    resetHeight();

    GroupContext (group5 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_analyse", TOGGLE_LEAVE_ACTIVE, group5);
    {
    GroupContext context (group5);
    place("nr_movies");

    GroupContext (group6 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_select", TOGGLE_LEAVE_ACTIVE, group6);
    {
    GroupContext context (group6);
    place("select_eigenval");
    place("eigenval_min");
    place("eigenval_max");
    }
    guientries["do_select"].cb_menu_i();  // This is to make the default effective

    }
    guientries["do_analyse"].cb_menu_i();  // This is to make the default effective

    }

    {
    GroupContext context (tabcarr[3]);
    tabcarr[3]->label("Compute");
    resetHeight();

    place("do_parallel_discio");
    place("nr_pool");
    place("do_pad1");
    place("skip_gridding");

    GroupContext (group7 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group7, true);
    {
    GroupContext context (group7);
    place("scratch_dir");
    }

    place("do_combine_thru_disc");

    current_y += STEPY / 2;  // Add a little spacer

    // Set up queue groups for running tab
    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));
    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group4);
    {
    GroupContext context (group4);
    place("gpu_ids");
    }
    guientries["use_gpu"].cb_menu_i();  // This is to make the default effective

    }
}

void JobWindow::initialiseMaskcreateWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_in", TOGGLE_DEACTIVATE);  //(current_y, "Input 3D map:", NODE::3DREF, "", "MRC map files (*.mrc)", "Provide an input MRC map from which to start binarizing the map.");
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Mask");
    resetHeight();

    place("lowpass_filter");
    place("angpix");

    current_y += STEPY / 2;  // Add a little spacer

    place("inimask_threshold");
    place("extend_inimask");
    place("width_mask_edge");
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Helix");
    resetHeight();

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_helix", TOGGLE_LEAVE_ACTIVE, group1);
    {
    GroupContext context (group1);
    place("helical_z_percentage");
    }
    guientries["do_helix"].cb_menu_i();  // to make default effective

    }
}


void JobWindow::initialiseJoinstarWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("Particles");
    resetHeight();

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_part", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("fn_part1", TOGGLE_DEACTIVATE);
    place("fn_part2", TOGGLE_DEACTIVATE);
    place("fn_part3", TOGGLE_DEACTIVATE);
    place("fn_part4", TOGGLE_DEACTIVATE);
    }
    guientries["do_part"].cb_menu_i();  // make default active

    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Micrographs");
    resetHeight();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_mic", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place("fn_mic1", TOGGLE_DEACTIVATE);
    place("fn_mic2", TOGGLE_DEACTIVATE);
    place("fn_mic3", TOGGLE_DEACTIVATE);
    place("fn_mic4", TOGGLE_DEACTIVATE);
    }
    guientries["do_mic"].cb_menu_i();  // make default active

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Movies");
    resetHeight();

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_mov", TOGGLE_DEACTIVATE, group3);  //(current_y, "Combine movie STAR files?", false, "", mov_group);
    {
    GroupContext context (group3);
    place("fn_mov1", TOGGLE_DEACTIVATE);
    place("fn_mov2", TOGGLE_DEACTIVATE);
    place("fn_mov3", TOGGLE_DEACTIVATE);
    place("fn_mov4", TOGGLE_DEACTIVATE);
    }
    guientries["do_mov"].cb_menu_i();  // make default active

    }
}

void JobWindow::initialiseSubtractWindow() {

    setupTabs(2);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_opt", TOGGLE_DEACTIVATE);
    place("fn_mask", TOGGLE_DEACTIVATE);

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_data", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("fn_data", TOGGLE_DEACTIVATE);
    }
    guientries["do_data"].cb_menu_i();  // make default active

    current_y += STEPY / 2;

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_fliplabel", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    place("fn_fliplabel", TOGGLE_DEACTIVATE);
    }
    guientries["do_fliplabel"].cb_menu_i();  // make default active

    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Centering");
    resetHeight();

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_center_mask", TOGGLE_DEACTIVATE, group3, true);
    {
    GroupContext context (group3);

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_center_xyz", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);
    place3("center_x", "center_y", "center_z", "Center coordinate - X, Y, Z (pix):", TOGGLE_DEACTIVATE);
    }
    guientries["do_center_xyz"].cb_menu_i();  // To make default effective

    }
    guientries["do_center_mask"].cb_menu_i();  // To make default effective

    current_y += STEPY / 2;

    place("new_box", TOGGLE_DEACTIVATE);
    }
}

void JobWindow::initialisePostprocessWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();
    place("fn_in", TOGGLE_DEACTIVATE);  //(current_y, "One of the 2 unfiltered half-maps:", NODE::HALFMAP, "", "MRC map files (*half1_class001_unfil.mrc)", "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
    place("fn_mask", TOGGLE_DEACTIVATE);  //(current_y, "Solvent mask:", NODE::MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");

    current_y += STEPY / 2;

    place("angpix");
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Sharpen");
    resetHeight();

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_auto_bfac", TOGGLE_LEAVE_ACTIVE, group1);
    {
    GroupContext context (group1);
    place("autob_lowres");
    }
    guientries["do_auto_bfac"].cb_menu_i();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_adhoc_bfac", TOGGLE_LEAVE_ACTIVE, group2);
    {
    GroupContext context (group2);
    place("adhoc_bfac");
    }
    guientries["do_adhoc_bfac"].cb_menu_i();

    current_y += STEPY / 2;

    place("fn_mtf");
    place("mtf_angpix");
    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Filter");
    resetHeight();

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_skip_fsc_weighting", TOGGLE_LEAVE_ACTIVE, group3);
    {
    GroupContext context (group3);
    place("low_pass");
    }
    guientries["do_skip_fsc_weighting"].cb_menu_i();

    }
}

void JobWindow::initialiseLocresWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    place("fn_in", TOGGLE_DEACTIVATE);
    place("fn_mask");

    current_y += STEPY / 2;

    place("angpix", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("ResMap");
    resetHeight();

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_resmap_locres", TOGGLE_DEACTIVATE, group1);
    {
    GroupContext context (group1);
    place("fn_resmap", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;

    place("pval", TOGGLE_DEACTIVATE);
    place("minres", TOGGLE_DEACTIVATE);
    place("maxres", TOGGLE_DEACTIVATE);
    place("stepres", TOGGLE_DEACTIVATE);
    }
    guientries["do_resmap_locres"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Relion");
    resetHeight();

    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_relion_locres", TOGGLE_DEACTIVATE, group2);
    {
    GroupContext context (group2);
    // place("locres_sampling", TOGGLE_DEACTIVATE);
    // place("randomize_at", TOGGLE_DEACTIVATE);
    // current_y += STEPY / 2;
    place("adhoc_bfac", TOGGLE_DEACTIVATE);
    place("fn_mtf", TOGGLE_DEACTIVATE);
    }
    guientries["do_relion_locres"].cb_menu_i();

    }
}

void JobWindow::initialiseMotionrefineWindow() {

    setupTabs(3);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    // I/O
    place("fn_mic", TOGGLE_DEACTIVATE);
    place("fn_data", TOGGLE_DEACTIVATE);
    place("fn_post", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;

    place("first_frame", TOGGLE_DEACTIVATE);
    place("last_frame", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;

    place("extract_size", TOGGLE_DEACTIVATE);
    place("rescale", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Train");
    resetHeight();

    // Train for optimal parameters
    GroupContext (group2 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_param_optim", TOGGLE_LEAVE_ACTIVE, group2);
    {
    GroupContext context (group2);
    place("eval_frac");
    place("optim_min_part");
    }
    guientries["do_param_optim"].cb_menu_i();

    }

    {
    GroupContext context (tabcarr[2]);
    tabcarr[2]->label("Polish");
    resetHeight();

    // Polishing
    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));
    place("do_polish", TOGGLE_DEACTIVATE, group1);

    current_y += STEPY / 2;

    {
    GroupContext context (group1);
    place("opt_params", TOGGLE_DEACTIVATE);

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_own_params", TOGGLE_DEACTIVATE, group4);
    {
    GroupContext context (group4);
    place("sigma_vel", TOGGLE_DEACTIVATE);
    place("sigma_div", TOGGLE_DEACTIVATE);
    place("sigma_acc", TOGGLE_DEACTIVATE);
    }
    guientries["do_own_params"].cb_menu_i();

    current_y += STEPY / 2;

    place("minres", TOGGLE_DEACTIVATE);
    place("maxres", TOGGLE_DEACTIVATE);
    }  // Is this the right place to call group1->end()?
    }
}

void JobWindow::initialiseCtfrefineWindow() {

    setupTabs(2);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("I/O");
    resetHeight();

    // I/O
    place("fn_data", TOGGLE_DEACTIVATE);
    place("fn_post", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Fit");
    resetHeight();

    GroupContext (group3 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));
    place("do_aniso_mag", TOGGLE_LEAVE_ACTIVE, group3, true);  //true means: activating aniso_mag will deactive higher-order aberrations

    current_y += STEPY / 2;

    {
    GroupContext context (group3);

    GroupContext (group1 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_ctf", TOGGLE_LEAVE_ACTIVE, group1);
    {
    GroupContext context (group1);
    place("do_defocus", TOGGLE_LEAVE_ACTIVE);
    place("do_astig", TOGGLE_LEAVE_ACTIVE);
    place("do_bfactor", TOGGLE_LEAVE_ACTIVE);
    place("do_phase", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["do_ctf"].cb_menu_i();

    current_y += STEPY / 2;

    GroupContext (group4 = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600 - MENUHEIGHT, ""));

    place("do_tilt", TOGGLE_LEAVE_ACTIVE, group4);
    {
    GroupContext context (group4);
    place("do_trefoil", TOGGLE_LEAVE_ACTIVE);
    }
    guientries["do_tilt"].cb_menu_i();

    current_y += STEPY / 2;

    place("do_4thorder", TOGGLE_LEAVE_ACTIVE);

    }

    guientries["do_aniso_mag"].cb_menu_i();
    current_y += STEPY / 2;

    place("minres", TOGGLE_DEACTIVATE);
    }
}

void JobWindow::initialiseExternalWindow() {

    setupTabs(2);

    {
    GroupContext context (tabcarr[0]);
    tabcarr[0]->label("Input");
    resetHeight();

    // I/O
    place("fn_exe", TOGGLE_DEACTIVATE);

    current_y += STEPY / 2;
    place("in_mov", TOGGLE_DEACTIVATE);
    place("in_mic", TOGGLE_DEACTIVATE);
    place("in_part", TOGGLE_DEACTIVATE);
    place("in_coords", TOGGLE_DEACTIVATE);
    place("in_3dref", TOGGLE_DEACTIVATE);
    place("in_mask", TOGGLE_DEACTIVATE);
    }

    {
    GroupContext context (tabcarr[1]);
    tabcarr[1]->label("Params");
    resetHeight();

    place2("param1_label", "param1_value", "Param1 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param2_label", "param2_value", "Param2 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param3_label", "param3_value", "Param3 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param4_label", "param4_value", "Param4 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param5_label", "param5_value", "Param5 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param6_label", "param6_value", "Param6 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param7_label", "param7_value", "Param7 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param8_label", "param8_value", "Param8 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param9_label", "param9_value", "Param9 label, value:", TOGGLE_LEAVE_ACTIVE);
    place2("param10_label", "param10_value", "Param10 label, value:", TOGGLE_LEAVE_ACTIVE);
    }
}
