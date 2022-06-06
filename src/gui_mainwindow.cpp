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

#include "src/gui_mainwindow.h"
#include "src/gui_background.xpm"

bool show_scheduler;
bool show_expand_stdout;

// The StdOutDisplay allows looking at the entire stdout or stderr file
int StdOutDisplay::handle(int ev) {

    if (ev == FL_PUSH && Fl::event_clicks()) {
        // double-click
        if (Fl::event_clicks()) {
            if (show_scheduler) {
                current_browse_directory = schedule.name;
            } else {
                if (current_job < 0) return 0;
                current_browse_directory = pipeline.processList[current_job].name;
            }
            FileName fn = current_browse_directory + fn_file;
            if (exists(fn)) {
                FileName fn_note;
                if (fn_file == "run.out") {
                    if (maingui_do_read_only) {
                        fn_note = fn.c_str();
                    } else {
                        // temp file
                        fn_note = ".gui_tmpstd";
                        std::string command = "awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' < " + fn + " > .gui_tmpstd";
                        system(command.c_str());
                    }
                } else {
                    fn_note = fn;
                }
                NoteEditorWindow *w = new NoteEditorWindow(800, 400, fn.c_str(), fn_note, !maingui_do_read_only);
                // allow_save (useful to remove past errors) unless maingui_do_read_only
                w->show();
                return 1;
            }
        }
    }
    return 0;
}

int SchedulerWindow::fill(FileName _pipeline_name, std::vector<FileName> _scheduled_jobs) {
    // color(GUI_BACKGROUND_COLOR);
    int y, ymax, ystep;
    y = 2; ymax = 2; ystep = 35;
    int xcol = w() - 120;

    // Scroll bars
    Fl_Scroll scroll(0, y, w(), h());
    scroll.type(Fl_Scroll::VERTICAL);

    my_jobs.clear();
    pipeline_name = _pipeline_name;
    for (int i = 0; i < _scheduled_jobs.size(); i++) {
        my_jobs.push_back(_scheduled_jobs[i]);
        int xcoor = i < 1 + _scheduled_jobs.size() / 2 ? 20 : w() - 170;
        if (i == 1 + _scheduled_jobs.size() / 2)
            y = 2;
        Fl_Check_Button *mycheck = new Fl_Check_Button(xcoor, y, ystep - 8, ystep - 8, _scheduled_jobs[i].c_str());
        mycheck->labelsize(ENTRY_FONTSIZE);
        check_buttons.push_back(mycheck);
        mycheck->value(1);
        y += ystep;
        if (y > ymax) { ymax = y; }
    }

    y = ymax;
    schedule_name = new Fl_Input(xcol, y, 100, ystep - 8, "Provide a name for this schedule: ");
    y += ystep;
    wait_before = new Fl_Input(xcol, y, 100, ystep - 8, "Wait this many minutes before starting?");
    y += ystep;
    repeat = new Fl_Input(xcol, y, 100, ystep - 8, "Run the jobs how many times?");
    y += ystep;
    wait = new Fl_Input(xcol, y, 100, ystep - 8, "Wait at least in between (in minutes)?");
    y += ystep;
    wait_after = new Fl_Input(xcol, y, 100, ystep - 8, "Wait at least after each job (in seconds)?");
    y += ystep;

    // Set the input value
    schedule_name->value("schedule1");
    schedule_name->color(GUI_INPUT_COLOR);
    schedule_name->textsize(ENTRY_FONTSIZE);
    schedule_name->labelsize(ENTRY_FONTSIZE);
    repeat->value("1");
    repeat->color(GUI_INPUT_COLOR);
    repeat->textsize(ENTRY_FONTSIZE);
    repeat->labelsize(ENTRY_FONTSIZE);
    wait->value("15");
    wait->color(GUI_INPUT_COLOR);
    wait->textsize(ENTRY_FONTSIZE);
    wait->labelsize(ENTRY_FONTSIZE);
    wait_before->value("0");
    wait_before->color(GUI_INPUT_COLOR);
    wait_before->textsize(ENTRY_FONTSIZE);
    wait_before->labelsize(ENTRY_FONTSIZE);
    wait_after->value("10");
    wait_after->color(GUI_INPUT_COLOR);
    wait_after->textsize(ENTRY_FONTSIZE);
    wait_after->labelsize(ENTRY_FONTSIZE);

    // Button to execute
    Fl_Button *execute_button = new Fl_Button(w() - 200, y, 80, 30, "Execute");
    execute_button->color(GUI_RUNBUTTON_COLOR);
    execute_button->labelsize(12);
    execute_button->callback(cb_execute, this);

    // Button to cancel
    Fl_Button *cancel_button = new Fl_Button(w() - 100, y, 80, 30, "Cancel");
    cancel_button->color(GUI_RUNBUTTON_COLOR);
    cancel_button->labelsize(12);
    cancel_button->callback(cb_cancel, this);

    resizable(*this);
    show();

    return Fl::run();
}

void SchedulerWindow::cb_cancel(Fl_Widget*, void* v) {
    SchedulerWindow* T = (SchedulerWindow*)v;
    T->hide();
}

void SchedulerWindow::cb_execute(Fl_Widget*, void* v) {
    SchedulerWindow* T = (SchedulerWindow*)v;
    T->cb_execute_i();
    T->hide();
}

void SchedulerWindow::cb_execute_i() {
    FileName fn_sched(schedule_name->value());
    FileName fn_check = "RUNNING_PIPELINER_" + pipeline_name + "_" + fn_sched;
    if (exists(fn_check)) {
        fl_message("%s", ((std::string) "ERROR: a file called " + fn_check + " already exists. \n This implies another set of scheduled jobs with this name is already running. \n Cancelling job execution...").c_str());
    } else {
        // Make a string with all job-ids to process
        std::string jobids = "";
        for (int i = 0; i < my_jobs.size(); i++) {
            if (check_buttons[i]->value())
                jobids += my_jobs[i] + " ";
        }
        jobids = "\"" + jobids + "\"";

        std::string command = (std::string) "relion_pipeliner"
            + " --pipeline " + pipeline_name
            + " --schedule " + fn_sched
            + " --repeat " + repeat->value()
            + " --min_wait " + wait->value()
            + " --min_wait_before " + wait_before->value()
            + " --sec_wait_after " + wait_after->value()
            + " --RunJobs " + jobids
            + " &";  // Run this in the background, so control returns to the window
        system(command.c_str());
        std::cout << " Launching: " << command << std::endl;
        std::cout << " Stop execution of this set of scheduled jobs by deleting file: " << fn_check << std::endl;
    }
}

NoteEditorWindow::NoteEditorWindow(int w, int h, const char* title, FileName _fn_note, bool _allow_save): Fl_Window(w, h, title) {
    allow_save = _allow_save;
    editor = new Fl_Text_Editor(0, 0, w, h - 50);
    editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,10);
    textbuff_note = new Fl_Text_Buffer;
    editor->buffer(textbuff_note);
    textbuff_note->transcoding_warning_action = nullptr;
    fn_note = _fn_note;
    if (exists(fn_note)) {
        int err = textbuff_note->loadfile(fn_note.c_str());
    } else {
        textbuff_note->text("Describe what this job or project is about here...");
    }
    editor->insert_position(editor->buffer()->length());
    editor->show_insert_position();

    if (allow_save) {
        // Button to save and exit
        Fl_Button *save_button = new Fl_Button(w - 200, h - 40, 80, 30, "Save");
        save_button->color(GUI_RUNBUTTON_COLOR);
        save_button->labelsize(12);
        save_button->callback(cb_save, this);
    }

    // Button to exit
    Fl_Button *cancel_button = new Fl_Button(w - 100, h - 40, 80, 30, "Cancel");
    cancel_button->color(GUI_RUNBUTTON_COLOR);
    cancel_button->labelsize(12);
    cancel_button->callback(cb_cancel, this);
    resizable(*this);
}

template <typename T>
T *creatio_ex_nihilo(void *v) {
    return (T*) v;
}

void NoteEditorWindow::cb_cancel(Fl_Widget*, void* v) {
    NoteEditorWindow *T = (NoteEditorWindow*) v;
    T->hide();
}

void NoteEditorWindow::cb_save(Fl_Widget*, void* v) {
    NoteEditorWindow *T = (NoteEditorWindow*) v;
    T->save();
    T->hide();
}

void NoteEditorWindow::save() {
    textbuff_note->savefile(fn_note.c_str());
}

GuiMainWindow::GuiMainWindow(
    int w, int h, const char* title,
    FileName fn_pipe, FileName fn_sched,
    int _update_every_sec, int _exit_after_sec, bool _do_read_only
): Fl_Window(w, h, title) {
    // Set initial Timer
    tickTimeLastChanged();

    show_expand_stdout = false;

    // Setup read_only
    maingui_do_read_only = _do_read_only;
    pipeline.do_read_only = _do_read_only;
    do_order_alphabetically = false;

    FileName fn_lock = ".gui_projectdir";
    if (!exists(fn_lock)) {
        int ret = fl_choice("Your current directory does not look like a RELION project directory.\nOnly run the RELION GUI from your project directory.\nDo you want to start a new project here?", "No", "Yes", 0);
        this->begin();  // apparently fl_choice changes Fl_Group::current. Thus, we have to reclaim it.
        if (ret != 1) {
            std::cout << " Exiting ... " << std::endl;
            exit(0);
        }
        touch(".gui_projectdir");
    }

    // First set up the old part of the GUI
    h = GUIHEIGHT_OLD;

    {
    /// TODO: control file location and use better figure
    GroupContext context (background_grp = new Fl_Group(WCOL0 - 10, 0, w - WCOL0, h - 55));

    // Initial screen picture with some explanation on how to use the GUI
    // image_box = new Fl_Box(WCOL0 - 8, 0, w - WCOL0, h - 35);  // widget that will contain image
    image_box = new Fl_Box(WCOL0 - 8, 50, w - WCOL0, h - 120);  // widget that will contain image
    xpm_image = new Fl_Pixmap(gui_background);
    image_box->image(xpm_image);  // attach xpm image to box
    }

    // Read in schedule, if it exists.
    // Otherwise, just initialise schedule with its name.
    if (fn_sched != "") {
        show_scheduler = true;
        create_scheduler_gui = true;
        schedule.do_read_only = _do_read_only;
        schedule.setName(fn_sched + "/");
        pipeline.name = fn_sched + "/schedule";
        if (exists(schedule.name + "schedule.star")) {
            schedule.read(DONT_LOCK);
            pipeline.name = fn_sched + "/schedule";
        } else {
            std::string command = "mkdir -p " + fn_sched;
            system(command.c_str());
            schedule.write(DONT_LOCK);  // empty write
        }
    } else {
        show_scheduler = false;
        create_scheduler_gui = false;
        // Read in the pipeline STAR file if it exists
        pipeline.name = fn_pipe;
    }
    if (exists(pipeline.name + "_pipeline.star")) {
        PipeLine::rwlock (pipeline, "mainGUI constructor");
        // With the locking system, each read needs to be followed soon with a write
    } else {
        pipeline.write();
    }

    color(GUI_BACKGROUND_COLOR);
    menubar = new Fl_Menu_Bar(-3, 0, WCOL0 - 7, MENUHEIGHT);
    menubar->add("File/Re-read pipeline", FL_ALT + 'r', cb_reread_pipeline, this);
    menubar->add("File/Edit project note", FL_ALT + 'e', cb_edit_project_note, this);
    if (!maingui_do_read_only)
        menubar->add("File/Print all notes",  0, cb_print_notes, this);
    if (!maingui_do_read_only)
        menubar->add("File/Remake .Nodes\\/", FL_ALT + 'n', cb_remake_nodesdir, this);
    menubar->add("File/Display", FL_ALT + 'd', cb_display, this);
    menubar->add("File/_Overwrite continue", FL_ALT + 'o', cb_toggle_overwrite_continue, this);
    menubar->add("File/_Show initial screen", FL_ALT + 'z', cb_show_initial_screen, this);
    if (!maingui_do_read_only)
        menubar->add("File/_Empty trash", FL_ALT + 't', cb_empty_trash, this);
    menubar->add("File/About", 0, cb_about, this);
    menubar->add("File/Quit", FL_ALT + 'q', cb_quit, this);
    if (!maingui_do_read_only) {
        menubar->add("Jobs/Save job settings", FL_ALT + 's', cb_save, this);
        menubar->add("Jobs/_Load job settings", FL_ALT + 'l', cb_load, this);
    }
    menubar->add("Jobs/Order alphabetically", FL_ALT + 'a', cb_order_jobs_alphabetically, this);
    menubar->add("Jobs/_Order chronologically", FL_ALT + 'c', cb_order_jobs_chronologically, this);
    if (!maingui_do_read_only) {
        menubar->add("Jobs/_Undelete job(s)", FL_ALT + 'u', cb_undelete_job, this);
        menubar->add("Jobs/Run scheduled jobs", 0, cb_start_pipeliner, this);
        menubar->add("Jobs/Stop running scheduled jobs", 0, cb_stop_pipeliner, this);
        menubar->add("Jobs/Export scheduled job(s)",  0, cb_export_jobs, this);
        menubar->add("Jobs/_Import scheduled job(s)",  0, cb_import_jobs, this);
        menubar->add("Jobs/Gently clean all jobs", FL_ALT + 'g', cb_gently_clean_all_jobs, this);
        menubar->add("Jobs/Harshly clean all jobs", FL_ALT + 'h', cb_harshly_clean_all_jobs, this);

        // See which schedules there are
        FileName schedule_wildcard = "Schedules/*";
        std::vector<FileName> schedules;
        schedule_wildcard.globFiles(schedules);
        for (const auto &schedule : schedules) {
            menubar->add(("Schedules/" + schedule).c_str(), 0, cb_toggle_schedule, this);
        }
        menubar->add("Schedules/_Copy schedule", 0, cb_copy_schedule, this);
        menubar->add("Schedules/_Show pipeline", FL_ALT + 'p', cb_toggle_pipeline, this);
    }
    current_y = MENUHEIGHT + 10;

    // Fill browser in the right order
    browser = new Fl_Hold_Browser(10, MENUHEIGHT + 5, WCOL0 - 20, h - MENUHEIGHT - 60);
    {
    browser->textsize(RLN_FONTSIZE - 1);
    current_job = -1;

    std::array<std::pair<const char*, Process::Type>, NR_BROWSE_TABS> specs {{
        {"Import",               Process::IMPORT},
        {"Motion correction",    Process::MOTIONCORR},
        {"CTF estimation",       Process::CTFFIND},
        {"Manual picking",       Process::MANUALPICK},
        {"Auto-picking",         Process::AUTOPICK},
        {"Particle extraction",  Process::EXTRACT},
        {"Subset selection",     Process::CLASSSELECT},
        {"2D classification",    Process::CLASS2D},
        {"3D initial model",     Process::INIMODEL},
        {"3D classification",    Process::CLASS3D},
        {"3D auto-refine",       Process::AUTO3D},
        {"3D multi-body",        Process::MULTIBODY},
        {"CTF refinement",       Process::CTFREFINE},
        {"Bayesian polishing",   Process::MOTIONREFINE},
        {"Mask creation",        Process::MASKCREATE},
        {"Join star files",      Process::JOINSTAR},
        {"Particle subtraction", Process::SUBTRACT},
        {"Post-processing",      Process::POST},
        {"Local resolution",     Process::RESMAP},
        {"External",             Process::EXTERNAL}
    }};

    for (int i = 0; i < specs.size(); ++i) {
        GroupContext context (browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615 - MENUHEIGHT));
        browser->add(specs[i].first);
        gui_jobwindows[i] = new JobWindow();
        gui_jobwindows[i]->initialise(specs[i].second);
    }

    browser->callback(cb_select_browsegroup, this);
    browser->end();
    browser->select(1);  // just start from the beginning
    }

    // Add run buttons on the menubar as well

    print_CL_button = new Fl_Button(GUIWIDTH - 215, h - 90, 100, 32, "Check command");
    print_CL_button->color(GUI_RUNBUTTON_COLOR);
    print_CL_button->labelsize(11);
    print_CL_button->callback(cb_print_cl, this);

    expand_stdout_button = new Fl_Button(XJOBCOL1, GUIHEIGHT_EXT_START, 85, MENUHEIGHT, "I/O view");
    expand_stdout_button->color(GUI_BUTTON_COLOR);
    expand_stdout_button->callback(cb_toggle_expand_stdout, this);

    // A) Pipeliner part of the GUI
    pipeliner_grp = new Fl_Group(0, 0, 2 * w, 2 * h);
    {
    GroupContext context (pipeliner_grp);

    run_button = new Fl_Button(GUIWIDTH - 110 , h - 90, 100, 32, "Run!");
    run_button->color(GUI_RUNBUTTON_COLOR);
    run_button->labelfont(FL_ITALIC);
    run_button->labelsize(14);
    run_button->callback(cb_run, this);
    if (maingui_do_read_only)
        run_button->deactivate();

    schedule_button = new Fl_Button(GUIWIDTH - 320 , h - 90, 100, 32, "Schedule");
    schedule_button->color(GUI_RUNBUTTON_COLOR);
    schedule_button->labelfont(FL_ITALIC);
    schedule_button->labelsize(14);
    schedule_button->callback(cb_schedule, this);
    if (maingui_do_read_only)
        schedule_button->deactivate();

    menubar2 = new Fl_Menu_Bar(XJOBCOL1 + 87, GUIHEIGHT_EXT_START, 95, MENUHEIGHT);
    menubar2->color(GUI_BUTTON_COLOR);
    menubar2->add("Job actions/Edit Note", 0, cb_edit_note, this);
    if (!maingui_do_read_only) {
        menubar2->add("Job actions/Alias", 0, cb_set_alias, this);
        menubar2->add("Job actions/Abort running", 0, cb_abort, this);
        menubar2->add("Job actions/Mark as finished", 0, cb_mark_as_finished, this);
        menubar2->add("Job actions/Mark as failed", 0, cb_mark_as_failed, this);
        menubar2->add("Job actions/Make flowchart", 0, cb_make_flowchart, this);
        menubar2->add("Job actions/Gentle clean", 0, cb_gentle_cleanup, this);
        menubar2->add("Job actions/Harsh clean", 0, cb_harsh_cleanup, this);
        menubar2->add("Job actions/Delete", 0, cb_delete, this);
    }

    // Fl_input with the alias of the new job (or the name of an existing one)
    alias_current_job = new Fl_Input(XJOBCOL2 , GUIHEIGHT_EXT_START + 3, JOBCOLWIDTH, MENUHEIGHT - 6, "Current:");

    // Left-hand side browsers for input/output nodes and processes
    display_io_node = new Fl_Choice(XJOBCOL3 + 50, GUIHEIGHT_EXT_START + 3, 200, MENUHEIGHT - 6);
    display_io_node->label("Display:");
    display_io_node->color(GUI_BUTTON_COLOR);
    display_io_node->callback(cb_display_io_node, this);

    pipeliner_jobs_grp = new Fl_Group(0, 0, 2 * w, 2 * h);
    {
    GroupContext context (pipeliner_jobs_grp);

    // Add browsers for finished and running jobs
    Fl_Text_Buffer *textbuff1 = new Fl_Text_Buffer();
    textbuff1->text("Finished jobs");
    Fl_Text_Display* textdisp1 = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
    textdisp1->buffer(textbuff1);
    textdisp1->color(GUI_BACKGROUND_COLOR);
    finished_job_browser = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHEIGHT + 25);
    finished_job_browser->callback(cb_select_finished_job, this);
    finished_job_browser->textsize(RLN_FONTSIZE - 1);
    finished_job_browser->end();

    Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
    textbuff2->text("Running jobs");
    Fl_Text_Display* textdisp2 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
    textdisp2->buffer(textbuff2);
    textdisp2->color(GUI_BACKGROUND_COLOR);
    running_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    running_job_browser->callback(cb_select_running_job, this);
    running_job_browser->textsize(RLN_FONTSIZE - 1);
    running_job_browser->end();

    Fl_Text_Buffer *textbuff3 = new Fl_Text_Buffer();
    textbuff3->text("Scheduled jobs");
    Fl_Text_Display* textdisp3 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
    textdisp3->buffer(textbuff3);
    textdisp3->color(GUI_BACKGROUND_COLOR);
    scheduled_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    scheduled_job_browser->callback(cb_select_scheduled_job, this);
    scheduled_job_browser->textsize(RLN_FONTSIZE - 1);

    Fl_Text_Buffer *textbuff4 = new Fl_Text_Buffer();
    textbuff4->text("Input to this job");
    Fl_Text_Display* textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
    textdisp4->buffer(textbuff4);
    textdisp4->color(GUI_BACKGROUND_COLOR);
    input_job_browser = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    input_job_browser->callback(cb_select_input_job, this);
    input_job_browser->textsize(RLN_FONTSIZE - 1);

    Fl_Text_Buffer *textbuff5 = new Fl_Text_Buffer();
    textbuff5->text("Output from this job");
    Fl_Text_Display* textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
    textdisp5->buffer(textbuff5);
    textdisp5->color(GUI_BACKGROUND_COLOR);
    output_job_browser = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    output_job_browser->callback(cb_select_output_job, this);
    output_job_browser->textsize(RLN_FONTSIZE - 1);

    // Display stdout and stderr of jobs
    textbuff_stdout = new Fl_Text_Buffer();
    textbuff_stderr = new Fl_Text_Buffer();
    // Disable warning message about UTF-8 transcoding
    textbuff_stdout->transcoding_warning_action = nullptr;
    textbuff_stderr->transcoding_warning_action = nullptr;
    disp_stdout = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDOUT_Y - 5, w - 20, 105);
    disp_stderr = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDERR_Y - 5, w - 20, 50);
    disp_stdout->fn_file = "run.out";
    disp_stderr->fn_file = "run.err";
    textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
    disp_stdout->buffer(textbuff_stdout);
    disp_stderr->buffer(textbuff_stderr);
    disp_stderr->textcolor(FL_RED);
    disp_stdout->textsize(RLN_FONTSIZE - 1);
    disp_stderr->textsize(RLN_FONTSIZE - 1);
    disp_stdout->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
    disp_stderr->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
    disp_stdout->scrollbar_width(0);
    disp_stderr->scrollbar_width(0);
    }

    }

    // B) Scheduler part of the GUI
    scheduler_grp = new Fl_Group(0, 0, 4 * w, 4 * h);
    {
    GroupContext context (scheduler_grp);

    scheduler_run_grp = new Fl_Group(0, 0, 4 * w, 4 * h);
    {
    GroupContext context (scheduler_run_grp);

    // Buttons for current_node and running/aborting the schedule
    scheduler_current_node = new Fl_Choice(XJOBCOL1 + 90 + 65, GUIHEIGHT_EXT_START + 3, 140, 23);
    scheduler_current_node->label("Current:");
    scheduler_current_node->textsize(RLN_FONTSIZE-2);
    scheduler_current_node->color(GUI_INPUT_COLOR);

    scheduler_set_current_button = new Fl_Button(XJOBCOL1 + 90 + 210, GUIHEIGHT_EXT_START + 3, 50, 23);
    scheduler_set_current_button->label("Set");
    scheduler_set_current_button->color(GUI_BUTTON_COLOR);
    scheduler_set_current_button->callback(cb_scheduler_set_current, this);

    scheduler_prev_button = new Fl_Button(XJOBCOL1 + 90 + 210 + 55, GUIHEIGHT_EXT_START + 3, 50, 23);
    scheduler_prev_button->label("Prev");
    scheduler_prev_button->color(GUI_BUTTON_COLOR);
    scheduler_prev_button->callback(cb_scheduler_prev, this);

    scheduler_next_button = new Fl_Button(XJOBCOL1 + 90 + 210 + 2 * 55, GUIHEIGHT_EXT_START + 3, 50, 23);
    scheduler_next_button->label("Next");
    scheduler_next_button->color(GUI_BUTTON_COLOR);
    scheduler_next_button->callback(cb_scheduler_next, this);

    scheduler_reset_button = new Fl_Button(XJOBCOL1 + 90 + 210 + 3 * 55, GUIHEIGHT_EXT_START + 3, 50, 23);
    scheduler_reset_button->label("Reset");
    scheduler_reset_button->color(GUI_BUTTON_COLOR);
    scheduler_reset_button->callback(cb_scheduler_reset, this);

    scheduler_run_button = new Fl_Button(GUIWIDTH - 90, GUIHEIGHT_EXT_START + 1, 80, 25);
    scheduler_run_button->label("Run!");
    scheduler_run_button->color(GUI_RUNBUTTON_COLOR);
    scheduler_run_button->labelfont(FL_ITALIC);
    scheduler_run_button->labelsize(14);
    scheduler_run_button->callback(cb_scheduler_run, this);
    }

    scheduler_unlock_button = new Fl_Button(GUIWIDTH - 256, GUIHEIGHT_EXT_START + 1, 80, 25);
    scheduler_unlock_button->label("Unlock");
    scheduler_unlock_button->labelfont(FL_ITALIC);
    scheduler_unlock_button->labelsize(14);
    scheduler_unlock_button->color(GUI_RUNBUTTON_COLOR);
    scheduler_unlock_button->callback(cb_scheduler_unlock, this);

    // Don't allow any changes on the GUI while a Schedule is running, i.e. its directory is locked for writing
    scheduler_abort_button = new Fl_Button(GUIWIDTH - 173, GUIHEIGHT_EXT_START + 1, 80, 25);
    scheduler_abort_button->label("Abort");
    scheduler_abort_button->labelfont(FL_ITALIC);
    scheduler_abort_button->labelsize(14);
    scheduler_abort_button->color(GUI_RUNBUTTON_COLOR);
    scheduler_abort_button->callback(cb_scheduler_abort, this);

    // scheduler_grp->end();

    scheduler_job_name = new Fl_Input(GUIWIDTH - 550, h - 83, 150, 25, "Name:");
    scheduler_job_name->color(GUI_INPUT_COLOR);

    add_job_button = new Fl_Button(GUIWIDTH - 110 , h - 90, 100, 32, "Add job");
    add_job_button->color(GUI_RUNBUTTON_COLOR);
    add_job_button->labelfont(FL_ITALIC);
    add_job_button->labelsize(14);
    add_job_button->callback(cb_scheduler_add_job, this);

    // Select one of three modes for adding a new job
    scheduler_job_mode = new Fl_Choice(GUIWIDTH - 400 , h - 83, 80, 25);
    scheduler_job_mode->label("");
    scheduler_job_mode->color(GUI_BUTTON_COLOR);
    scheduler_job_mode->textsize(12);
    scheduler_job_mode->menu(job_mode_options);

    scheduler_job_has_started = new Fl_Choice(GUIWIDTH - 320 , h - 83, 100, 25);
    scheduler_job_has_started->label("");
    scheduler_job_has_started->color(GUI_BUTTON_COLOR);
    scheduler_job_has_started->textsize(12);
    scheduler_job_has_started->menu(job_has_started_options);
    /// TODO: fill options for this choice!

    scheduler_jobs_grp = new Fl_Group(0, 0, 4 * w, 4 * h);
    {
    GroupContext context (scheduler_jobs_grp);

    // Scheduler variables
    int height_var = 35;
    Fl_Text_Buffer *textbuffvar = new Fl_Text_Buffer();
    textbuffvar->text("Variables");
    Fl_Text_Display* textdispvar = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START + height_var, JOBCOLWIDTH - 105, 24);
    textdispvar->buffer(textbuffvar);
    textdispvar->textsize(12);
    textdispvar->color(GUI_BACKGROUND_COLOR);
    scheduler_variable_name = new Fl_Input(XJOBCOL1, GUIHEIGHT_EXT_START + height_var + 23, JOBCOLWIDTH * 0.4, 21);
    scheduler_variable_name->color(GUI_INPUT_COLOR);
    scheduler_variable_name->textsize(RLN_FONTSIZE - 2);
    scheduler_variable_value = new Fl_Input(XJOBCOL1 + JOBCOLWIDTH * 0.4, GUIHEIGHT_EXT_START + height_var + 23, JOBCOLWIDTH * 0.6, 21);
    scheduler_variable_value->color(GUI_INPUT_COLOR);
    scheduler_variable_value->textsize(RLN_FONTSIZE-2);
    delete_scheduler_variable_button = new Fl_Button(XJOBCOL1 + JOBCOLWIDTH - 105, GUIHEIGHT_EXT_START + height_var, 50, 23);
    delete_scheduler_variable_button->color(GUI_BUTTON_COLOR);
    delete_scheduler_variable_button->labelfont(FL_ITALIC);
    delete_scheduler_variable_button->labelsize(RLN_FONTSIZE);
    delete_scheduler_variable_button->label("Del");
    delete_scheduler_variable_button->callback(cb_delete_scheduler_variable, this);
    set_scheduler_variable_button = new Fl_Button(XJOBCOL1 + JOBCOLWIDTH - 50, GUIHEIGHT_EXT_START + height_var, 50, 23);
    set_scheduler_variable_button->color(GUI_BUTTON_COLOR);
    set_scheduler_variable_button->labelfont(FL_ITALIC);
    set_scheduler_variable_button->labelsize(RLN_FONTSIZE);
    set_scheduler_variable_button->label("Set");
    set_scheduler_variable_button->callback(cb_set_scheduler_variable, this);
    scheduler_variable_browser = new Fl_Hold_Browser(XJOBCOL1, GUIHEIGHT_EXT_START + height_var + 44, JOBCOLWIDTH, 182);
    scheduler_variable_browser->callback(cb_select_scheduler_variable, this);
    scheduler_variable_browser->textsize(RLN_FONTSIZE - 2);
    scheduler_variable_browser->end();

    // Scheduler operators
    Fl_Text_Buffer *textbuffnode = new Fl_Text_Buffer();
    textbuffnode->text("Operators");
    Fl_Text_Display* textdispnode = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START + height_var, JOBCOLWIDTH - 105, 24);
    textdispnode->buffer(textbuffnode);
    textdispnode->textsize(12);
    textdispnode->color(GUI_BACKGROUND_COLOR);
    scheduler_operator_type = new Fl_Choice(XJOBCOL2, GUIHEIGHT_EXT_START + 23 + height_var, JOBCOLWIDTH / 2 + 10, 21);
    scheduler_operator_type->color(GUI_INPUT_COLOR);
    scheduler_operator_type->menu(operator_type_options);
    scheduler_operator_type->textsize(RLN_FONTSIZE-2);
    scheduler_operator_output = new Fl_Choice(XJOBCOL2 + 34 + JOBCOLWIDTH / 2, GUIHEIGHT_EXT_START + 23 + height_var, JOBCOLWIDTH / 2 - 34, 21);
    scheduler_operator_output->label("->");
    scheduler_operator_output->color(GUI_INPUT_COLOR);
    scheduler_operator_output->textsize(RLN_FONTSIZE-2);
    scheduler_operator_input1 = new Fl_Choice(XJOBCOL2 + 20, GUIHEIGHT_EXT_START + 44 + height_var, JOBCOLWIDTH / 2 - 20, 21);
    scheduler_operator_input1->label("i1:");
    scheduler_operator_input1->color(GUI_INPUT_COLOR);
    scheduler_operator_input1->textsize(RLN_FONTSIZE-2);
    scheduler_operator_input2 = new Fl_Choice(XJOBCOL2 + 34 + JOBCOLWIDTH / 2, GUIHEIGHT_EXT_START + 44 + height_var, JOBCOLWIDTH / 2 - 34, 21);
    scheduler_operator_input2->label("i2:");
    scheduler_operator_input2->textsize(RLN_FONTSIZE-2);
    scheduler_operator_input2->color(GUI_INPUT_COLOR);
    delete_scheduler_operator_button = new Fl_Button(XJOBCOL2 + JOBCOLWIDTH - 105, GUIHEIGHT_EXT_START + height_var, 50, 23);
    delete_scheduler_operator_button->color(GUI_BUTTON_COLOR);
    delete_scheduler_operator_button->labelfont(FL_ITALIC);
    delete_scheduler_operator_button->labelsize(RLN_FONTSIZE);
    delete_scheduler_operator_button->label("Del");
    delete_scheduler_operator_button->callback(cb_delete_scheduler_operator, this);
    add_scheduler_operator_button = new Fl_Button(XJOBCOL2 + JOBCOLWIDTH - 50, GUIHEIGHT_EXT_START + height_var, 50, 23);
    add_scheduler_operator_button->color(GUI_BUTTON_COLOR);
    add_scheduler_operator_button->labelfont(FL_ITALIC);
    add_scheduler_operator_button->labelsize(RLN_FONTSIZE);
    add_scheduler_operator_button->label("Add");
    add_scheduler_operator_button->callback(cb_add_scheduler_operator, this);
    scheduler_operator_browser = new Fl_Hold_Browser(XJOBCOL2, GUIHEIGHT_EXT_START + height_var + 65, JOBCOLWIDTH, 161);
    scheduler_operator_browser->callback(cb_select_scheduler_operator, this);
    scheduler_operator_browser->textsize(RLN_FONTSIZE - 2);
    scheduler_operator_browser->end();
    int height_ops = height_var + 134;

    // Scheduler jobs
    Fl_Text_Buffer *textbuff3s = new Fl_Text_Buffer();
    textbuff3s->text("Jobs");
    Fl_Text_Display* textdisp3s = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT - 160, JOBCOLWIDTH - 50, 24);
    textdisp3s->buffer(textbuff3s);
    textdisp3s->textsize(12);
    textdisp3s->color(GUI_BACKGROUND_COLOR);
    scheduler_delete_job_button = new Fl_Button(XJOBCOL1 + JOBCOLWIDTH - 50, GUIHEIGHT_EXT - 160, 50, 23);
    scheduler_delete_job_button->color(GUI_BUTTON_COLOR);
    scheduler_delete_job_button->labelfont(FL_ITALIC);
    scheduler_delete_job_button->labelsize(RLN_FONTSIZE);
    scheduler_delete_job_button->label("Del");
    scheduler_delete_job_button->callback(cb_delete_scheduler_job, this);
    scheduler_job_browser = new Fl_Hold_Browser(XJOBCOL1, GUIHEIGHT_EXT - 160 + 23, JOBCOLWIDTH, 128);
    scheduler_job_browser->callback(cb_select_scheduled_job, this);
    scheduler_job_browser->textsize(RLN_FONTSIZE - 1);

    Fl_Text_Buffer *textbuff4s = new Fl_Text_Buffer();
    textbuff4s->text("Input to this job");
    Fl_Text_Display* textdisp4s = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT - 160, JOBCOLWIDTH, 24);
    textdisp4s->buffer(textbuff4s);
    textdisp4s->textsize(12);
    textdisp4s->color(GUI_BACKGROUND_COLOR);
    scheduler_input_job_browser = new Fl_Hold_Browser(XJOBCOL2, GUIHEIGHT_EXT - 160 + 24, JOBCOLWIDTH, 50);
    scheduler_input_job_browser->callback(cb_select_input_job, this);
    scheduler_input_job_browser->textsize(RLN_FONTSIZE - 1);

    Fl_Text_Buffer *textbuff5s = new Fl_Text_Buffer();
    textbuff5s->text("Output from this job");
    Fl_Text_Display* textdisp5s = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT - 160 + 76, JOBCOLWIDTH, 24);
    textdisp5s->buffer(textbuff5s);
    textdisp5s->textsize(12);
    textdisp5s->color(GUI_BACKGROUND_COLOR);
    scheduler_output_job_browser = new Fl_Hold_Browser(XJOBCOL2, GUIHEIGHT_EXT - 160 + 100, JOBCOLWIDTH, 50);
    scheduler_output_job_browser->callback(cb_select_output_job, this);
    scheduler_output_job_browser->textsize(RLN_FONTSIZE - 1);

    // Scheduler edges
    Fl_Text_Buffer *textbuffedge = new Fl_Text_Buffer();
    textbuffedge->text("Edges");
    Fl_Text_Display* textdispedge = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START + height_var, JOBCOLWIDTH - 105, 24);
    textdispedge->buffer(textbuffedge);
    textdispedge->textsize(12);
    textdispedge->color(GUI_BACKGROUND_COLOR);
    scheduler_edge_input= new Fl_Choice(XJOBCOL3, GUIHEIGHT_EXT_START + height_var + 23, JOBCOLWIDTH / 2 + 10, 21);
    scheduler_edge_input->color(GUI_INPUT_COLOR);
    scheduler_edge_input->textsize(RLN_FONTSIZE - 2);
    scheduler_edge_output = new Fl_Choice(XJOBCOL3 + 34 + JOBCOLWIDTH / 2, GUIHEIGHT_EXT_START + height_var + 23, JOBCOLWIDTH / 2 - 34, 21);
    scheduler_edge_output->label("->");
    scheduler_edge_output->color(GUI_INPUT_COLOR);
    scheduler_edge_output->textsize(RLN_FONTSIZE - 2);
    scheduler_edge_boolean = new Fl_Choice(XJOBCOL3 + 20, GUIHEIGHT_EXT_START + height_var + 44, JOBCOLWIDTH / 2 - 20, 21);
    scheduler_edge_boolean->label("if:");
    scheduler_edge_boolean->color(GUI_INPUT_COLOR);
    scheduler_edge_boolean->textsize(RLN_FONTSIZE-2);
    scheduler_edge_outputtrue = new Fl_Choice(XJOBCOL3 + 34 + JOBCOLWIDTH / 2, GUIHEIGHT_EXT_START + height_var + 44, JOBCOLWIDTH / 2 - 34, 21);
    scheduler_edge_outputtrue->label(":");
    scheduler_edge_outputtrue->textsize(RLN_FONTSIZE - 2);
    scheduler_edge_outputtrue->color(GUI_INPUT_COLOR);
    delete_scheduler_edge_button = new Fl_Button(XJOBCOL3 + JOBCOLWIDTH - 105, GUIHEIGHT_EXT_START + height_var, 50, 23);
    delete_scheduler_edge_button->color(GUI_BUTTON_COLOR);
    delete_scheduler_edge_button->labelfont(FL_ITALIC);
    delete_scheduler_edge_button->labelsize(RLN_FONTSIZE);
    delete_scheduler_edge_button->label("Del");
    delete_scheduler_edge_button->callback(cb_delete_scheduler_edge, this);
    add_scheduler_edge_button = new Fl_Button(XJOBCOL3 + JOBCOLWIDTH - 50, GUIHEIGHT_EXT_START + height_var, 50, 23);
    add_scheduler_edge_button->color(GUI_BUTTON_COLOR);
    add_scheduler_edge_button->labelfont(FL_ITALIC);
    add_scheduler_edge_button->labelsize(RLN_FONTSIZE);
    add_scheduler_edge_button->label("Add");
    add_scheduler_edge_button->callback(cb_add_scheduler_edge, this);
    scheduler_edge_browser = new Fl_Hold_Browser(XJOBCOL3, GUIHEIGHT_EXT_START + height_var + 65, JOBCOLWIDTH, 320);
    scheduler_edge_browser->callback(cb_select_scheduler_edge, this);
    scheduler_edge_browser->textsize(RLN_FONTSIZE - 2);
    scheduler_edge_browser->end();
    }

    scheduler_run_grp->end();  // scheduler_run_grp->end() already called

    }

    if (show_scheduler) {
        pipeliner_grp->hide();
        scheduler_grp->show();
        fillSchedulerNodesAndVariables();
        if (schedule.isWriteLocked()) {
            scheduler_run_grp->deactivate();
        }
    } else {
        scheduler_grp->hide();
        pipeliner_grp->show();
    }

    // B) Scheduler part of the GUI
    expand_stdout_grp = new Fl_Group(0, 0, 4 * w, 4 * h);
    {
    GroupContext context (expand_stdout_grp);

    disp_expand_stdout = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 - 5, w - 20, 300);
    disp_expand_stderr = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 - 5 + 305, w - 20, 85);
    disp_expand_stdout->fn_file = "run.out";
    disp_expand_stderr->fn_file = "run.err";
    textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
    disp_expand_stdout->buffer(textbuff_stdout);
    disp_expand_stderr->buffer(textbuff_stderr);
    disp_expand_stderr->textcolor(FL_RED);
    disp_expand_stdout->textsize(RLN_FONTSIZE - 1);
    disp_expand_stderr->textsize(RLN_FONTSIZE - 1);
    disp_expand_stdout->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
    disp_expand_stderr->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
    disp_expand_stdout->scrollbar_width(0);
    disp_expand_stderr->scrollbar_width(0);
    }

    if (!show_expand_stdout) expand_stdout_grp->hide();

    // Fill the actual browsers
    fillRunningJobLists();
    fillStdOutAndErr();

    // Mechanism to update stdout and stderr continuously and also update the JobLists
    // Also exit the GUI if it has been idle for too long
    update_every_sec = _update_every_sec;
    exit_after_sec = (float)_exit_after_sec;
    if (update_every_sec > 0)
        Fl::add_timeout(update_every_sec, Gui_Timer_CB, (void*)this);

    cb_show_initial_screen_i();

    // Set and activate current selection from side-browser
    cb_select_browsegroup_i(true);  // make default active; true is used to show_initial_screen
    is_main_continue = false;      // default is a new run
}

static void Gui_Timer_CB(void *userdata) {
    GuiMainWindow *o = (GuiMainWindow*)userdata;

    time_t now = time(nullptr);

    double dif = difftime(now, time_last_change);
    // If the GUI has been idle for too long, then exit
    if (dif > o->exit_after_sec) {
        std::cout << " The relion GUI has been idle for more than " << o->exit_after_sec << " seconds. Exiting... " << std::endl;
        exit(0);
    }

    if (show_scheduler) {
        // Always refill the stdout and stderr windows for scheduler
        o->fillStdOutAndErr();
        FileName mychanged = schedule.name + SCHEDULE_HAS_CHANGED;
        if (exists(mychanged)) {
            schedule.read(DONT_LOCK);
               o->fillSchedulerNodesAndVariables();
               std::remove(mychanged.c_str());
        }
    } else {
        // Update the stdout and stderr windows if we're currently pointing at a running job
        if (current_job >= 0 && pipeline.processList[current_job].status == Process::RUNNING)
            o->fillStdOutAndErr();

        // Check for job completion if the pipeline has been changed

        if (exists(PIPELINE_HAS_CHANGED))
            o->updateJobLists();
    }

    // Refresh every so many seconds
    Fl::repeat_timeout(o->update_every_sec, Gui_Timer_CB, userdata);
}

static void delete_menubar(Fl_Menu_Bar *menubar) {
    delete menubar;
    menubar = nullptr;
}

void GuiMainWindow::clear() {
    delete_menubar(menubar);
    delete_menubar(menubar2);
}

std::string GuiMainWindow::getJobNameForDisplay(Process &job) {
    FileName fn_pre, fn_jobnr, fn_post;

    if (show_scheduler) {
        return ((FileName) job.name).afterFirstOf(schedule.name).beforeLastOf("/");
    } else if (!decomposePipelineFileName(job.name, fn_pre, fn_jobnr, fn_post)) {
        return job.name;
    } else {
        return fn_jobnr.afterFirstOf("b").beforeFirstOf("/") + ": " + (job.alias == "None" ? job.name : job.alias);
    }
}

std::string decorate_name_or_alias(std::string name_or_alias, int status) {
    switch (status) {
        case Process::FINISHED_ABORTED:
        return "@C1@-@." + name_or_alias;
        case Process::FINISHED_FAILURE:
        return "@C1@." + name_or_alias;
        // case Process::FINISHED_SUCCESS:
        default:
        return name_or_alias;
    }
}

// Update the content of the finished, running and scheduled job lists
void GuiMainWindow::fillRunningJobLists() {
    // Go back to the same positions in the vertical scroll bars of the job lists after updating...
    int mypos_running = running_job_browser->position();
    int mypos_scheduled = scheduled_job_browser->position();
    int mypos_finished = finished_job_browser->position();
    int myhpos_running = running_job_browser->hposition();
    int myhpos_scheduled = scheduled_job_browser->hposition();
    int myhpos_finished = finished_job_browser->hposition();

    // Clear whatever was in there
    finished_job_browser->clear();
    finished_processes.clear();
    running_job_browser->clear();
    running_processes.clear();
    scheduled_job_browser->clear();
    scheduled_processes.clear();

    // Fill the finished Jobs browsers
    if (do_order_alphabetically) {
        // Only re-order the finished jobs!
        std::vector<std::pair<std::string, long int> > enumerate_jobs;
        for (long int i = pipeline.processList.size() - 1; i >= 0; i--) {
            enumerate_jobs.push_back(std::make_pair((
                pipeline.processList[i].alias != "None" ?
                pipeline.processList[i].alias : pipeline.processList[i].name
            ), i));
        }
        // Sort the pairs
        // (the first element of each pair will be used for the comparison)
        std::sort(enumerate_jobs.begin(), enumerate_jobs.end());

        for (auto &job : enumerate_jobs) {
            long int i = job.second;
            const auto &process = pipeline.processList[i];
            if (
                process.status == Process::FINISHED_SUCCESS ||
                process.status == Process::FINISHED_FAILURE ||
                process.status == Process::FINISHED_ABORTED
            ) {
                finished_processes.push_back(i);
                finished_job_browser->add(decorate_name_or_alias(job.first, process.status).c_str());
            }
        }
    } else {
        // For finished jobs, search backwards, so that last jobs are at the top
        for (long int i = pipeline.processList.size() - 1; i >= 0; i--) {
            if (
                pipeline.processList[i].status == Process::FINISHED_SUCCESS ||
                pipeline.processList[i].status == Process::FINISHED_FAILURE ||
                pipeline.processList[i].status == Process::FINISHED_ABORTED
            ) {
                finished_processes.push_back(i);
                finished_job_browser->add(decorate_name_or_alias(getJobNameForDisplay(pipeline.processList[i]), pipeline.processList[i].status).c_str());
            }
        }
    }

    // For running and scheduled jobs, search forwards, so that last jobs are at the bottom
    for (long int i = 0; i < pipeline.processList.size(); i++) {
        std::string jobname = getJobNameForDisplay(pipeline.processList[i]);
        switch (pipeline.processList[i].status) {

            case Process::RUNNING:
            running_processes.push_back(i);
            running_job_browser->add(jobname.c_str());
            break;

            case Process::SCHEDULED:
            scheduled_processes.push_back(i);
            scheduled_job_browser->add(jobname.c_str());
            break;

        }
    }

    running_job_browser  -> position(mypos_running);
    scheduled_job_browser-> position(mypos_scheduled);
    finished_job_browser -> position(mypos_finished);
    running_job_browser  ->hposition(myhpos_running);
    scheduled_job_browser->hposition(myhpos_scheduled);
    finished_job_browser ->hposition(myhpos_finished);
}

void GuiMainWindow::fillToAndFromJobLists() {
    display_io_node->clear();
    input_job_browser->clear();
    output_job_browser->clear();
    scheduler_input_job_browser->clear();
    scheduler_output_job_browser->clear();
    io_nodes.clear();
    input_processes.clear();
    output_processes.clear();

    if (current_job >= 0) {
        // Where do the input nodes come from?
        for (long int node : pipeline.processList[current_job].inputNodeList) {

            // no display for movie rootname
            if (pipeline.nodeList[node].type != Node::MOVIES) {
                FileName fnt = pipeline.nodeList[node].name;
                if (exists(fnt)) {
                    fnt = "in: " + fnt.afterLastOf("/");
                    display_io_node->add(fnt.c_str());
                    io_nodes.push_back(node);
                }
            }

            long int proc = pipeline.nodeList[node].outputFromProcess;
            if (proc >= 0) {
                // Check if this process was already there
                if (std::find(
                    input_processes.begin(), input_processes.end(), proc
                ) == input_processes.end()) {
                    input_processes.push_back(proc);
                    std::string jobname = getJobNameForDisplay(pipeline.processList[proc]);
                    if (show_scheduler) {
                        Fl_Hold_Browser *browser = scheduler_input_job_browser;
                    } else {
                        Fl_Select_Browser *browser = input_job_browser;
                    }
                    browser->add(jobname.c_str());
                }
            }
        }
        // Where do the output nodes lead to?
        for (long int node : pipeline.processList[current_job].outputNodeList) {
            FileName fnt = pipeline.nodeList[node].name;
            if (exists(fnt)) {
                fnt = "out: " + fnt.afterLastOf("/");
                display_io_node->add(fnt.c_str());
                io_nodes.push_back(node);
            }

            for (long int proc : pipeline.nodeList[node].inputForProcessList) {
                // Check if this process was already there
                if (std::find(
                    output_processes.begin(), output_processes.end(), proc
                ) == output_processes.end()) {
                    output_processes.push_back(proc);
                    std::string jobname = getJobNameForDisplay(pipeline.processList[proc]);
                    if (show_scheduler) {
                        scheduler_output_job_browser->add(jobname.c_str());
                    } else {
                        output_job_browser->add(jobname.c_str());
                    }
                }
            }
        }
    }
}

// Update the content of the finished, running and scheduled job lists
void GuiMainWindow::fillSchedulerNodesAndVariables() {
    // Go back to the same positions in the vertical scroll bars of the job lists after updating...
    int mypos_scheduler_variable = scheduler_variable_browser->value();
    int mypos_scheduler_operator = scheduler_operator_browser->value();
    int mypos_scheduler_edge = scheduler_edge_browser->value();
    int mypos_scheduler_job = scheduler_job_browser->value();

    // Clear whatever was in there
    scheduler_variable_browser->clear();
    scheduler_operator_browser->clear();
    scheduler_operator_output->clear();
    scheduler_operator_input1->clear();
    scheduler_operator_input2->clear();
    scheduler_edge_browser->clear();
    scheduler_edge_input->clear();
    scheduler_edge_output->clear();
    scheduler_edge_boolean->clear();
    scheduler_edge_outputtrue->clear();
    scheduler_job_browser->clear();
    scheduled_processes.clear();
    scheduler_current_node->clear();
    operators_list.clear();

    // Fill jobs browser
    for (long int i = 0; i < pipeline.processList.size(); i++) {
        scheduler_job_browser->add((getJobNameForDisplay(pipeline.processList[i])).c_str());
        scheduled_processes.push_back(i);
    }
    // Also get input/output

    // Fill edges browser
    for (const auto &edge : schedule.edges) {
        std::string label = edge.is_fork ?
            edge.inputNode + " -> (" + edge.myBooleanVariable + ") ? " + edge.outputNodeTrue + " : " + edge.outputNode :
            edge.inputNode + " -> " + edge.outputNode;
        scheduler_edge_browser->add(label.c_str());
    }

    // Fill variables browser, and pull-down menus for operator input/output
    for (const auto &x : schedule.getCurrentFloatVariables()) {
        std::string label = x.first
            + " = " + floatToString(x.second.value)
            + " (" + floatToString(x.second.original_value) + ")";
        scheduler_variable_browser->add(label.c_str());
        scheduler_operator_output->add(x.first.c_str());
        scheduler_operator_input1->add(x.first.c_str());
        scheduler_operator_input2->add(x.first.c_str());
    }

    for (const auto &p : schedule.getCurrentBooleanVariables()) {
        std::string label = p.first
            + " = " + (p.second.value ? "True" : "False")
            + " (" + (p.second.original_value ? "True" : "False") + ")";
        scheduler_variable_browser->add(label.c_str());
        scheduler_operator_output->add(p.first.c_str());
        scheduler_operator_input1->add(p.first.c_str());
        scheduler_operator_input2->add(p.first.c_str());
        scheduler_edge_boolean->add(p.first.c_str());
    }

    for (const auto &str : schedule.getCurrentStringVariables()) {
        std::string label = str.first
            + " = " + str.second.value
            + " (" + str.second.original_value + ")";
        scheduler_variable_browser->add(label.c_str());
        scheduler_operator_output->add(str.first.c_str());
        scheduler_operator_input1->add(str.first.c_str());
        scheduler_operator_input2->add(str.first.c_str());
    }

    // Fill operator browser
    for (const auto &opr : schedule.getCurrentOperators()) {
        std::string label = opr.first;
        scheduler_operator_browser->add(label.c_str());
        operators_list.push_back(opr.first);
        scheduler_edge_input->add(opr.first.c_str());
        scheduler_edge_output->add(opr.first.c_str());
        scheduler_edge_outputtrue->add(opr.first.c_str());
        scheduler_current_node->add(opr.first.c_str());
    }

    // Also add jobnames to the input/output nodes of the edges
    for (auto &process : pipeline.processList) {
        if (process.status == Process::SCHEDULED) {
            std::string jobname = getJobNameForDisplay(process);
            scheduler_edge_input->add(jobname.c_str());
            scheduler_edge_output->add(jobname.c_str());
            scheduler_edge_outputtrue->add(jobname.c_str());
            scheduler_current_node->add(jobname.c_str());
        }
    }

    // Set the value of the current_node
    // Set the current_node
    if (schedule.current_node != "undefined") {
        scheduler_current_node->value(scheduler_current_node->find_item(schedule.current_node.c_str()));
    }

    scheduler_operator_output->add("");
    scheduler_operator_input1->add("");
    scheduler_operator_input2->add("");
    scheduler_edge_outputtrue->add("");
    scheduler_edge_boolean->add("");

    if (mypos_scheduler_variable >= 0) {
        scheduler_variable_browser->value(mypos_scheduler_variable);
        select_scheduler_variable();
    }
    if (mypos_scheduler_operator >= 0) {
        scheduler_operator_browser->value(mypos_scheduler_operator);
        cb_select_scheduler_operator_i();
    }
    if (mypos_scheduler_edge >= 0) {
        scheduler_edge_browser->value(mypos_scheduler_edge);
        cb_select_scheduler_edge_i();
    }
    if (mypos_scheduler_job > 0) {
        scheduler_job_browser->value(mypos_scheduler_job);
    }

    if (schedule.isWriteLocked()) {
        scheduler_run_grp->deactivate();
    } else {
        scheduler_run_grp->activate();
    }
}

void GuiMainWindow::fillStdOutAndErr() {
    FileName fn_out = "";
    FileName fn_err = "";
    FileName fn_outtail, fn_errtail;
    if (current_job >= 0 || show_scheduler) {
        std::string myroot = show_scheduler ? schedule.name : pipeline.processList[current_job].name;
        fn_out = myroot + "run.out";
        fn_err = myroot + "run.err";
        fn_outtail = myroot + ".run.out.tail";
        fn_errtail = myroot + ".run.err.tail";
    }

    if (exists(fn_out)) {
        if (maingui_do_read_only) {
            int err = textbuff_stdout->loadfile(fn_out.c_str());
        } else {
            // Remove annoying carriage returns
            std::string command = "tail -n 20 < " + fn_out + " | awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' > " + fn_outtail;
            int res = system(command.c_str());
            std::ifstream in(fn_outtail.c_str(), std::ios_base::in);
            if (in.fail())
                REPORT_ERROR((std::string) "MetaDataTable::read: File " + fn_outtail + " does not exist");
            int err = textbuff_stdout->loadfile(fn_outtail.c_str());
        }
        // Scroll to the bottom
        disp_stdout->insert_position(textbuff_stdout->length()-1);
        disp_stdout->show_insert_position();
        disp_expand_stdout->insert_position(textbuff_stdout->length()-1);
        disp_expand_stdout->show_insert_position();
    } else {
        textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    }

    if (exists(fn_err)) {
        if (maingui_do_read_only) {
            int err = textbuff_stderr->loadfile(fn_err.c_str());
        } else {
            std::string command = "tail -10 " + fn_err + " > " + fn_errtail;
            int res = system(command.c_str());
            std::ifstream in(fn_errtail.c_str(), std::ios_base::in);
            if (in.fail())
                REPORT_ERROR((std::string) "MetaDataTable::read: File " + fn_errtail + " does not exist");
            int err = textbuff_stderr->loadfile(fn_errtail.c_str());
        }
        // Scroll to the bottom
        disp_stderr->insert_position(textbuff_stderr->length()-1);
        disp_stderr->show_insert_position();
        disp_expand_stderr->insert_position(textbuff_stderr->length()-1);
        disp_expand_stderr->show_insert_position();
    } else {
        textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
    }
}

void GuiMainWindow::tickTimeLastChanged() {
    time(&time_last_change);
}

void GuiMainWindow::updateJobLists() {
    pipeline.checkProcessCompletion();
    if (show_scheduler)
        fillSchedulerNodesAndVariables();
    else
        fillRunningJobLists();
    fillToAndFromJobLists();
}

void GuiMainWindow::loadJobFromPipeline(int this_job) {
    // Set the "static int" to which job we're currently pointing
    current_job = this_job;
    int itype = pipeline.processList[current_job].type;
    // The following line allows certain browse buttons to only open the current directory (using CURRENT_ODIR)
    current_browse_directory = pipeline.processList[current_job].name;

    // What type of job is this?
    for (int t = 0; t < NR_BROWSE_TABS; t++) {
        if (gui_jobwindows[t]->myjob.type == itype)
            browser->value(t + 1);
    }

    // change GUI to the corresponding jobwindow
    cb_select_browsegroup_i();

    // Re-read the settings for this job and update the values inside the GUI
    int iwin = browser->value() - 1;
    gui_jobwindows[iwin]->myjob.read(pipeline.processList[current_job].name, is_main_continue);
    gui_jobwindows[iwin]->updateMyGui();

    // If a finished or running job was loaded from the pipeline: set this to be a continuation job
    // If a scheduled job was loaded, only set is_main_continue to true when it is Process::SCHEDULED
    //if (pipeline.processList[current_job].status == Process::SCHEDULED && !gui_jobwindows[iwin]->myjob.is_continue)
    //	is_main_continue = false;
    //else
    //	is_main_continue = true;

    // Any job loaded from the pipeline will initially be set as a continuation job
    // but for show-scheduler, no job should be a continuation
    if (show_scheduler) {
        is_main_continue = false;
        do_overwrite_continue = true;
    } else {
        is_main_continue = true;
    }
    cb_toggle_continue_i();

    // Set the alias in the window
    alias_current_job->value((getJobNameForDisplay(pipeline.processList[current_job])).c_str());
    alias_current_job->position(0);  //left-centered text in box

    // Update all job lists in the main GUI
    updateJobLists();

    // File the out and err windows
    fillStdOutAndErr();
}

void GuiMainWindow::cb_select_browsegroup(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    // When clicking the job browser on the left: reset current_job to -1 (i.e. a new job, not yet in the pipeline)
    current_job = -1;
    T->cb_select_browsegroup_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_browsegroup_i(bool show_initial_screen) {
    // Update timer
    tickTimeLastChanged();

    // Hide the initial screen
    if (show_initial_screen)
        background_grp->show();
    else
        background_grp->hide();

    int iwin = browser->value() - 1;
    if (iwin < 0 || iwin >= NR_BROWSE_TABS) return;
    // Show the 'selected' group. Hide the others.
    for (int t = 0; t < NR_BROWSE_TABS; t++) {
        // During the initial screen: show a nice picture with some explanations
        // browser starts counting at 1...
        if (t == iwin && !show_initial_screen) {
            browse_grp[t]->show();
        } else {
            browse_grp[t]->hide();
        }
    }

    // Update all job lists in the main GUI
    updateJobLists();

    is_main_continue = false;
    do_overwrite_continue = false;

    // If the GUI got changed, put that change into the joboption now
    gui_jobwindows[iwin]->updateMyJob();

    // toggle the continue status of this job
    cb_toggle_continue_i();

    alias_current_job->value("Give_alias_here");

    scheduler_job_name->value("");
    scheduler_job_name->activate();
    scheduler_job_has_started->deactivate();
    scheduler_job_has_started->picked(&job_has_started_options[1]);  // initialise to has_not_started

    textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
}

void GuiMainWindow::cb_select_finished_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_finished_job_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_finished_job_i() {
    // Update timer
    tickTimeLastChanged();

    // Show the 'selected' group. Hide the others.
    int idx = finished_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
        loadJobFromPipeline(finished_processes[idx]);
}

void GuiMainWindow::cb_select_running_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_running_job_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_running_job_i() {
    // Update timer
    tickTimeLastChanged();

    // Show the 'selected' group. Hide the others.
    int idx = running_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
        loadJobFromPipeline(running_processes[idx]);

}

void GuiMainWindow::cb_select_scheduled_job(Fl_Widget* o, void* v) {
//	std::cout << "v = " << v << std::endl;
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_scheduled_job_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_scheduled_job_i() {
    // Update timer
    tickTimeLastChanged();

    // Show the 'selected' group. Hide the others.
    int idx = (show_scheduler ? scheduler_job_browser->value() :
                                scheduled_job_browser->value() ) - 1;
    if (idx < 0) return;

    // only if a non-empty line was selected
    loadJobFromPipeline(scheduled_processes[idx]);

    if (show_scheduler) {
        FileName jobname = getJobNameForDisplay(pipeline.processList[current_job]);
        scheduler_job_name->value(jobname.c_str());
        scheduler_job_name->deactivate();
        bool found = false;
        for (int i = 0; i < 3; i++) {
            if (schedule.jobs[jobname].mode == job_mode_options[i].label()) {
                found = true;
                scheduler_job_mode->value(i);
            }
        }
        scheduler_job_has_started->value(schedule.jobs[jobname].job_has_started);
        scheduler_job_has_started->activate();

        if (!found) REPORT_ERROR("ERROR: unrecognised job_mode ...");
    }
}

void GuiMainWindow::cb_select_input_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_input_job_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_input_job_i() {
    // Update timer
    tickTimeLastChanged();

    // Show the 'selected' group. Hide the others.
    if (show_scheduler) {
        Fl_Hold_Browser *browser = scheduler_input_job_browser;
    } else {
        Fl_Select_Browser *browser = input_job_browser;
    }
    // only if a non-empty line was selected
    int i = browser->value() - 1;
    if (i >= 0)
        loadJobFromPipeline(input_processes[i]);
}

void GuiMainWindow::cb_select_output_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_output_job_i();
    run_button->activate();
}

void GuiMainWindow::cb_select_output_job_i() {
    // Update timer
    tickTimeLastChanged();

    // Show the 'selected' group. Hide the others.
    if (show_scheduler) {
        Fl_Hold_Browser *browser = scheduler_output_job_browser;
    } else {
        Fl_Select_Browser *browser = output_job_browser;
    }

    // only if a non-empty line was selected
    int i = browser->value() - 1;
    if (i >= 0)
        loadJobFromPipeline(output_processes[i]);
}

void GuiMainWindow::cb_display_io_node(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_display_io_node_i();
    run_button->activate();
}

void GuiMainWindow::cb_display_io_node_i() {
    // Run relion_display on the output node
    long int mynode = io_nodes[display_io_node->value()];
    std::string command;

    if (pipeline.nodeList[mynode].type == Node::MIC_COORDS) {
        // A manualpicker jobwindow for display of micrographs....
        RelionJob manualpickjob;
        FileName fn_job = ".gui_manualpick";
        bool iscont = false;
        if (exists(fn_job + "job.star") || exists(fn_job + "run.job")) {
            manualpickjob.read(fn_job.c_str(), iscont, true);  // true means do initialise
        } else {
            fl_message("ERROR: Save a Manual picking job parameter file (using the Save jobs settings option from the Jobs menu) before displaying coordinate files. ");
            return;
        }

        // Get the name of the micrograph STAR file from reading the suffix file
        FileName fn_suffix = pipeline.nodeList[mynode].name;
        if (fn_suffix.getExtension() == "star") {
            std::ifstream in (fn_suffix.data(), std::ios_base::in);
            FileName fn_star;
            in >> fn_star ;
            in.close();
            if (!fn_star.empty()) {
                FileName fn_dirs = fn_suffix.beforeLastOf("/") + "/";
                fn_suffix = fn_suffix.afterLastOf("/").without("coords_suffix_");
                fn_suffix = fn_suffix.withoutExtension();
                // Launch the manualpicker...
                command = "`which relion_manualpick` --i " + fn_star
                        + " --odir " + fn_dirs
                        + " --pickname " + fn_suffix
                        + " --scale " + manualpickjob.joboptions["micscale"].getString()
                        + " --sigma_contrast " + manualpickjob.joboptions["sigma_contrast"].getString()
                        + " --black " + manualpickjob.joboptions["black_val"].getString()
                        + " --white " + manualpickjob.joboptions["white_val"].getString();

                try {
                    float mylowpass = manualpickjob.joboptions["lowpass"].getNumber();
                    if (mylowpass > 0.0)
                    command += " --lowpass " + manualpickjob.joboptions["lowpass"].getString();
                } catch (std::string errmsg) {
                    fl_message("joboption[\"lowpass\"] %s", errmsg.c_str());
                    return;
                }

                try {
                    float myhighpass = manualpickjob.joboptions["highpass"].getNumber();
                    if (myhighpass > 0.0)
                    command += " --highpass " + manualpickjob.joboptions["highpass"].getString();
                } catch (std::string errmsg) {
                    fl_message("joboption[\"highpass\"] %s", errmsg.c_str());
                    return;
                }

                try {
                    float myangpix = manualpickjob.joboptions["angpix"].getNumber();
                    if (myangpix > 0.0)
                    command += " --angpix " + manualpickjob.joboptions["angpix"].getString();
                } catch (std::string errmsg) {
                    fl_message("joboption[\"angpix\"] %s", errmsg.c_str());
                    return;
                }

                command += " --ctf_scale " + manualpickjob.joboptions["ctfscale"].getString();

                command += " --particle_diameter " + manualpickjob.joboptions["diameter"].getString();

                if (manualpickjob.joboptions["do_color"].getBoolean()) {
                    command += " --color_label " + manualpickjob.joboptions["color_label"].getString();
                    command += " --blue " + manualpickjob.joboptions["blue_value"].getString();
                    command += " --red " + manualpickjob.joboptions["red_value"].getString();
                    if (!manualpickjob.joboptions["fn_color"].getString().empty())
                    command += " --color_star " + manualpickjob.joboptions["fn_color"].getString();
                }

                // Other arguments for extraction
                command += " " + manualpickjob.joboptions["other_args"].getString() + " &";
            } else {
                fl_message("Only coordinates in .star format, generated in the pipeline, can be displayed here.");
            }
        } else {
            fl_message("Only coordinates in .star format, generated in the pipeline, can be displayed here.");
        }
    } else if (pipeline.nodeList[mynode].type == Node::PDF_LOGFILE) {
        const char *default_pdf_viewer = getenv ("RELION_PDFVIEWER_EXECUTABLE");
        if (!default_pdf_viewer) { default_pdf_viewer = DEFAULT::PDFVIEWER; }
        std::string myviewer(default_pdf_viewer);
        command = myviewer + " " + pipeline.nodeList[mynode].name + "&";
    } else if (pipeline.nodeList[mynode].type == Node::POLISH_PARAMS) {
        command = "cat " + pipeline.nodeList[mynode].name;
    } else if (pipeline.nodeList[mynode].type != Node::POST) {
        command = "relion_display --gui --i " + pipeline.nodeList[mynode].name + " &";
    }
    // std::cerr << " command= " << command << std::endl;
    int res = system(command.c_str());
}

std::string scheduler_get_io(Fl_Choice *scheduler_edge_io, std::string io) {
    int i = scheduler_edge_io->value();
    if (i < 0 || i >= scheduler_edge_io->size()) {
        std::cerr << " Error getting " + io + " from scheduler edge window, please try again ..." << std::endl;
        return "";  // The error handling here depends on the empty string not being a normal return value.
    } else {
        return scheduler_edge_io->text(i);
    }
}


void GuiMainWindow::cb_add_scheduler_edge(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    // Get input
    std::string input = scheduler_get_io(scheduler_edge_input, "input");
    if (input.empty()) return;

    // Get output
    std::string output = scheduler_get_io(scheduler_edge_output, "output");
    if (output.empty()) return;

    /// TODO: Move schedule.read(DO_LOCK) and schedule.write(DO_LOCK) out of the if-else block.
    int i = scheduler_edge_boolean->value();
    if (i >= 0) {
        // Get outputtrue
        std::string mybool = scheduler_edge_boolean->text(i);
        int j = scheduler_edge_outputtrue->value();
        if (j < 0 || j >= scheduler_edge_outputtrue->size()) {
            std::cerr << " Error getting outputtrue from scheduler edge window, please try again ..." << std::endl;
            return;
        } else {
            std::string outputtrue = scheduler_edge_outputtrue->text(j);
            Schedule::rwlock lock (schedule);
            schedule.addFork(input, mybool, outputtrue, output);
        }
    } else {
        Schedule::rwlock lock (schedule);
        schedule.addEdge(input, output);
    }

    T->fillSchedulerNodesAndVariables();
}

static bool user_wants_to(std::string action, std::string describe_action) {
    return fl_choice("%s", "Cancel", action.c_str(), nullptr, ("Are you sure you want to " + describe_action + "?").c_str());
}

void GuiMainWindow::cb_delete_scheduler_edge(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    int i = scheduler_edge_browser->value();
    if (i <= 0) {
        fl_message("Please select a job.");
        return;
    }

    if (!user_wants_to("Delete!", "delete this edge")) {
        do_overwrite_continue = false;
        return;
    }

    {
        Schedule::rwlock lock (schedule);
        schedule.removeEdge(i - 1);
    }
    // Also reset entry fields
    scheduler_edge_input->value(-1);
    scheduler_edge_output->value(-1);
    scheduler_edge_outputtrue->value(-1);
    scheduler_edge_boolean->value(-1);
    scheduler_edge_browser->value(-1);
    T->fillSchedulerNodesAndVariables();
}

void GuiMainWindow::cb_select_scheduler_edge(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_scheduler_edge_i();
}

void GuiMainWindow::cb_select_scheduler_edge_i() {
    // Get position of the browser:
    int i = scheduler_edge_browser->value();
    if (i >= 1) {

        FileName mytext = scheduler_edge_browser->text(i);

        scheduler_edge_input->value(scheduler_edge_input->find_item(schedule.edges[i - 1].inputNode.c_str()));
        scheduler_edge_output->value(scheduler_edge_output->find_item(schedule.edges[i - 1].outputNode.c_str()));
        if (schedule.edges[i - 1].is_fork) {
            scheduler_edge_boolean->value(scheduler_edge_boolean->find_item(schedule.edges[i - 1].myBooleanVariable.c_str()));
            scheduler_edge_outputtrue->value(scheduler_edge_outputtrue->find_item(schedule.edges[i - 1].outputNodeTrue.c_str()));
        } else {
            scheduler_edge_boolean->value(scheduler_edge_boolean->find_item(""));
            scheduler_edge_outputtrue->value(scheduler_edge_outputtrue->find_item(""));
        }

    } else {
        scheduler_edge_input->value(-1);
        scheduler_edge_output->value(-1);
        scheduler_edge_boolean->value(-1);
        scheduler_edge_outputtrue->value(-1);
    }
}

void GuiMainWindow::cb_set_scheduler_variable(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    std::string sched_var_name = scheduler_variable_name->value();
    std::string sched_var_val = scheduler_variable_value->value();

    if (sched_var_name.empty()) return;

    Schedule::rwlock lock (schedule);
    schedule.setVariable(sched_var_name, sched_var_val);
    schedule.setOriginalVariable(sched_var_name, sched_var_val);
    // Also reset entry fields
    scheduler_variable_name->value("");
    scheduler_variable_value->value("");
    T->fillSchedulerNodesAndVariables();
}

void GuiMainWindow::cb_delete_scheduler_variable(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    std::string myname = scheduler_variable_name->value();
    if (myname.empty()) return;

    if (!user_wants_to("Delete!", "delete this variable, and all operators or edges that use it")) {
        do_overwrite_continue = false;
        return;
    }

    {
        Schedule::rwlock lock (schedule);
        schedule.removeVariable(myname);
    }
    // Also reset entry fields
    scheduler_variable_name->value("");
    scheduler_variable_value->value("");
    T->fillSchedulerNodesAndVariables();
}

void GuiMainWindow::cb_add_scheduler_operator(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    int i = scheduler_operator_type->value();
    if (i < 0) {
        std::cerr << "ERROR: select an operator type. Try again... " << std::endl;
        return;
    }
    std::string type = scheduler_operator_type->text(i);
    i = scheduler_operator_output->value();
    std::string output = i < 0 || i >= scheduler_operator_output->size() ? "" : scheduler_operator_output->text(i);
    i = scheduler_operator_input1->value();
    std::string input1 = i < 0 || i >= scheduler_operator_input1->size() ? "" : scheduler_operator_input1->text(i);
    i = scheduler_operator_input2->value();
    std::string input2 = i < 0 || i >= scheduler_operator_input2->size() ? "" : scheduler_operator_input2->text(i);
    SchedulerOperator op;
    try {
        op = schedule.initialiseOperator(type, input1, input2, output);
    } catch (std::string errmsg) {
        fl_message("%s", errmsg.c_str());
        return;
    }

    if (schedule.isOperator(op.getName())) {
        fl_message("ERROR: this operator already exists...");
        return;
    }

    {
        Schedule::rwlock lock (schedule);
        schedule.addOperator(op);
    }
    // Also reset entry fields
    scheduler_operator_type->value(-1);
    scheduler_operator_output->value(-1);
    scheduler_operator_input1->value(-1);
    scheduler_operator_input2->value(-1);
    T->fillSchedulerNodesAndVariables();
}

void GuiMainWindow::cb_delete_scheduler_operator(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (!user_wants_to("Delete!", "delete this operator and its connecting edges")) {
        do_overwrite_continue = false;
        return;
    }

    const std::string type = scheduler_operator_type->text(scheduler_operator_type->value());
    std::string output = "", input1 = "", input2 = "";

    // Some operators do not have these arguments.
    if (scheduler_operator_output->value() >= 0)
        output = scheduler_operator_output->text(scheduler_operator_output->value());
    if (scheduler_operator_input1->value() >= 0)
        input1 = scheduler_operator_input1->text(scheduler_operator_input1->value());
    if (scheduler_operator_input2->value() >= 0)
        input2 = scheduler_operator_input2->text(scheduler_operator_input2->value());

    {
        Schedule::rwlock lock (schedule);
        schedule.removeOperator(schedule.getOperatorName(type, input1, input2, output));
    }

    // Also reset entry fields
    scheduler_operator_type->value(-1);
    scheduler_operator_output->value(-1);
    scheduler_operator_input1->value(-1);
    scheduler_operator_input2->value(-1);

    T->fillSchedulerNodesAndVariables();
}

void GuiMainWindow::cb_select_scheduler_variable(Fl_Widget* o, void*v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->select_scheduler_variable();
}

void GuiMainWindow::select_scheduler_variable() {
    // Get position of the browser:
    int i = scheduler_variable_browser->value();
    if (i >= 1) {
        FileName mytext = scheduler_variable_browser->text(i);
        FileName myname = mytext.beforeFirstOf(" = ");
        FileName myval = mytext.afterFirstOf(" = ");
        myval = myval.beforeFirstOf(" (");
        scheduler_variable_name->value(myname.c_str());
        scheduler_variable_value->value(myval.c_str());
    } else {
        scheduler_variable_name->value("");
        scheduler_variable_value->value("");
    }
}

void GuiMainWindow::cb_select_scheduler_operator(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_select_scheduler_operator_i();
}

void GuiMainWindow::cb_select_scheduler_operator_i() {
    // Get position of the browser:
    int i = scheduler_operator_browser->value();
    if (i >= 1) {
        FileName myname = scheduler_operator_browser->text(i);
        std::string type, input1, input2, output;
        schedule.getOperatorParameters(myname, type, input1, input2, output);

        scheduler_operator_type->value(scheduler_operator_type->find_item(type.c_str()));
        scheduler_operator_output->value(scheduler_operator_output->find_item(
            scheduler_operator_output->find_item(output.c_str()) ? output.c_str() : ""
        ));
        scheduler_operator_input1->value(scheduler_operator_input1->find_item(
            scheduler_operator_input1->find_item(input1.c_str()) ? input1.c_str() : ""
        ));
        scheduler_operator_input2->value(scheduler_operator_input2->find_item(
            scheduler_operator_input2->find_item(input2.c_str()) ? input2.c_str() : ""
        ));
    } else {
        scheduler_operator_type->value(-1);
        scheduler_operator_output->value(-1);
        scheduler_operator_input1->value(-1);
        scheduler_operator_input2->value(-1);
    }
}

void GuiMainWindow::cb_scheduler_set_current(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->scheduler_set_current();
}

void GuiMainWindow::scheduler_set_current() {
    if (scheduler_current_node->value() < 0) {
        std::cerr << " ERROR: scheduler_current_node->value()= " << scheduler_current_node->value() << std::endl;
        return;
    }

    {
        Schedule::rwlock lock (schedule);
        schedule.current_node = std::string(
            scheduler_current_node->text(scheduler_current_node->value())
        );
    }

    // If a schedule has finished: activate the GUI again
    if (schedule.current_node == "EXIT") {
        scheduler_run_grp->activate();
    }
    if (schedule.isJob(schedule.current_node)) {
        for (long int id : scheduled_processes) {
            if (schedule.current_node == getJobNameForDisplay(pipeline.processList[id])) {
                scheduler_job_browser->value(id + 1);
                cb_select_scheduled_job_i();
            }
        }
    } else {
        for (int i = 0; i < operators_list.size(); i++) {
            if (schedule.current_node == operators_list[i]) {
                scheduler_operator_browser->value(i + 1);
                cb_select_scheduler_operator_i();
            }
        }
    }

    // Also set the edge from this node to the next one!
    for (int i = 0; i < schedule.edges.size(); i++) {
        if (schedule.edges[i].inputNode == schedule.current_node) {
            scheduler_edge_browser->value(i + 1);
            cb_select_scheduler_edge_i();
        }
    }

}

void GuiMainWindow::cb_scheduler_next(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    std::string mycurrent = schedule.current_node;
    if (schedule.current_node == "undefined") {
        if (schedule.edges.size() > 0) {
            schedule.current_node = schedule.edges[0].inputNode;
            scheduler_current_node->value(scheduler_current_node->find_item(schedule.current_node.c_str()));
            T->scheduler_set_current();
        }
        return;
    }
    for (const auto &edge : schedule.edges) {
        if (edge.inputNode == mycurrent) {
            std::string nextnode;
            if (edge.is_fork) {
                std::string question = "Fork on " + edge.myBooleanVariable + ". Do you want this to be True or False?";
                nextnode = fl_choice(
                    "%s", "False", "True", nullptr, question.c_str()
                ) ? edge.outputNodeTrue :
                    edge.outputNode;
            } else {
                nextnode = edge.outputNode;
            }
            const Fl_Menu_Item *myitem = scheduler_current_node->find_item(nextnode.c_str());
            if (!myitem) {
                fl_message("ERROR: next node is undefined");
            } else {
                scheduler_current_node->value(myitem);
                T->scheduler_set_current();
            }
            return;
        }
    }
}

void GuiMainWindow::cb_scheduler_prev(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    // If already at the beginning, just return
    if (schedule.current_node == "undefined") return;

    std::string myprev = schedule.getPreviousNode();
    std::cerr << " myprev= " << myprev << std::endl;
    if (myprev == "undefined") {
        fl_message("ERROR: previous node is undefined");
    } else {
        scheduler_current_node->value(scheduler_current_node->find_item(myprev.c_str()));
        T->scheduler_set_current();
    }
}

void GuiMainWindow::cb_scheduler_reset(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    if (user_wants_to("Reset!", "reset all variables to their initial state, in order to start over from scratch")) {
        {
            Schedule::rwlock lock (schedule);
            schedule.reset();
        }
        T->fillSchedulerNodesAndVariables();
        T->scheduler_set_current();
    }
}

void GuiMainWindow::cb_scheduler_unlock(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    schedule.unlock();
    show_expand_stdout = true;
    T->toggle_expand_stdout();
    scheduler_run_grp->activate();
}

void GuiMainWindow::cb_scheduler_abort(Fl_Widget *o, void* v) {
    if (user_wants_to("Abort!", "abort this schedule"))
        schedule.abort();
}

void GuiMainWindow::cb_scheduler_run(Fl_Widget *o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    FileName name_wo_dir = schedule.name;

    std::string command = " relion_scheduler --schedule " + name_wo_dir.afterFirstOf("Schedules/") + " --run --pipeline_control " + schedule.name + " >> "
        + schedule.name + "run.out 2>> " + schedule.name + "run.err &";
    int res = system(command.c_str());
    scheduler_run_grp->deactivate();

    show_expand_stdout = false;
    T->toggle_expand_stdout();
}

void GuiMainWindow::cb_display(Fl_Widget* o, void* v) {
    system("relion_display --gui &");
}

void GuiMainWindow::cb_toggle_continue_i() {

    if (is_main_continue || do_overwrite_continue) {
        if (do_overwrite_continue) {
            run_button->label("Overwrite!");
            add_job_button->label("Save");
            add_job_button->color(GUI_BUTTON_COLOR);
        } else {
            run_button->label("Continue!");
        }
        run_button->color(GUI_BUTTON_COLOR);
        run_button->labelfont(FL_ITALIC);
        run_button->labelsize(13);
        alias_current_job->deactivate();
    } else {
        run_button->label("Run!");
        add_job_button->label("Add job");
        add_job_button->color(GUI_RUNBUTTON_COLOR);
        run_button->color(GUI_RUNBUTTON_COLOR);
        run_button->labelfont(FL_ITALIC);
        run_button->labelsize(16);
        alias_current_job->activate();
    }

    int my_window = browser->value() - 1;
    gui_jobwindows[my_window]->toggle_new_continue(is_main_continue && !do_overwrite_continue);
}

void GuiMainWindow::cb_print_cl(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    int iwin = browser->value() - 1;
    // And update the job inside it
    gui_jobwindows[iwin]->updateMyJob();

    try {
        pipeline.getCommandLineJob(
            gui_jobwindows[iwin]->myjob, current_job, is_main_continue, false,
            DONT_MKDIR, do_overwrite_continue, T->commands, T->final_command
        );
    } catch (std::string errmsg) {
        fl_message("%s", errmsg.c_str());
        return;
    }

    fl_input("%s", join(T->commands, " && ").c_str(), " The command is: ");
    // Don't free the returned string! It comes from Fl_Input::value(), which returns
    // "pointer to an internal buffer - do not free() this".
}

// Run button callback functions
void GuiMainWindow::cb_run(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    // Deactivate Run button to prevent the user from accidentally submitting many jobs
    run_button->deactivate();
    // Run the job
    T->run(false, false);  // 1st false means dont only_schedule, 2nd false means dont open the note editor window
}

// Run button callback functions
void GuiMainWindow::cb_schedule(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->run(true, false);  // 1st true means only_schedule, do not run, 2nd false means dont open the note editor window
}

void GuiMainWindow::run(bool only_schedule, bool do_open_edit) {
    if (do_overwrite_continue && user_wants_to("Overwrite!", "overwrite this job")) {
        do_overwrite_continue = false;
        return;
    }

    // Get which jobtype the GUI is on now
    int iwin = browser->value() - 1;
    // And update the job inside it
    gui_jobwindows[iwin]->updateMyJob();

    // Update timer
    tickTimeLastChanged();

    try {
        pipeline.runJob(
            gui_jobwindows[iwin]->myjob, current_job, only_schedule,
            is_main_continue, false, do_overwrite_continue
        );
    } catch (const std::string &errmsg) {
        fl_message("%s", errmsg.c_str());
        // Allow the user to fix the error and submit this job again
        run_button->activate();
        return;
    }

    // Update all job lists in the main GUI
    updateJobLists();

    // Open the edit note window
    if (do_open_edit) {
        // Open the note editor window
        cb_edit_note_i();
    }

    // Also set alias from the alias_current_job input
    if (!is_main_continue && !do_overwrite_continue) {
        std::string alias = alias_current_job->value();
        if (alias != "Give_alias_here" && alias != pipeline.processList[current_job].name)
            set_alias(alias);
    }

    do_overwrite_continue = false;

    // Select this job now
    loadJobFromPipeline(current_job);
}

void GuiMainWindow::cb_delete_scheduler_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    std::vector<bool> deleteProcesses, deleteNodes;
    pipeline.deleteJobGetNodesAndProcesses(current_job, true, deleteNodes, deleteProcesses);

    // Before we do anything: confirm this is really what the user wants to do....
    std::string describe_action = "delete the following jobs, and their connecting edges? \n";
    for (size_t i = 0; i < deleteProcesses.size(); i++) {
        if (deleteProcesses[i])
            describe_action += " - " + T->getJobNameForDisplay(pipeline.processList[i]) + "\n";
    }
    if (user_wants_to("Move", describe_action)) {

        // Remove the jobs from the schedule itself
        {
        Schedule::rwlock lock (schedule);
        for (int i = 0; i < deleteProcesses.size(); i++)
            if (deleteProcesses[i])
                schedule.removeJob(T->getJobNameForDisplay(pipeline.processList[i]));
        }

        // And remove from the local pipeliner
        pipeline.deleteNodesAndProcesses(deleteNodes, deleteProcesses);

        // Reset current_job
        current_job = -1;
        scheduler_job_name->value("");
        T->fillStdOutAndErr();

        // Update all job lists in the main GUI
        T->updateJobLists();
    }
}

void GuiMainWindow::cb_scheduler_add_job(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    // Get which jobtype the GUI is on now
    int iwin = browser->value() - 1;
    // And update the job inside it
    gui_jobwindows[iwin]->updateMyJob();

    std::string mode = job_mode_options[scheduler_job_mode->value()].label();
    std::string jobname = scheduler_job_name->value();

    if (do_overwrite_continue) {
        // Write job settings, which might have been updated.
        gui_jobwindows[iwin]->myjob.write(pipeline.processList[current_job].name);

        Schedule::rwlock lock (schedule);
        // Write job_mode, which might have been updated.
        schedule.jobs[jobname].mode = job_mode_options[scheduler_job_mode->value()].label();
        schedule.jobs[jobname].job_has_started = (job_has_started_options[scheduler_job_has_started->value()].label() == "has started");
    } else {
        // Add job to the schedule
        // Get the mode, and the jobname
        if (jobname.empty()) {
            fl_message("%s","You need to provide a Name for this job in the scheduler.");
            return;
        }

        /// TODO: test the command line
        try {
            std::string dummy;
            T->final_command = gui_jobwindows[iwin]->myjob.getCommands(dummy, T->commands, false, 1);
        } catch (const std::string &errmsg) {
            fl_message("%s", errmsg.c_str());
            return;
        }

        {
        Schedule::rwlock lock (schedule);
        schedule.addJob(gui_jobwindows[iwin]->myjob, jobname, mode);
        }

        scheduler_job_name->value("");
        T->updateJobLists();
    }
}

// Run button callback functions
void GuiMainWindow::cb_delete(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (current_job < 0) {
        fl_message("Please select a job.");
        return;
    }

    bool ask = true, recursively = true;
    std::vector<bool> deleteProcesses, deleteNodes;
    pipeline.deleteJobGetNodesAndProcesses(current_job, recursively, deleteNodes, deleteProcesses);

    // Before we do anything: confirm this is really what the user wants to do....
    std::string describe_action = "";
    if (ask) {
        describe_action = "move the following processes to Trash? \n";
        for (size_t i = 0; i < deleteProcesses.size(); i++) {
            if (deleteProcesses[i]) {
                Process job = pipeline.processList[i];
                describe_action += " - " + (job.alias == "None" ? job.name : job.alias) + "\n";
            }
        };
    }

    if (!ask || user_wants_to("Move", describe_action)) {

        pipeline.deleteNodesAndProcesses(deleteNodes, deleteProcesses);

        // Reset current_job
        current_job = -1;
        T->fillStdOutAndErr();

        // Update all job lists in the main GUI
        T->updateJobLists();
    }
}

// Run button callback functions
void GuiMainWindow::cb_gently_clean_all_jobs(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->clean_all_jobs(false);
}

// Run button callback functions
void GuiMainWindow::cb_harshly_clean_all_jobs(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->clean_all_jobs(true);
}

void GuiMainWindow::clean_all_jobs(bool harshly) {
    std::string describe_action = "clean up intermediate files from the entire pipeline";
    if (harshly) {
        describe_action = "harshly " + describe_action + " \n\n\
Harsh cleaning will remove micrographs, movies and particle stacks from all MotionCorr, Extract, \n\
Polish and Subtract directories. This means you will NOT be able to use those images in subsequent runs anymore, \n\
although you could always recreate the data by continuing the job (possibly at considerable computing costs).\n \n \
You can protect specific jobs from harsh cleaning by creating a file called \"NO_HARSH_CLEAN\" inside their directory,\n\
e.g. by using \"touch Polish/job045/NO_HARSH_CLEAN\". Below is a list of currently protected jobs:\n \n";
        for (const auto &process : pipeline.processList) {
            if (process.status == Process::FINISHED_SUCCESS && (
                process.type == Process::MOTIONCORR || 
                process.type == Process::EXTRACT || 
                process.type == Process::SUBTRACT
            )) {
                if (exists(process.name + "NO_HARSH_CLEAN"))
                    describe_action += process.name + " \n";
            }
        }
    } else {
        describe_action = "gently " + describe_action;
    }

    if (user_wants_to("Clean up", describe_action)) {
        std::cout << (harshly ? "Harshly" : "Gently") << " cleaning all finished jobs ..." << std::endl;

        try {
            pipeline.cleanupAllJobs(harshly);
        } catch (std::string errmsg) {
            fl_message("%s", errmsg.c_str());
        }

        fl_message("Done cleaning! Don't forget the files are all still in the Trash folder. Use the \"Empty Trash\" option from the File menu to permanently delete them.");
    }
}

// Run button callback functions
void GuiMainWindow::cb_gentle_cleanup(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cleanup(-1, true, false);
}

void GuiMainWindow::cb_harsh_cleanup(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cleanup(-1, true, true);
}

void GuiMainWindow::cleanup(int jobindex, bool ask, bool harshly) {
    // Allow cleaning the currently selected job from the GUI
    if (jobindex < 0) {
        if (current_job >= 0) {
            jobindex = current_job;
        } else {
            fl_message("Please select a job.");
            return;
        }
    }

    if (!ask || user_wants_to("Clean up", "clean up intermediate files from " + pipeline.processList[current_job].name)) {
        try {
            pipeline.cleanupJob(jobindex, harshly);
        } catch (std::string errmsg) {
            fl_message("%s", errmsg.c_str());
        }
    }
}

// Run button callback functions
void GuiMainWindow::cb_set_alias(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->set_alias();
}

void GuiMainWindow::set_alias(std::string alias) {
    if (current_job < 0) {
        fl_message("Please select a job.");
        return;
    }

    FileName fn_pre, fn_jobnr, fn_post, fn_dummy, default_ask;
    if (!decomposePipelineFileName(pipeline.processList[current_job].name, fn_pre, fn_jobnr, fn_post))
        REPORT_ERROR("GuiMainWindow::set_alias ERROR: invalid pipeline process name: " + pipeline.processList[current_job].name);

    // Start the asking window with the current alias
    FileName fn_alias = pipeline.processList[current_job].alias;
    if (fn_alias != "None") {
        default_ask = fn_alias.without(fn_pre);
        if (default_ask[default_ask.length() - 1] == '/')
            default_ask = default_ask.beforeLastOf("/");
    } else {
        default_ask = fn_jobnr.beforeLastOf("/");
    }

    while (true) {
        // If the alias already contains a uniquedate string it may be a continuation of a relion_refine job
        // (where alias_current_job contains a different uniqdate than the outputname of the job)
        if (alias.empty() || decomposePipelineFileName(alias, fn_dummy, fn_dummy, fn_dummy)) {
            // if an alias is provided, just check it is unique, otherwise ask
            const char *palias = fl_input("Rename to: ", default_ask.c_str());
            if (!palias) return;
            std::string al2 (palias);  // Direct initialisation
            alias = al2;
        }

        try {
            pipeline.setAliasJob(current_job, alias);
            break;
        } catch (std::string errmsg) {
            alias = "";
            fl_message("%s", errmsg.c_str());
        }
    }
}

// Run button callback functions
void GuiMainWindow::cb_abort(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (current_job < 0) {
        fl_message("Please select a job.");
        return;
    }

    if (pipeline.processList[current_job].status != Process::RUNNING) {
        fl_message("You can only abort running jobs ... ");
        return;
    }

    if (user_wants_to("Abort", "abort job: " + pipeline.processList[current_job].name))
        touch(pipeline.processList[current_job].name + RELION_JOB_ABORT_NOW);
}

// Run button callback functions
void GuiMainWindow::cb_mark_as_finished(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->mark_as_finished(false);
}

// Run button callback functions
void GuiMainWindow::cb_mark_as_failed(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->mark_as_finished(true);
}

void GuiMainWindow::mark_as_finished(bool is_failed) {

    if (current_job < 0) {
        fl_message("You can only mark existing jobs as finished!");
        return;
    }

    try {
        pipeline.markAsFinishedJob(current_job, is_failed);
    } catch (std::string errmsg) {
        fl_message("%s", errmsg.c_str());
        return;
    }

    updateJobLists();
}

// Run button callback functions
void GuiMainWindow::cb_make_flowchart(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (current_job < 0) {
        fl_message("Please select a job.");
        return;
    }

    try {
        pipeline.makeFlowChart(current_job, true);
    } catch (std::string errmsg) {
        fl_message("%s", errmsg.c_str());
        return;
    }

    T->updateJobLists();
}

void GuiMainWindow::cb_edit_note(Fl_Widget*, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_edit_note_i(false);
}

void GuiMainWindow::cb_edit_project_note(Fl_Widget*, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_edit_note_i(true);
}

void GuiMainWindow::cb_edit_note_i(bool is_project_note) {
    FileName fn_note;
    std::string title;
    if (is_project_note) {
        fn_note = "project_note.txt";
        title = "Overall project notes";
    } else {
        if (current_job < 0) {
            fl_message(" You can only edit the note for existing jobs ... ");
            return;
        }
        Process job = pipeline.processList[current_job];
        fn_note = job.name + "note.txt";
        title = job.alias == "None" ? job.name : job.alias;
    }
    NoteEditorWindow* w = new NoteEditorWindow(660, 400, title.c_str(), fn_note, !maingui_do_read_only);
    w->show();
}

// Save button callback function
void GuiMainWindow::cb_save(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->save();
}

void GuiMainWindow::save() {
    // Get which job we're dealing with, and update it from the GUI
    int iwin = browser->value() - 1;
    gui_jobwindows[iwin]->updateMyJob();

    // For scheduled jobs, also allow saving the .job file in the output directory
    if (current_job >= 0 && pipeline.processList[current_job].status == Process::SCHEDULED) {
        gui_jobwindows[iwin]->myjob.write(pipeline.processList[current_job].name);
    }
    // Write the hidden file
    gui_jobwindows[iwin]->myjob.write("");
}

// Load button callback function
void GuiMainWindow::cb_load(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->load();
}

void GuiMainWindow::load() {
    int iwin = browser->value() - 1;
    gui_jobwindows[iwin]->myjob.read("", is_main_continue);
    alias_current_job->value("Give_alias_here");
    gui_jobwindows[iwin]->updateMyGui();

    // Make the current continue-setting active
    cb_toggle_continue_i();
}

// Load button callback function
void GuiMainWindow::cb_undelete_job(Fl_Widget* o, void* v) {
    std::string fn_dir = "./Trash/.";
    std::string fn_filter = "Pipeline STAR files (job_pipeline.star)";
    Fl_File_Chooser chooser (
        fn_dir.c_str(), fn_filter.c_str(), 
        Fl_File_Chooser::SINGLE, "Choose pipeline STAR file to import"
    );
    chooser.show();
    // Block until user picks something.
    while (chooser.shown()) Fl::wait();

    // User hit cancel?
    if (!chooser.value()) return;

    char relname[FL_PATH_MAX];
    fl_filename_relative(relname, sizeof(relname), chooser.value());
    FileName fn_pipe(relname);

    pipeline.undeleteJob(fn_pipe);
}

void GuiMainWindow::cb_export_jobs(Fl_Widget* o, void* v) {
    try {
        // Get the name of this block of exported jobs and make the corresponding directory
        std::string dir = fl_input(std::string("Name of the exported block of jobs? ").c_str(), std::string("export1").c_str());
        pipeline.exportAllScheduledJobs(dir);
    } catch (const std::string &errmsg) {
        fl_message("%s", errmsg.c_str());
    }
}

void GuiMainWindow::cb_import_jobs(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    // Get the directory with the Exported jobs
    std::string fn_dir = ".";
    std::string fn_filter = "Export STAR file (exported.star)";
    Fl_File_Chooser chooser(fn_dir.c_str(),  fn_filter.c_str(), Fl_File_Chooser::SINGLE, "Choose pipeline STAR file to import");
    chooser.show();
    // Block until user picks something.
    while (chooser.shown()) Fl::wait();

    // User hit cancel?
    if (!chooser.value()) return;

    FileName fn_export(chooser.value());

    pipeline.importJobs(fn_export);

    // refresh the joblists
    T->updateJobLists();
}

// Re-order running and finished job lists
void GuiMainWindow::cb_order_jobs_alphabetically(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    do_order_alphabetically = true;
    T->fillRunningJobLists();
}

// Re-order running and finished job lists
void GuiMainWindow::cb_order_jobs_chronologically(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    do_order_alphabetically = false;
    T->fillRunningJobLists();
}

// Empty-trash button callback function
void GuiMainWindow::cb_empty_trash(Fl_Widget* o, void* v) {
    if (user_wants_to("Empty Trash", "remove the entire Trash folder")) {
        std::string command = "rm -rf Trash";
        std::cout << " Executing: " << command << std::endl;
        system(command.c_str());
    }
}

void GuiMainWindow::cb_print_notes(Fl_Widget*, void* v) {

    FileName fn_tmp = pipeline.name + "_all_notes.txt";
    std::ofstream fh (fn_tmp.c_str(), std::ios::out);

    for (const auto &process : pipeline.processList) {
        FileName fn_note = process.name+"note.txt";
        fh << " ################################################################ " << std::endl;
        fh << " # Job= " << process.name;
        if (process.alias != "None")
            fh << " alias: " << process.alias;
        fh << std::endl;
        if (exists(fn_note)) {
            std::ifstream in (fn_note.data(), std::ios_base::in);
            if (in.fail())
            REPORT_ERROR((std::string) "ERROR: cannot read file " + fn_note);
            in.seekg(0);
            std::string line;
            while (getline(in, line, '\n'))
                fh << line << std::endl;
        }
    }

    fl_message("Done writing all notes into file: %s", fn_tmp.c_str());
}

void GuiMainWindow::cb_remake_nodesdir(Fl_Widget*, void* v) {
    pipeline.remakeNodeDirectory();
}

void GuiMainWindow::cb_reread_pipeline(Fl_Widget*, void* v) {
    PipeLine::rwlock lock (pipeline, " mainGUI reread_pipeline_i");
    // With the locking system, each read needs to be followed soon with a write
}

void GuiMainWindow::cb_reactivate_runbutton(Fl_Widget* o, void* v) {
    run_button->activate();
}

void GuiMainWindow::cb_toggle_overwrite_continue(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    do_overwrite_continue = !do_overwrite_continue;
    T->cb_toggle_continue_i();
}

void GuiMainWindow::cb_show_initial_screen(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_show_initial_screen_i();
}

void GuiMainWindow::cb_show_initial_screen_i() {
    run_button->deactivate();
    cb_select_browsegroup_i(true);
}

void GuiMainWindow::cb_toggle_pipeliner_scheduler(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_toggle_pipeliner_scheduler_i();
}

void GuiMainWindow::cb_toggle_pipeliner_scheduler_i() {
    if (show_scheduler) {
        pipeliner_grp->hide();
        scheduler_grp->show();

        // If this schedule is running, then use the I/O viewer, otherwise use Jobs viewer
        show_expand_stdout = !schedule.isWriteLocked();

    } else {
        scheduler_grp->hide();
        pipeliner_grp->show();
        // After toggling, always go back to non-expanded view
        show_expand_stdout = true;
    }
    toggle_expand_stdout();
}

void GuiMainWindow::cb_create_schedule(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (!show_scheduler) return;

    FileName fn_new = (std::string) fl_input("%s", "", std::string(" Name of the new schedule: ").c_str());
    if (fn_new.length() > 0) {
        T->cb_toggle_schedule_i(false, fn_new);
    }
}

void GuiMainWindow::cb_copy_schedule(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;

    if (!show_scheduler) return;

    FileName fn_copy = std::string(fl_input(
        "%s", "", std::string(" Name of the copy schedule: ").c_str()
    ));

    if (fn_copy.length() == 0) return;

    schedule.copy("Schedules/" + fn_copy);
    T->cb_toggle_schedule_i(false, fn_copy);
}

void GuiMainWindow::cb_toggle_schedule(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_toggle_schedule_i(false);
}

void GuiMainWindow::cb_toggle_pipeline(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->cb_toggle_schedule_i(true);
}

void GuiMainWindow::cb_toggle_schedule_i(bool do_pipeline, FileName fn_new_schedule) {
    if (do_pipeline) {
        // Read in the pipeline STAR file if it exists
        pipeline.name = "default";
        show_scheduler = false;
    } else {
        show_scheduler = true;
        FileName fn_sched;
        if (fn_new_schedule != "") {
            fn_sched = "Schedules/" + fn_new_schedule;
            // Also add entry to the menu
            std::string mylabel = "Schedules/Schedules/" + fn_new_schedule;
            menubar->add(mylabel.c_str(), 0, cb_toggle_schedule, this);
        } else {
            fn_sched = "Schedules/" + std::string(menubar->text());
        }
        schedule.setName(fn_sched + "/");
        pipeline.name = fn_sched + "/schedule";

        // Read in scheduler or create new one if it did not exist
        if (exists(schedule.name+"schedule.star")) {
            schedule.read(DONT_LOCK);
            pipeline.name = fn_sched + "/schedule";
        } else {
            system(("mkdir -p " + fn_sched).c_str());
            schedule.write(DONT_LOCK);  // empty write
        }
        fillStdOutAndErr();
    }

    if (exists(pipeline.name + "_pipeline.star")) {
        PipeLine::rwlock lock (pipeline, "mainGUI constructor");
        // With the locking system, each read needs to be followed soon with a write
    } else {
        pipeline.write();
    }

    cb_toggle_pipeliner_scheduler_i();

    updateJobLists();
}

void GuiMainWindow::cb_start_pipeliner(Fl_Widget* o, void* v) {
    // GuiMainWindow *T = (GuiMainWindow*) v;
    // GuiMainWindow *T = (GuiMainWindow*) v;
    std::vector<FileName> job_names;

    for (long int id : scheduled_processes) {
        job_names.push_back(pipeline.processList[id].name);
    }
    SchedulerWindow *w = new SchedulerWindow(400, 300, "Select which jobs to execute");
    w->fill(pipeline.name, job_names);
}

void GuiMainWindow::cb_stop_pipeliner(Fl_Widget* o, void* v) {
    // GuiMainWindow *T = (GuiMainWindow*) v;
    Fl_File_Chooser chooser(
        ".",
        std::string("Pipeline scheduled file (RUNNING_PIPELINER_" + pipeline.name + "_*)").c_str(),
        Fl_File_Chooser::SINGLE, "Choose which scheduler to stop"
    );
    chooser.show();
    // Block until user picks something.
    while (chooser.shown()) Fl::wait();

    // User hit cancel?
    if (!chooser.value()) return;

    FileName fn_del(chooser.value());
    std::cout << " Deleting file : " << fn_del << std::endl;
    std::remove(fn_del.c_str());
}

void GuiMainWindow::cb_toggle_expand_stdout(Fl_Widget* o, void* v) {
    GuiMainWindow *T = (GuiMainWindow*) v;
    T->toggle_expand_stdout();
}

void GuiMainWindow::toggle_expand_stdout() {
    Fl_Group *jobs_grp = show_scheduler ? scheduler_jobs_grp : pipeliner_jobs_grp;
    if (show_expand_stdout) {
        expand_stdout_grp->hide();
        jobs_grp->show();
        expand_stdout_button->label("I/O view");
    } else {
        expand_stdout_grp->show();
        jobs_grp->hide();
        expand_stdout_button->label("Job view");
    }
    show_expand_stdout = !show_expand_stdout;
}

constexpr const char* about() {
    // Compile-time #include hack
    return "RELION " RELION_SHORT_VERSION "\n\n"
    #include "src/help.txt"
    ;
}

void GuiMainWindow::cb_about(Fl_Widget* o, void* v) {
    ShowHelpText help = ShowHelpText(about());
}

void GuiMainWindow::cb_quit(Fl_Widget* o, void* v) {
    exit(0);
}
