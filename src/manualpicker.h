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

#ifndef MANUALPICKER_H_
#define MANUALPICKER_H_

// This #define / #undef pair protects against another Complex definition in fltk.
#define Complex
#include <FL/Fl.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_Box.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Text_Display.H>
#undef Complex

#include "src/metadata_table.h"
#include "src/args.h"
#include "src/funcs.h"
#include "src/filename.h"
#include "src/gui_entries.h"
#include "src/jaz/obs_model.h"

const int MWCOL1 = 300;
const int MWCOL2 = 60;
const int MWCOL3 = 60;
const int MWCOL4 = 80;
const int MXCOL0 = 30;
const int MXCOL1 = MXCOL0 + MWCOL1 + 10;
const int MXCOL2 = MXCOL1 + MWCOL2 + 10;
const int MXCOL3 = MXCOL2 + MWCOL3 + 10;
const int MXCOL4 = MXCOL3 + MWCOL4 + 10;
const int TOTALWIDTH = MWCOL1 + MWCOL2 + MWCOL3 + MWCOL4 + MWCOL4 + 100;
const int TOTALHEIGHT = 500;

// The button for picking particles
void cb_viewmic(Fl_Widget *w, void *data);
// The button for viewing the CTF
void cb_viewctf(Fl_Widget *w, void *data);
// The selection button
void cb_selectmic(Fl_Widget *w, void *data);

// This class only puts scrollbars around the resizable canvas
class manualpickerGuiWindow: public Fl_Window {

    public:

    // Input, picking & output names
    FileName fn_in, fn_sel;

    // Allow saving selected micrographs?
    bool do_allow_save;

    // Save default selection immediately? (useful for always generating output files in pipeline)
    bool do_fast_save;

    // MetaDataTable of input micrographs
    MetaDataTable MDin;

    // Observation model of input micrographs
    ObservationModel obsModel;

    // Constructor with w x h size of the window and a title
    manualpickerGuiWindow(int W, int H, const char *title = 0): Fl_Window(W, H, title){}

    // Fill the window with all entries
    int fill();

    private:

    static void cb_menubar_save(Fl_Widget*, void*);
    inline void cb_menubar_save_i();

    static void cb_menubar_invert_selection(Fl_Widget*, void*);
    inline void cb_menubar_invert_selection_i();

    static void cb_menubar_quit(Fl_Widget*, void*);
    inline void cb_menubar_quit_i();

    static void cb_menubar_recount(Fl_Widget*, void*);
    inline void cb_menubar_recount_i();

    void readOutputStarfile();
    void writeOutputStarfile();

};

class ManualPicker {

    public:

    // I/O Parser
    IOParser parser;

    // The input micrographs
    MetaDataTable MDin;

    // Observation model for the input mirographs
    ObservationModel obsModel;

    // Input, picking & output names
    FileName fn_in, fn_sel;

    // Allow save selected micrographs?
    bool do_allow_save;

    // Save an output selection file immediately (with all micrographs selected)
    bool do_fast_save;

    public:

    // Read command line arguments
    void read(int argc, char **argv);

    // Print usage instructions
    void usage();

    // Initialise some general stuff after reading
    void initialise();

    // General function to decide what to do
    void run();

    private:

    void writeOutput();

};

#endif /* MANUALPICKER_H_ */
