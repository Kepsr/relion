import os
import math
import traceback
try:
    import Tkinter as tk
    import tkMessageBox
    import tkFileDialog
except ImportError:
    # The GUI is optional.
    # If the user requests it,
    # it will fail when it tries to open
    # so we can ignore the error for now.
    pass

OPTIONS_FILE = 'relion_it_options.py'


def prefix_RELION_IT(msg):
    return ' RELION_IT: ' + msg


def captureException(method):
    def newmethod():
        try:
            return method()
        except Exception as exc:
            tkMessageBox.showerror("Error", exc.message)
            traceback.print_exc()
    return newmethod


class RelionItGui(object):

    @classmethod
    def launch(cls, options):
        print(prefix_RELION_IT('launching GUI...'))
        tk_root = tk.Tk()
        tk_root.title("relion_it.py setup")
        cls(tk_root, options)
        tk_root.mainloop()

    def __init__(self, main_window, options):
        self.main_window = main_window
        self.options = options
        self.warnings = []

        ### Create GUI

        main_frame = tk.Frame(main_window)
        main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame, right_frame = 2 * (tk.Frame(main_frame),)
        for frame, side in (
            (left_frame,  tk.LEFT),
            (right_frame, tk.RIGHT),
        ):
            frame.pack(side=side, anchor=tk.N, fill=tk.X, expand=1)

        def makeLabelFrame(text, side):
            frame = tk.LabelFrame(side, text=text, padx=5, pady=5)
            frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
            tk.Grid.columnconfigure(frame, 1, weight=1)
            return frame

        def makeEntry(frame, option, labeltext, row):
            tk.Label(frame, text=labeltext).grid(row=row, sticky=tk.W)
            entry = tk.Entry(frame)
            entry.grid(row=row, column=1, sticky=tk.EW)
            entry.insert(0, str(option))
            return entry

        def makeEntryWithVar(frame, option, labeltext, row, vartype=tk.StringVar, wide=False):
            '''
            Return `entry, var`.
            '''
            tk.Label(frame, text=labeltext).grid(row=row, sticky=tk.W)
            var = vartype()  # for data binding
            entry = tk.Entry(frame, textvariable=var)
            if wide:
                entry.grid(row=row, column=1, sticky=tk.EW, columnspan=2)
            else:
                entry.grid(row=row, column=1, sticky=tk.EW)
            entry.insert(0, str(option))
            return entry, var

        def makeEntryWithVarAndAttr(frame, option, labeltext, row, attrtext, vartype=tk.StringVar):
            '''
            Return `entry, var, attr`.
            '''
            tk.Label(frame, text=labeltext).grid(row=row, sticky=tk.W)
            var = vartype()  # for data binding
            entry = tk.Entry(frame, textvariable=var)
            entry.grid(row=row, column=1, sticky=tk.EW)
            entry.insert(0, str(option))
            attr = tk.Label(frame, text=attrtext)
            attr.grid(row=row, column=2, sticky=tk.W)
            return entry, var, attr

        # Convenience function for making file browser buttons
        def makeBrowseButton(frame, var, filetypes):
            def command():
                chosen_file = tkFileDialog.askopenfilename(filetypes=filetypes)
                if chosen_file is not None:
                    # Make path relative if it's in the current directory
                    if chosen_file.startswith(os.getcwd()):
                        chosen_file = os.path.relpath(chosen_file)
                    var.set(chosen_file)
            return tk.Button(frame, text="Browse...", command=command)

        def makeButtonWithVar(frame, option, labeltext, row):
            tk.Label(frame, text=labeltext).grid(row=row, sticky=tk.W)
            var = tk.BooleanVar()
            button = tk.Checkbutton(frame, var=var)
            button.grid(row=row, column=1, sticky=tk.W)
            if option:
                button.select()
            return button, var

        ### Experiment frame

        expt_frame = makeLabelFrame("Experimental details", left_frame)

        self.voltage_entry = makeEntry(
            expt_frame, options.voltage, "Voltage (kV):", 0,
        )

        self.cs_entry = makeEntry(
            expt_frame, options.Cs, "Cs (mm):", 1,
        )

        self.phaseplate_button, self.phaseplate_var = makeButtonWithVar(
            expt_frame, options.ctffind_do_phaseshift, "Phase plate?", 2,
        )

        self.angpix_entry, self.angpix_var = makeEntryWithVar(
            expt_frame, options.angpix, u"Pixel size (\u212B):", 3, vartype=tk.DoubleVar,
        )

        self.exposure_entry = makeEntry(
            expt_frame, options.motioncor_doseperframe,
            u"Exposure rate (e\u207B / \u212B\u00B2 / frame):", 4,
        )

        ### Particle frame

        particle_frame = makeLabelFrame("Particle details", left_frame)

        (
            (self.particle_max_diam_entry, self.particle_max_diam_var),
            (self.particle_min_diam_entry, self.particle_min_diam_var),
        ) = (
            makeEntryWithVar(
                particle_frame, option, text, row, vartype=tk.DoubleVar, wide=True,
            ) for row, (option, text) in enumerate((
                (options.autopick_LoG_diam_max, u"Longest diameter (\u212B):"),
                (options.autopick_LoG_diam_min, u"Shortest diameter (\u212B):"),
            ))
        )

        self.ref_3d_entry, self.ref_3d_var = makeEntryWithVar(
            particle_frame, options.autopick_3dreference, "3D reference (optional):", 2,
        )
        makeBrowseButton(
            particle_frame, self.ref_3d_var,
            filetypes=(('MRC file', '*.mrc'), ('All files', '*'))
        ).grid(row=2, column=2)

        (
            (self.mask_diameter_entry,         self.mask_diameter_var,         self.mask_diameter_px),
            (self.box_size_entry,              self.box_size_var,              self.box_size_in_angstrom),
            (self.extract_small_boxsize_entry, self.extract_small_boxsize_var, self.extract_angpix),
            #                                                                  ^ will be used if `self.angpix_entry.get()` fails
        ) = (
            makeEntryWithVarAndAttr(
                particle_frame, option, labeltext, row+3, attrtext, vartype=tk.DoubleVar,
            ) for row, (option, labeltext, attrtext) in enumerate((
                (options.mask_diameter,            u"Mask diameter (\u212B):", "= NNN px"        ),
                (options.extract_boxsize,          "Box size (px):",           u"= NNN \u212B"   ),
                (options.extract_small_boxsize[0], "Down-sample to (px):",     u"= NNN \u212B/px"),
            ))
        )

        self.auto_boxsize_button, self.auto_boxsize_var = makeButtonWithVar(
            particle_frame, True, "Calculate for me:", 6
        )

        ### Project frame

        project_frame = makeLabelFrame("Project details", right_frame)

        tk.Label(project_frame, text="Project directory:").grid(
            row=0, sticky=tk.W
        )
        tk.Label(project_frame, text=os.getcwd(), anchor=tk.W).grid(
            row=0, column=1, sticky=tk.W, columnspan=2
        )

        self.import_images_entry, self.import_images_var = makeEntryWithVar(
            project_frame, self.options.import_images, "Pattern for movies:", 1
        )
        makeBrowseButton(
            project_frame, self.import_images_var,
            filetypes=(('Image file', '{*.mrc, *.mrcs, *.tif, *.tiff}'), ('All files', '*'))
        ).grid(row=1, column=2)

        self.gainref_entry, self.gainref_var = makeEntryWithVar(
            project_frame, self.options.motioncor_gainreference, "Gain reference (optional):", 2
        )
        makeBrowseButton(
            project_frame, self.gainref_var,
            filetypes=(('MRC file', '*.mrc'), ('All files', '*'))
        ).grid(row=2, column=2)

        ### Pipeline frame

        pipeline_frame = makeLabelFrame("Pipeline control", right_frame)

        (
            (self.stop_after_ctf_button, self.stop_after_ctf_var),
            (self.class2d_button,        self.class2d_var),
            (self.class3d_button,        self.class3d_var),
            (self.second_pass_button,    self.second_pass_var),
            (self.class2d_pass2_button,  self.class2d_pass2_var),
            (self.class3d_pass2_button,  self.class3d_pass2_var),
        ) = (
            makeButtonWithVar(
                pipeline_frame, option, text, row
            ) for row, (option, text) in enumerate((
                (options.stop_after_ctf_estimation, "Stop after CTF estimation?"),
                (options.do_class2d[0],             "Do 2D classification?"),
                (options.do_class3d[0],             "Do 3D classification?"),
                (options.do_pass[1],                "Do second pass (only if no 3D ref)?"),
                (options.do_class2d[1],             "Do 2D classification (2nd pass)?"),
                (options.do_class3d[1],             "Do 3D classification (2nd pass)?"),
            ))
        )

        ### Add logic to the box size boxes

        for var in (self.box_size_var, self.extract_small_boxsize_var):
            var.trace('w', self.update_box_size_labels)

        for var in (self.angpix_var, self.particle_max_diam_var):
            var.trace('w', self.update_box_sizes)

        self.auto_boxsize_button.config(command=self.update_box_sizes)

        ### Add logic to the check boxes

        for button in (
            self.stop_after_ctf_button, self.class3d_button, self.second_pass_button,
        ):
            button.config(command=self.update_pipeline_control_state)

        self.ref_3d_var.trace('w', self.update_pipeline_control_state)

        ###

        padding = 5

        button_frame = tk.Frame(right_frame)
        button_frame.pack(padx=padding, pady=padding, fill=tk.X, expand=1)

        self.run_button = tk.Button(button_frame, text="Save & run", command=self.run_pipeline)
        self.run_button.pack(padx=padding, pady=padding, side=tk.RIGHT)

        self.save_button = tk.Button(button_frame, text="Save options", command=self.save_options)
        self.save_button.pack(padx=padding, pady=padding, side=tk.RIGHT)

        # Show initial pixel sizes
        self.update_box_sizes()

    def update_box_size_labels(self, *args, **kwargs):
        def handleboxsize():
            self.box_size_in_angstrom.config(text=u"= NNN \u212B")
        def handlemaskdiam():
            self.mask_diameter_px.config(text="= NNN px")
        def handleextractangpix():
            self.extract_angpix.config(text=u"= NNN \u212B/px")
        try:
            angpix = self.angpix_entry.get()
        except ValueError:
            # Can't update any of the labels without angpix
            handlemaskdiam()
            handleboxsize()
            handleextractangpix()
            return
        try:
            mask_diameter_px = self.mask_diameter_entry.get() / angpix
            self.mask_diameter_px.config(text="= {:.1f} px".format(mask_diameter_px))
        except (ValueError, ZeroDivisionError):
            handlemaskdiam()
            # Don't return - an error here doesn't stop us calculating the other labels
        try:
            box_angpix = angpix * self.box_size_entry.get()
            self.box_size_in_angstrom.config(text=u"= {:.1f} \u212B".format(box_angpix))
        except ValueError:
            # Can't update these without the box size
            handleboxsize()
            handleextractangpix()
            return
        try:
            small_box_angpix = box_angpix / self.extract_small_boxsize_entry.get()
            self.extract_angpix.config(text=u"= {:.3f} \u212B/px".format(small_box_angpix))
        except (ValueError, ZeroDivisionError):
            # Can't update the downscaled pixel size unless the downscaled box size is valid
            handleextractangpix()

    def update_box_sizes(self, *args, **kwargs):
        # Always activate entry boxes - either we're activating them anyway, or we need to edit the text.
        # For text editing we need to activate the box first then deactivate again afterwards.
        entries = (
            self.mask_diameter_entry,
            self.box_size_entry,
            self.extract_small_boxsize_entry,
        )
        for entry in entries:
            entry.config(state=tk.NORMAL)
        if self.auto_boxsize_var.get():
            try:
                angpix = self.angpix_entry.get()
                particle_size_angstroms = self.particle_max_diam_entry.get()
                particle_size_pixels = particle_size_angstroms / angpix
                mask_diameter = 1.1 * particle_size_angstroms
                box_size = calculate_box_size(particle_size_pixels)
                small_boxsize = calculate_downscaled_box_size(box_size, angpix)

                for entry, datum in zip(entries, (
                    mask_diameter,
                    box_size,
                    small_boxsize,
                )):
                    entry.delete(0, tk.END)
                    entry.insert(0, str(datum))

            except:
                # Ignore errors - they will be picked up if the user tries to save the options
                pass
            for entry in entries:
                entry.config(state=tk.DISABLED)
        self.update_box_size_labels()

    def update_pipeline_control_state(self, *args, **kwargs):
        new_state = tk.DISABLED if self.stop_after_ctf_var.get() else tk.NORMAL
        for button_or_entry in (
            self.class2d_button,
            self.class3d_button,
            self.particle_max_diam_entry,
            self.particle_min_diam_entry,
            self.ref_3d_entry,
            self.auto_boxsize_button,
        ):
            button_or_entry.config(state=new_state)
        # Update the box size controls with care to avoid activating them when we shouldn't
        if new_state == tk.DISABLED:
            for entry in (
                self.mask_diameter_entry, self.box_size_entry, self.extract_small_boxsize_entry,
            ):
                entry.config(state=new_state)
        else:
            self.update_box_sizes()
        can_do_second_pass = (
            self.class3d_var.get()
            and not self.ref_3d_var.get()
            and not self.stop_after_ctf_var.get()
        )
        self.second_pass_button.config(state=tk.NORMAL if can_do_second_pass else tk.DISABLED)
        will_do_second_pass = can_do_second_pass and self.second_pass_var.get()
        self.class2d_pass2_button.config(state=tk.NORMAL if will_do_second_pass else tk.DISABLED)
        self.class3d_pass2_button.config(state=tk.NORMAL if will_do_second_pass else tk.DISABLED)

    def fetch_options_from_gui(self):
        '''
        Fetch the current values from the GUI widgets and store them in the options object.

        Returns:
            A list of warning messages about possible incorrect option values.

        Raises:
            ValueError: If an option value is invalid.
        '''

        self.warnings = []
        opts = self.options

        opts.stop_after_ctf_estimation = self.stop_after_ctf_var.get()
        opts.do_class2d                = self.class2d_var.get(), self.class2d_pass2_var.get()
        opts.do_class3d                = self.class3d_var.get(), self.class3d_pass2_var.get()
        opts.do_pass[1]                = self.second_pass_var.get()

        def try_get(numtype, attribute, name):
            '''
            Try to return `numtype(attribute.get())`
            (where `numtype` is either `float` or `int`).
            '''
            try:
                return numtype(attribute.get())
            except ValueError:
                raise ValueError("{} must be a number".format(name))

        def check_positive(option, name):
            '''
            Add warning if `option` is not positive.
            '''
            if option <= 0.0:
                self.warnings.append("- {} should be a positive number".format(name))

        def try_get_diam(numtype, attribute, name):
            attr = attribute.get()
            try:
                return numtype(attr)
            except ValueError:
                if len(attr) == 0 and opts.stop_after_ctf_estimation:
                    # This was left blank and won't be used, set to zero to avoid errors in calculations later
                    return 0.0
                else:
                    raise ValueError("{} must be a number".format(name))

        def check_file_exists(filename, name):
            if filename and not os.path.isfile(filename):
                self.warnings.append("- {} '{}' does not exist".format(name, filename))

        opts.voltage = try_get(float, self.voltage_entry, "Voltage")
        check_positive(opts.voltage, "Voltage")

        opts.Cs = try_get(float, self.cs_entry, "Cs")

        opts.ctffind_do_phaseshift = self.phaseplate_var.get()

        opts.angpix = try_get(float, self.angpix_entry, "Pixel size")  # XXX This might be redundant, given `self.angpix_entry.get()` should now have type `float`.
        check_positive(opts.angpix, "Pixel size")

        opts.motioncor_doseperframe = try_get(float, self.exposure_entry, "Exposure rate")
        check_positive(opts.motioncor_doseperframe, "Exposure rate")

        opts.autopick_LoG_diam_max = try_get_diam(float, self.particle_max_diam_entry, "Particle longest diameter")
        opts.autopick_LoG_diam_min = try_get_diam(float, self.particle_min_diam_entry, "Particle shortest diameter")

        opts.autopick_3dreference = self.ref_3d_entry.get()
        check_file_exists(opts.autopick_3dreference, "3D reference file")

        opts.mask_diameter = try_get(float, self.mask_diameter_entry, "Mask diameter")
        check_positive(opts.mask_diameter, "Mask diameter")

        opts.extract_boxsize = try_get(int, self.box_size_entry, "Box size")
        check_positive(opts.extract_boxsize, "Box size")

        opts.extract_small_boxsize = 2 * (try_get(int, self.extract_small_boxsize_entry, "Down-sampled box size"),)
        opts.extract_downscale = True, True
        check_positive(opts.extract_small_boxsize[0], "Down-sampled box size")

        opts.import_images = self.import_images_entry.get()
        if opts.import_images.startswith(('/', '..')):
            self.warnings.append("- Movies should be located inside the project directory")
        if '*' not in opts.import_images:
            self.warnings.append("- Pattern for input movies should normally contain a '*' to select more than one file")

        opts.motioncor_gainreference = self.gainref_entry.get()
        check_file_exists(opts.motioncor_gainreference, "Gain reference file")

    def calculate_full_options(self):
        '''
        Update the options from the values that have been fetched from the GUI.

        This method uses the values that the user has set in the GUI to calculate a number of other options for the
        script.
        '''
        opts = self.options

        # No 3D reference - do LoG autopicking in the first pass
        opts.autopick_do_LoG = not opts.autopick_3dreference
        opts.class3d_reference = opts.autopick_3dreference
        # If we have a 3D reference, do a single pass with a large batch size
        if opts.autopick_3dreference:
            opts.autopick_refs_min_distance = opts.autopick_LoG_diam_max * 0.7
            opts.do_pass[1] = False

        # Now set a sensible batch size (leaving the batch size for the second pass at its default 100,000)
        opts.batchsize[0] = 10000 if opts.do_pass[1] else 100000

    @captureException
    def save_options(self):
        '''
        Update the full set of options from the values in the GUI, and save them to a file.

        Returns:
            True if the options were valid and saved successfully, otherwise False.
        '''
        self.fetch_options_from_gui()
        if not self.warnings or tkMessageBox.askokcancel(
            'Warning', '\n'.join(self.warnings),
            icon='warning', default=tkMessageBox.CANCEL,
        ):
            self.calculate_full_options()
            print(prefix_RELION_IT("Writing all options to {}".format(OPTIONS_FILE)))
            self.options.write_to(OPTIONS_FILE)
            return
        raise Exception()

    def run_pipeline(self):
        # Update options from GUI, close GUI, run pipeline.
        try:
            self.save_options()
            self.main_window.destroy()
            self.options.run_pipeline()
        except Exception:
            pass


def calculate_box_size(particlesize):
    '''Return an even side length 20% larger than particle.'''
    boxsize = int(math.ceil(1.2 * particlesize))
    return boxsize + boxsize % 2


def calculate_downscaled_box_size(box_size_pix, angpix):
    for small_box_pix in (
        48, 64, 96, 128, 160, 192, 256, 288, 300, 320, 360,
        384, 400, 420, 450, 480, 512, 640, 768, 896, 1024,
    ):
        # Don't go larger than the original box
        if small_box_pix > box_size_pix:
            return box_size_pix
        # If Nyquist freq. is better than 8.5 A, use this downscaled box,
        # otherwise continue to next size up
        if angpix * box_size_pix / small_box_pix < 4.25:
            return small_box_pix
    # Fall back to a warning message
    raise Exception("Box size is too large!")

