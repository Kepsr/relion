#!/usr/bin/env python2.7
'''
relion_it.py
============

Script for automated, on-the-fly single-particle analysis in RELION 3

Authors: Sjors H.W. Scheres, Takanori Nakane & Colin M. Palmer

Usage:
    relion_it.py  [extra_options.py [extra_options2.py ....] ] [--gui] [--continue]

To get started, go to the intended location of your RELION project directory and make sure your micrographs are
accessible from within it (e.g. in a subdirectory called `Movies/' - use a symlink if necessary). Then run this
script, providing the names of files containing options if needed. (To call the script, you'll need to enter the full
path to it, put the directory containing it on your PATH environment variable, or put a copy of the script in the
current directory.)

Run with the `--gui' option to launch a simple GUI which will set up a run from a few basic options. (The GUI can
also be used to save a complete options file that you can then edit as required.)

Once the script is running, open a normal RELION GUI to see what's happening and visualise the results.

See below for full instructions including how to handle errors. If you have any problems, please edit the script as
needed, call on your local Python expert or email the CCP-EM mailing list (https://www.jiscmail.ac.uk/ccpem).


Overview
--------

relion_it.py creates a number of RELION jobs and then runs one or more `relion_pipeliner' processes to schedule them
(exactly like using the "Schedule" button in the RELION GUI). Instructions and information are printed to the terminal
by relion_it.py as it runs.

relion_it.py uses a large number of options to control how the jobs are run. It's designed to be very flexible and so
these options can be changed in a number of ways:

- The easiest way is to use the simple GUI (enabled by passing the `--gui' argument), which allows you to set a few
  simple options. These are then used to calculate appropriate values for the complete set of options. (See "Using the
  GUI" below for more information on this.)

- For more control, options can be put into one or more Python files (with a simple "option_name = value" format or
  with more complicated calculations - see "Options files" below for more information). The names of these options
  files can passed as command line arguments to relion_it.py.

- For maximum control, you can make your own copy of this script and change the option values and the code itself
  however you want.

Before running relion_it.py, you need to make sure you're in your intended RELION project directory, and that your
movie files are accessible by relative paths within that directory (as usual for a RELION project). You could do this
by moving the files from the microscope straight into the project directory, using a symlink from your project
directory to the real location of the data, or running a script to create a new symlink to each micrograph as it is
collected.


Options files
-------------

relion_it.py uses a large number of options for controlling both the flow of the script and the parameters for
individual jobs. These options can be read from Python script files when relion_it.py is started.

The options are all listed the body of the script below, with a comment to explain each option. One way to use this
script is to copy it in its entirety into your project directory, edit the options directly in the script and then
run it (with no command line arguments). However, it's often better to keep the script in the RELION source directory
(where it can be updated easily) and use options files to configure it.

An example of a simple options file is:

angpix = 1.06

This would override the default pixel size value, but leave all other options at their defaults.

The options files are read and interpreted as Python scripts. A simple list of "option_name = value" lines is all
that is needed, though you can also use any Python commands you like to do more complex calculations. To generate
an example file containing all of the options, run "relion_it.py --gui" and then click the "Save options" button,
which will save all the current options to a file called `relion_it_options.py' in the working directory.

The options are named descriptively so you can probably understand what most of them do quite easily. For more help on
any particular option, look at the comment above its definition in this script, or search the script's code to see
how it is used.

Options files can be useful as templates. As an example, at Diamond Light Source's eBIC facility, we have a template
file called `dls_cluster_options.py' that contains the necessary settings to make relion_it.py submit most of its jobs
to run on the DLS GPU cluster. You could also set up standard templates for a particular microscope (say, voltage and
Cs settings) or for a particular project or computer configuration.

When relion_it.py starts, it reads all options files in the order they are given on the command line. Subsequent files
will override earlier ones, so the last value given for any particular option will be the value that is used.

If you start relion_it.py with the `--continue' argument, it will automatically add `relion_it_options.py' to the end
of the list of options files. This means that if you are in a project directory where the relion_it.py GUI has
previously been used, all options will be defined in the relion_it_options.py file and they will override any other
options files given on the command line. (This is very useful for restarting the script after a problem, but it would
be pointless to combine `--continue' with any options template files.)

Note that if relion_it.py finds option names that it doesn't recognise while it's reading an options file, it will
print a warning (but continue anyway). If you've been editing options files by hand, you should check the output from
relion_it.py when it starts to make sure there are no typos in the options you wanted to set. (If you're using local
variables for intermediate Python calculations in an options file, it's a good idea to use names starting with a
leading underscore so you can immediately tell them apart from warnings about genuine spelling mistakes.)


Using the GUI
-------------

The GUI provides a simple way to start new projects with relion_it.py. If you want to use it, prepare your project
directory as described above, then start the GUI with "relion_it.py --gui". (If you're using any template options
files, you can give those too, for example "relion_it.py /path/to/site/options.py --gui".)

The window that appears should be self-explanatory. Fill in the options as needed for your project, and use the check
boxes on the right to control what processing steps will be done. When you're ready, click either "Save options" or
"Save & run". The program will check the values you've entered and then use them to calculate a few extra options for
relion_it.py. The options will then be saved to a file called `relion_it_options.py', and if you clicked "Save & run"
the processing run will start immediately.

If any of the entered values are invalid (for example, if there are letters in a field which should be a number), the
GUI will display a message box with an error when you click one of the buttons. It will also display a warning if any
values appear to be incorrect (but you can choose to ignore the warning by clicking "OK").

The GUI will try to calculate some extra options from the values you enter using the following rules:

1. If a 3D reference is given, use a single pass with reference-based autopicking, minimum distance between particles
   of 0.7 times the particle size, and a batch size of 100,000 particles.

2. If no 3D reference is given, run a first pass with reference-free LoG autopicking and a batch size of 10,000, and
   then a second pass with reference-based autopicking and a batch size of 100,000.

These options should be sensible in many cases, but if you'd like to change them, save the options from the GUI using
the "Save options" button, close the GUI, and edit the `relion_it_options.py' file to change the option values as
needed. You can then start the processing run with "relion_it.py --continue".


Running the pipelines
---------------------

relion_it.py uses several different scheduling pipelines to run its jobs. While each one is running, a file is
created in the project directory called `RUNNING_PIPELINER_<name>'. A log of the jobs run by that pipeline is stored
in `pipeline_<name>.log'.

If you want to stop one of the pipelines for any reason, delete its `RUNNING_' file and within a minute or two the
pipeliner will notice that the file has been removed and stop.

relion_it.py itself uses a similar file called `RUNNING_RELION_IT', and you can delete this to stop the script (which
will not affect any pipelines that are already running). It keeps a list of all of the jobs it has submitted in a
file called `RELION_IT_SUBMITTED_JOBS'. This file can be edited manually if necessary (but not while the script is
running!) Most of the jobs are run by the `preprocessing' pipeline. This will do the following:

  1. Import movies
  2. Motion correction
  3. CTF estimation
  4. Particle auto-picking
  5. Particle extraction
  6. Batch selection

After a number of particles have been extracted (1,000 by default), a 2D classification job will be run to provide
feedback on the quality of the data collection and particle picking.

Particles are split into batches of a fixed size (default 10,000 for the first pass with no reference, or 100,000
otherwise). The first batch is special: as it grows, the 2D classification job is re-run repeatedly to provide early
feedback on the quality of the data. For subsequent batches, the script waits for each batch to be complete before
running 2D classification on it.

You can provide reference structures for auto-picking and 3D classification. (If you provide a 3D reference in the
GUI it will automatically be used for both tasks.)

If you do not provide a reference for auto-picking, reference-free LoG picking will be used. If you do not provide a
reference for classification, relion_it.py will run the preprocessing pipeline twice. In the first pass, an initial
model will be generated, and then a second pass of preprocessing will be done using the initial model as a reference
for auto-picking and classification.

relion_it.py makes an effort to try to identify a suitable reference to use from the classes produced by the
InitialModel job, but if it selects an inappropriate reference, you can change it by stopping the pipelines and
script ("rm RUNNING_*"), updating the reference filename stored in the file named `RELION_IT_2NDPASS_3DREF', deleting
the relevant jobs (`autopick2_job' and those following) from the `RELION_IT_SUBMITTED_JOBS' file, then restarting the
pipeline with "relion_it.py --continue".


Fixing problems
---------------

One-off job failure
```````````````````
Occasionally, a single job can fail with an isolated error, for example if there are temporary network problems while
working on a remote filesystem. If this happens, RELION will wait forever for the files to appear that would indicate
the job has finished. In the meantime, no new jobs will be run, which can cause a backlog of micrographs to build up.

To fix this (for a preprocessing job), you can just try to re-run the job from the RELION GUI. Select the job in the
"Running jobs" list, then click "Job actions" -> "Mark as finished". Select the job again in the "Finished jobs"
list, then click "Continue!" to re-start the job.

That approach should work for preprocessing jobs, but probably won't work for classification or inital model
generation jobs, since those cannot be continued and must instead be restarted from the beginning. The best way to do
that is to restart the job manually, outside the RELION GUI, and then when the job finishes RELION should continue as
if the job had never failed.

For example, with a failed local job:

    ps -e | grep relion                            # to check if the job is still active
    kill <process_id>                              # to stop the job
    # now re-run the commands from the job's `note.txt' file

or with a job that was submitted to an SGE cluster queue:

    qstat                                          # to check if the job is still active in the queue
    qdel <job_id>                                  # to remove the job from the queue
    qsub job_type/job_directory/run_submit.script  # to re-submit the job

The other option is to just run a new job from the RELION GUI in the normal way (select the job you want to "copy" in
the jobs list, make a "new" job by clicking on the job type in the list in the top-left of the GUI, then click
"Run!"). However, if you do this, relion_it.py will not know about the new job and will not run any further
downstream processing based on it. In this situation, you can either continue to process your data manually in RELION,
or you could edit the `RELION_IT_SUBMITTED_JOBS' file to replace the failed job with the manual one, and delete the
jobs that followed the original one. After that, if you re-run the script it should continue as normal from that
job onwards.


Repeated job failure
````````````````````
If a job fails repeatedly, it usually indicates that there is some problem with the job parameters or the files that
the job needs to access.

In favourable cases, it's possible you could fix the problem by selecting the job in the RELION GUI, changing one of
the parameters that is not greyed out, then clicking "Continue!". Often, though, the problem will be with one of the
parameters that can't be changed for a job that already exists, so the job will need to be deleted and recreated with
a different set of parameters.

To handle this situation, stop all of the pipelines and the relion_it.py script ("rm RUNNING_*"), then identify and
fix the problem. Often, the problem will be an error in one of the job parameters, which can usually be fixed by
changing one of the script options (for example by changing the settings in `relion_it_options.py', if you originally
used the GUI to start the run).

If the problem is caused by missing files from an upstream job, you might need to check the output of previous jobs
and look in the job directories to figure out what the problem is. Again, if it's an error in the parameters for a
job, you can probably fix it by editing `relion_it_options.py'.

After changing any script options, you'll need to use the RELION GUI to delete the affected job and all jobs
downstream of it, and also remove them from the list in the `RELION_IT_SUBMITTED_JOBS' file. Then you should be able
to restart the pipelines by running "relion_it.py --continue".

If you still can't get a particular job to run without errors, you can at least continue to run the upstream jobs
that are working properly. You can do this either by changing the options for relion_it.py (there are options to
switch off 2D or 3D classification, or to stop after CTF estimation), or by manually scheduling the jobs you want
using the RELION GUI. Remember that after running relion_it.py, you have a normal RELION project, so if the script
can't do what you want, you can simply stop it and then use all of RELION's normal job management and scheduling
abilities.


Advanced usage
--------------

It's possible to customise many aspects of the way relion_it.py works, but the details go beyond the scope of this
introduction. Simple customisation can be done by setting appropriate option values (see "Option files" above). For
more substantial changes, you might need to edit the script's Python code to get the behaviour you want. Most of the
important logic is in the `RelionItOptions.run_pipeline' function so that's a good place to start. Good luck!

'''

from __future__ import print_function
from __future__ import division  # always use float division

import argparse
import glob
import inspect
import math
import os
import runpy
import time
import traceback
import re
import star
import itertools

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

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING_RELION_IT'
SECONDPASS_REF3D_FILE = 'RELION_IT_2NDPASS_3DREF'
SETUP_CHECK_FILE = 'RELION_IT_SUBMITTED_JOBS'
PREPROCESS_SCHEDULE_PASS1 = 'PREPROCESS'
PREPROCESS_SCHEDULE_PASS2 = 'PREPROCESS_PASS2'
OPTIONS_FILE = 'relion_it_options.py'

# Generators save memory.
range = xrange


def is_dunder_name(name):
    return name.startswith('__') and name.endswith('__')


def prefix_RELION_IT(msg):
    return ' RELION_IT: ' + msg


def prefix_ERROR(msg):
    return ' ERROR: ' + msg


def bool_to_word(x):
    return 'Yes' if x else 'No'


class RelionItOptions(object):
    '''
    Options for the relion_it pipeline setup script.

    When initialised, this contains default values for all options. 
    Call `update_from` to override the defaults with a dictionary of new values.
    '''
    #############################################################################
    # Change the parameters below to reflect your experiment                    #
    # Current defaults reflect cryo-ARM betagal data set of RELION-3.0 tutorial #
    #############################################################################

    ### General parameters
    # Pixel size in Angstroms in the input movies
    angpix = 0.885
    # Acceleration voltage (in kV)
    voltage = 200
    # Polara = 2.0; Talos/Krios = 2.7; some Cryo-ARM = 1.4
    Cs = 1.4

    ### Import images (Linux wild card; movies as *.mrc, *.mrcs, *.tiff or *.tif; single-frame micrographs as *.mrc)
    import_images = 'Movies/*.tiff'
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    images_are_movies = True

    ### MotionCorrection parameters
    # Dose in electrons per squared Angstrom per fraction
    motioncor_doseperframe = 1.277
    # Gain-reference image in MRC format (only necessary if input movies are not yet gain-corrected, e.g. compressed TIFFs from K2)
    motioncor_gainreference = 'Movies/gain.mrc'
    # EER upsampling (1 = 4K, 2 = 8K). If you use 8K rendering, the pixel size (angpix) MUST be the half of the physical pixel size and the motioncor_binning should be 2.
    eer_upsampling = 1
    # EER fractionation. The dose rate (motioncor_doseperframe) is e/A2/fraction after this fractionation.
    eer_grouping = 20

    ### CTF estimation parameters
    # Most cases won't need changes here...

    ### Autopick parameters
    # Use reference-free Laplacian-of-Gaussian picking (otherwise use reference-based template matching instead)
    autopick_do_LoG = True
    # Minimum and maximum diameter in Angstrom for the LoG filter
    autopick_LoG_diam_min = 150
    autopick_LoG_diam_max = 180
    # Use positive values (0-1) to pick fewer particles; use negative values (-1-0) to pick more particles
    autopick_LoG_adjust_threshold = 0.0
    autopick_LoG_upper_threshold = 999.0
    #
    # OR:
    #
    # References for reference-based picking (when autopick_do_LoG = False)
    autopick_2dreferences = ''
    # OR: provide a 3D references for reference-based picking (when autopick_do_LoG = False)
    autopick_3dreference = ''

    # Threshold for reference-based autopicking (threshold 0 will pick too many particles. Default of 0.4 is hopefully better. Ultimately, just hope classification will sort it all out...)
    autopick_refs_threshold = 0.4
    # Minimum inter-particle distance for reference-based picking (~70% of particle diameter often works well)
    autopick_refs_min_distance = 120
    #
    # For both LoG and refs:
    #
    # Use this to remove false positives from carbon edges (useful range: 1.0-1.2, -1 to switch off)
    autopick_stddev_noise = -1
    # Use this to remove false positives from carbon edges (useful range: -0.5-0.0; -999 to switch off)
    autopick_avg_noise = -999

    ### Extract parameters
    # Box size of particles in the averaged micrographs (in pixels)
    extract_boxsize = 256
    # Down-scale the particles upon extraction?
    extract_downscale = False
    # Box size of the down-scaled particles (in pixels)
    extract_small_boxsize = 64
    # In second pass, down-scale the particles upon extraction?
    extract2_downscale = False
    # In second pass, box size of the down-scaled particles (in pixels)
    extract2_small_boxsize = 128

    ### Now perform 2D and/or 3D classification with the extracted particles?
    do_class2d = True
    # And/or perform 3D classification?
    do_class3d = True
    # Repeat 2D and/or 3D-classification for batches of this many particles
    batch_size = 10000
    # Number of 2D classes to use
    class2d_nr_classes  = 50
    # Diameter of the mask used for 2D/3D classification (in Angstrom)
    mask_diameter = 190
    # Symmetry group (when using SGD for initial model generation, C1 may work best)
    symmetry = 'C1'
    #
    ### 3D-classification parameters
    # Number of 3D classes to use
    class3d_nr_classes = 4
    # Have initial 3D model? If not, calculate one using SGD initial model generation
    have_3d_reference = False
    # Initial reference model
    class3d_reference = ''
    # Is reference on correct greyscale?
    class3d_ref_is_correct_greyscale = False
    # Has the initial reference been CTF-corrected?
    class3d_ref_is_ctf_corrected = True
    # Initial lowpass filter on reference
    class3d_ini_lowpass = 40

    ### Use the largest 3D class from the first batch as a 3D reference for a second pass of autopicking? (only when do_class3d is True)
    do_second_pass = True
    # Only move on to template-based autopicking if the 3D references achieves this resolution (in A)
    minimum_resolution_3dref_2ndpass = 20
    # In the second pass, perform 2D classification?
    do_class2d_pass2 = True
    # In the second pass, perform 3D classification?
    do_class3d_pass2 = False
    # Batch size in the second pass
    batch_size_pass2 = 100000

    ###################################################################################
    ############ Often the parameters below can be kept the same for a given set-up
    ###################################################################################

    ### Repeat settings for entire pipeline
    # Repeat the pre-processing runs this many times (or until RUNNING_PIPELINER_default_PREPROCESS file is deleted)
    preprocess_repeat_times = 999
    # Wait at least this many minutes between each repeat cycle
    preprocess_repeat_wait = 1
    ### Stop after CTF estimation? I.e., skip autopicking, extraction, 2D/3D classification, etc?
    stop_after_ctf_estimation = False
    # Check every this many minutes if enough particles have been extracted for a new batch of 2D-classification
    batch_repeat_time = 1

    ### MotionCorrection parameters
    # Use RELION's own implementation of motion-correction (CPU-only) instead of the UCSF implementation?
    motioncor_do_own = True
    # The number of threads (only for RELION's own implementation) is optimal when nr_movie_frames/nr_threads = integer
    motioncor_threads = 6
    # Exectutable of UCSF MotionCor2
    motioncor_exe = '/public/EM/MOTIONCOR2/MotionCor2'
    # On which GPU(s) to execute UCSF MotionCor2
    motioncor_gpu = '0'
    # How many MPI processes to use for running motion correction?
    motioncor_mpi = 4
    # Local motion-estimation patches for MotionCor2
    motioncor_patches_x = 4
    motioncor_patches_y = 4
    # B-factor in A^2 for downweighting of high-spatial frequencies
    motioncor_bfactor = 150
    # Use binning=2 for super-resolution movies
    motioncor_binning = 1
    # Provide a defect file for your camera if you have one
    motioncor_defectfile = ''
    # orientation of the gain-reference w.r.t your movies (if input movies are not yet gain-corrected, e.g. TIFFs)
    motioncor_gainflip = 'No flipping (0)'
    motioncor_gainrot = 'No rotation (0)'
    # Other arguments for MotionCor2
    motioncor_other_args = ''
    # Submit motion correction job to the cluster?
    motioncor_submit_to_queue = False

    ### CTF estimation parameters
    # Amplitude contrast (Q0)
    ampl_contrast = 0.1
    # CTFFIND-defined parameters
    ctffind_boxsize = 512
    ctffind_astigmatism = 100
    ctffind_maxres = 5
    ctffind_minres = 30
    ctffind_defocus_max = 50000
    ctffind_defocus_min = 5000
    ctffind_defocus_step = 500
    # For Gctf: ignore parameters on the 'Searches' tab?
    ctffind_do_ignore_search_params = True
    # For Gctf: perform equi-phase averaging?
    ctffind_do_EPA = True
    # Also estimate phase shifts (for VPP data)
    ctffind_do_phaseshift = False
    # Executable to Kai Zhang's Gctf
    gctf_exe = '/public/EM/Gctf/bin/Gctf'
    # On which GPU(s) to execute Gctf
    gctf_gpu = '0'
    # Use Alexis Rohou's CTFFIND4 (CPU-only) instead?
    use_ctffind = True
    # Executable for Alexis Rohou's CTFFIND4
    ctffind4_exe = '/public/EM/ctffind/ctffind.exe'
    # How many MPI processes to use for running CTF estimation?
    ctffind_mpi = 8
    # Submit CTF estimation job to the cluster?
    ctffind_submit_to_queue = False

    ### Autopick parameters
    # Use GPU-acceleration for autopicking?
    autopick_do_gpu = True
    # Which GPU(s) to use for autopicking
    autopick_gpu = '0'
    # Low-pass filter for auto-picking the micrographs
    autopick_lowpass = 20
    # Shrink factor for faster picking (0 = fastest; 1 = slowest)
    autopick_shrink_factor = 0
    # How many MPI processes to use for running auto-picking?
    autopick_mpi = 1
     # Additional arguments for autopicking
    autopick_other_args = ''
    # Submit Autopick job to the cluster?
    autopick_submit_to_queue = False
    # Are the references CTF-corrected?
    autopick_refs_are_ctf_corrected = True
    # Do the references have inverted contrast wrt the micrographs?
    autopick_refs_have_inverted_contrast = True
    # Ignore CTFs until the first peak
    autopick_refs_ignore_ctf1stpeak = False
    # Diameter of mask for the references (in A; negative value for automated detection of mask diameter)
    autopick_refs_mask_diam = -1
    # In-plane angular sampling interval
    autopick_inplane_sampling = 10
    # Symmetry of the 3D reference for autopicking
    autopick_3dref_symmetry = 'C1'
    # 3D angular sampling for generating projections of the 3D reference for autopicking (30 degrees is usually enough)
    autopick_3dref_sampling = '30 degrees'
    # Pixel size in the provided 2D/3D references (negative for same as in motion-corrected movies)
    autopick_ref_angpix = -1

    ### Extract parameters
    # Diameter for background normalisation (in pixels; negative value: default is 75% box size)
    extract_bg_diameter = -1
    # How many MPI processes to use for running particle extraction?
    extract_mpi = 1
    # Submit Extract job to the cluster?
    extract_submit_to_queue = False

    ## Discard particles based on average/stddev values? (this may be important for SGD initial model generation)
    do_discard_on_image_statistics = False
    # Discard images that have average/stddev values that are more than this many sigma away from the ensemble average
    discard_sigma = 4
    # Submit discard job to the cluster?
    discard_submit_to_queue = False

    #### Common relion_refine paremeters used for 2D/3D classification and initial model generation
    # Read all particles in one batch into memory?
    refine_preread_images = False
    # Or copy particles to scratch disk?
    refine_scratch_disk = ''
    # Number of pooled particles?
    refine_nr_pool = 10
    # Use GPU-acceleration?
    refine_do_gpu = True
    # Which GPU to use (different from GPU used for pre-processing?)
    refine_gpu = '1'
    # How many MPI processes to use
    refine_mpi = 1
    # How many threads to use
    refine_threads = 6
    # Skip padding?
    refine_skip_padding = False
    # Submit jobs to the cluster?
    refine_submit_to_queue = False
    # Use fast subsets in 2D/3D classification when batch_size is bigger than this
    refine_batchsize_for_fast_subsets = 10000

    ### 2D classification parameters
    # Wait with the first 2D classification batch until at least this many particles are extracted
    minimum_batch_size = 10000
    # Number of iterations to perform in 2D classification
    # Must be at least 20 for fast subsets
    class2d_nr_iter = 20
    # Rotational search step (in degrees)
    class2d_angle_step = 6
    # Offset search range (in pixels)
    class2d_offset_range = 5
    # Offset search step (in pixels)
    class2d_offset_step = 1
    # Option to ignore the CTFs until their first peak (try this if all particles go into very few classes)
    class2d_ctf_ign1stpeak = False
    # Additional arguments to pass to relion-refine
    class2d_other_args = ''

    ### 3D classification parameters
    # Number of iterations to perform in 3D classification
    # Must be at least 20 for fast subsets
    class3d_nr_iter = 20
    # Reference mask
    class3d_reference_mask = ''
    # Option to ignore the CTFs until their first peak (try this if all particles go into very few classes)
    class3d_ctf_ign1stpeak = False
    # Regularisation parameter (T)
    class3d_T_value = 4
    # Angular sampling step
    class3d_angle_step = '7.5 degrees'
    # Offset search range (in pixels)
    class3d_offset_range = 5
    # Offset search step (in pixels)
    class3d_offset_step = 1
    # Additional arguments to pass to relion-refine
    class3d_other_args = ''

    ## SGD initial model generation
    # Number of models to generate simulatenously (K>1 may be useful for getting rid of outliers in the particle images)
    inimodel_nr_classes = 4
    # Ignore CTFs until first peak?
    inimodel_ctf_ign1stpeak = False
    # Enforce non-negative solvent?
    inimodel_solvent_flatten = True
    # Initial angular sampling
    inimodel_angle_step = '15 degrees'
    # Initial search range (in pixels)
    inimodel_offset_range = 6
    # Initial offset search step (in pixels)
    inimodel_offset_step = 2
    # Number of initial iterations
    inimodel_nr_iter_initial = 50
    # Number of in-between iterations
    inimodel_nr_iter_inbetween = 200
    # Number of final iterations
    inimodel_nr_iter_final = 50
    # Frequency to write out information
    inimodel_freq_writeout = 10
    # Initial resolution (in A)
    inimodel_resol_ini = 35
    # Final resolution (in A)
    inimodel_resol_final = 15
    # Initial mini-batch size
    inimodel_batchsize_ini = 100
    # Final mini-batch size
    inimodel_batchsize_final = 500
    # Increased noise variance half-life (off, i.e. -1, by default; values of ~1000 have been observed to be useful in difficult cases)
    inimodel_sigmafudge_halflife = -1
    # Additional arguments to pass to relion_refine (skip annealing to get rid of outlier particles)
    inimodel_other_args = ' --sgd_skip_anneal '

    ### Cluster submission settings
    # Name of the queue to which to submit the job
    queue_name = 'openmpi'
    # Name of the command used to submit scripts to the queue
    queue_submit_command = 'qsub'
    # The template for your standard queue job submission script
    queue_submission_template = '/public/EM/RELION/relion/bin/qsub.csh'
    # Minimum number of dedicated cores that need to be requested on each node
    queue_minimum_dedicated = 1

    ### End of options

    #######################################################################
    ############ typically no need to change anything below this line
    #######################################################################

    def update_from(self, d):
        '''
        Update this RelionItOptions object from a dictionary.

        Special values (with names like '__xxx__') are removed, allowing this
        method to be given a dictionary containing the namespace from a script
        run with ``runpy``.
        '''
        for key, value in d.items():
            if not is_dunder_name(key):  # exclude __name__, __builtins__ etc.
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(prefix_RELION_IT("Unrecognised option '{}'".format(key)))

    def print_options(self, filename):
        '''
        Write the current options to `filename`.

        This method writes the options in the same format as they are read,
        allowing options to be written to a file and re-used.

        Args:
            out_file: A file object (optional). If supplied, options will be
                written to this file, otherwise they will be printed to
                sys.stdout.

        Raises:
            ValueError: If there is a problem printing the options.
        '''
        # NOTE The writing to stdout mentioned in the above docstring is not implemented.
        with open(filename, 'w') as out_file:
            out_file.write("# Options file for relion_it.py\n\n")
            option_names = [
                key for key in dir(self)
                if not is_dunder_name(key) and not callable(getattr(self, key))
            ]

            def relevant(lines):
                return itertools.takewhile(
                    lambda x: x != '### End of options',
                    itertools.dropwhile(
                        lambda x: x != '### General parameters', 
                        map(str.strip, lines)
                ))

            assignmentpattern = re.compile(r'(.*)=(.*)')

            # Parse the source code for this class,
            # and write out all comments along with option lines containing new values
            for line in relevant(inspect.getsourcelines(RelionItOptions)[0]):
                if line.startswith('#') or not line:
                    # Print comments and blank lines as-is
                    out_file.write(line + "\n")
                else:
                    # Assume all other lines define an option name and value.
                    # Replace with new value.
                    m = re.match(assignmentpattern, line)
                    if m:
                        option_name = m.group(1).strip()
                        if option_name in option_names:
                            out_file.write('{} = {}\n'.format(option_name, repr(getattr(self, option_name))))
                            option_names.remove(option_name)
                        else:
                            # This error should not occur.
                            # If it does, there is probably a programming error.
                            raise ValueError("Unrecognised option name '{}'".format(option_name))
            if option_names:
                # This error should not occur.
                # If it does, there is probably a programming error.
                raise ValueError("Some options were not written to the output file: {}".format(option_names))

    def run_pipeline(self):
        '''
        Configure and run the RELION 3 pipeline with the given options.
        '''
        # Is this really necessary?
        # Don't think so...
        if not os.path.isfile(PIPELINE_STAR):
            with open(PIPELINE_STAR, 'w') as writefile:
                for line in ('data_pipeline_general', '_rlnPipeLineJobCounter 1'):
                    writefile.write(line + '\n')

        # Write RUNNING_RELION_IT file
        # When deleted, this script will stop
        with open(RUNNING_FILE, 'w'):
            pass

        # Write main GUI project file, so GUI won't ask to set up a project
        with open('.gui_projectdir', 'w'):
            pass

        # Set up GUI file for Manualpick job to allow easy viewing of autopick results
        writeManualPickingGuiFile(
            self.autopick_LoG_diam_min if self.autopick_do_LoG else 
            self.autopick_refs_min_distance
        )

        ### Prepare the list of queue arguments for later use
        queue_options = [
            '{} == {}'.format(question, answer) for question, answer in [
                ('Submit to queue?',                  'Yes'),
                ('Queue name:',                       self.queue_name),
                ('Queue submit command:',             self.queue_submit_command),
                ('Standard submission script:',       self.queue_submission_template),
                ('Minimum dedicated cores per node:', self.queue_minimum_dedicated),
            ]
        ]

        # If we're only doing motioncorr and ctf estimation,
        # forget about the second pass and the batch processing
        if self.stop_after_ctf_estimation:
            self.do_class2d = False
            self.do_class3d = False
            self.do_second_pass = False

        nr_passes = 2 if self.do_second_pass else 1

        # If SECONDPASS_REF3D_FILE exists,
        # go straight into the second pass
        first_pass = 0
        if self.do_second_pass:
            secondpass_ref3d, secondpass_ref3d_angpix = getSecondPassReference()
            if secondpass_ref3d != '':
                for msg in [
                    'found {} with angpix= {} as a 3D reference for second pass in file {}'.format(
                        secondpass_ref3d, secondpass_ref3d_angpix, SECONDPASS_REF3D_FILE
                    ),
                    'if the automatic selection of the reference turned out to be unsatisfactory,',
                    'you can re-run the second pass with another reference by:',
                    ' stopping the pipeline by deleting RUNNING_*',
                    ' updating the reference filename in {}'.format(SECONDPASS_REF3D_FILE),
                    ' deleting relevant jobs (autopick2_job and followings) in {}'.format(SETUP_CHECK_FILE),
                    ' and restarting the pipeline.',
                ]:
                    print(prefix_RELION_IT(msg))
                first_pass = 1
                self.autopick_3dreference = secondpass_ref3d
                self.autopick_ref_angpix = secondpass_ref3d_angpix
                self.autopick_2dreferences = ''
                self.autopick_do_LoG = False
                self.class3d_reference = secondpass_ref3d
                self.have_3d_reference = True

        # Allow to perform two passes through the entire pipeline (PREPROCESS and CLASS2D/3D batches)
        # The second pass, a 3D reference generated in the first pass will be used for template-based autopicking
        for ipass in range(first_pass, nr_passes):

            #### Set up the Import job
            import_options = [
                '{} == {}'.format(question, answer) for question, answer in [
                    ('Raw input files:',               self.import_images),
                    ('Import raw movies/micrographs?', 'Yes'),
                    ('Pixel size (Angstrom):',         self.angpix),
                    ('Voltage (kV):',                  self.voltage),
                    ('Spherical aberration (mm):',     self.Cs),
                    ('Amplitude contrast:',            self.ampl_contrast),
                    ('Are these multi-frame movies?',  bool_to_word(self.images_are_movies)),
                ]
            ]

            import_job, already_had_it = addJob('Import','import_job', SETUP_CHECK_FILE, import_options)

            if self.images_are_movies:
                #### Set up the MotionCor job
                motioncorr_options = [
                    '{} == {}'.format(question, answer) for question, answer in [
                        ('Input movies STAR file:',    str(import_job) + 'movies.star'),
                        ('MOTIONCOR2 executable:',     self.motioncor_exe),
                        ('Defect file:',               self.motioncor_defectfile),
                        ('Gain-reference image:',      self.motioncor_gainreference),
                        ('Gain flip:',                 self.motioncor_gainflip),
                        ('Gain rotation:',             self.motioncor_gainrot),
                        ('Do dose-weighting?',         'Yes'),
                        ('Dose per frame (e/A2):',     self.motioncor_doseperframe),
                        ('Number of patches X:',       self.motioncor_patches_x),
                        ('Number of patches Y:',       self.motioncor_patches_y),
                        ('Bfactor:',                   self.motioncor_bfactor),
                        ('Binning factor:',            self.motioncor_binning),
                        ('Which GPUs to use:',         self.motioncor_gpu),
                        ('Other MOTIONCOR2 arguments', self.motioncor_other_args),
                        ('Number of threads:',         self.motioncor_threads),
                        ('Number of MPI procs:',       self.motioncor_mpi),
                        ('Additional arguments:',      ' '.join((
                            '--eer_upsampling', str(self.eer_upsampling),
                            '--eer_grouping', str(self.eer_grouping),
                        ))),
                        ('Use RELION\'s own implementation?', bool_to_word(self.motioncor_do_own)),
                    ]
                ]

                if self.motioncor_do_own:
                    motioncorr_options.append('Save sum of power spectra? == {}'.format(
                        bool_to_word(self.use_ctffind)
                    ))

                if self.motioncor_submit_to_queue:
                    motioncorr_options.extend(queue_options)

                motioncorr_job, already_had_it  = addJob(
                    'MotionCorr', 'motioncorr_job',
                    SETUP_CHECK_FILE, motioncorr_options,
                )

            # Set up the CtfFind job
            ctffind_options = [
                '{} == {}'.format(question, answer) for question, answer in [
                    ('Amount of astigmatism (A):', self.ctffind_astigmatism),
                    ('FFT box size (pix):',        self.ctffind_boxsize),
                    ('Maximum defocus value (A):', self.ctffind_defocus_max),
                    ('Minimum defocus value (A):', self.ctffind_defocus_min),
                    ('Defocus step size (A):',     self.ctffind_defocus_step),
                    ('Maximum resolution (A):',    self.ctffind_maxres),
                    ('Minimum resolution (A):',    self.ctffind_minres),
                    ('Gctf executable:',           self.gctf_exe),
                    ('Which GPUs to use:',         self.gctf_gpu),
                    ('CTFFIND-4.1 executable:',    self.ctffind4_exe),
                    ('Number of MPI procs:',       self.ctffind_mpi),
                    ('Input micrographs STAR file:', (
                        str(motioncorr_job) + 'corrected_micrographs.star' if self.images_are_movies else
                        str(import_job) + 'micrographs.star'
                    )),
                    ('Use CTFFIND-4.1?',           bool_to_word(self.use_ctffind)),
                    ('Use Gctf instead?',          bool_to_word(not self.use_ctffind)),
                    ('Use power spectra from MotionCorr job?', bool_to_word(self.use_ctffind)),
                ]
            ]

            if not self.use_ctffind:
                ctffind_options.append('Ignore \'Searches\' parameters? == {}'.format(
                    bool_to_word(self.ctffind_do_ignore_search_params)
                ))
                ctffind_options.append('Perform equi-phase averaging? == {}'.format(
                    bool_to_word(self.ctffind_do_EPA)
                ))
            ctffind_options.append('Estimate phase shifts? == {}'.format(
                bool_to_word(self.ctffind_do_phaseshift)
            ))

            if self.ctffind_submit_to_queue:
                ctffind_options.extend(queue_options)

            ctffind_job, already_had_it = addJob(
                'CtfFind', 'ctffind_job', SETUP_CHECK_FILE, ctffind_options
            )

            runjobs = [import_job]
            if self.images_are_movies:
                runjobs.append(motioncorr_job)
            runjobs.append(ctffind_job)

            do_2d_classification = (
                ipass == 0 and self.do_class2d
            ) or (
                ipass == 1 and self.do_class2d_pass2
            )

            do_3d_classification = (
                ipass == 0 and self.do_class3d
            ) or (
                ipass == 1 and self.do_class3d_pass2
            )

            do_classification = do_2d_classification or do_3d_classification

            downscale = (
                ipass == 0 and self.extract_downscale
            ) or (
                ipass == 1 and self.extract2_downscale
            )

            # There is an option to stop on-the-fly processing after CTF estimation
            if not self.stop_after_ctf_estimation:
                autopick_options = [
                    '{} == {}'.format(question, answer) for question, answer in [
                        ('Input micrographs for autopick:',      ctffind_job + 'micrographs_ctf.star'),
                        ('Min. diameter for LoG filter (A)',     self.autopick_LoG_diam_min),
                        ('Max. diameter for LoG filter (A)',     self.autopick_LoG_diam_max),
                        ('Maximum resolution to consider (A)',   self.autopick_lowpass),
                        ('Adjust default threshold (stddev):',   self.autopick_LoG_adjust_threshold),
                        ('Upper threshold (stddev):',            self.autopick_LoG_upper_threshold),
                        ('2D references:',                       self.autopick_2dreferences),
                        ('3D reference:',                        self.autopick_3dreference),
                        ('Symmetry:',                            self.autopick_3dref_symmetry),
                        ('Pixel size in references (A)',         self.autopick_ref_angpix),
                        ('3D angular sampling:',                 self.autopick_3dref_sampling),
                        ('In-plane angular sampling (deg)',      self.autopick_inplane_sampling),
                        ('Picking threshold:',                   self.autopick_refs_threshold),
                        ('Minimum inter-particle distance (A):', self.autopick_refs_min_distance),
                        ('Mask diameter (A)',                    self.autopick_refs_mask_diam),
                        ('Maximum stddev noise:',                self.autopick_stddev_noise),
                        ('Minimum avg noise:',                   self.autopick_avg_noise),
                        ('Shrink factor:',                       self.autopick_shrink_factor),
                        ('Which GPUs to use:',                   self.autopick_gpu),
                        ('Additional arguments:',                self.autopick_other_args),
                        ('Number of MPI procs:',                 self.autopick_mpi),
                        ('OR: provide a 3D reference?',          bool_to_word(self.autopick_3dreference != '')),
                        ('OR: use Laplacian-of-Gaussian?',       bool_to_word(self.autopick_do_LoG)),
                        ('Are References CTF corrected?',        bool_to_word(self.autopick_refs_are_ctf_corrected)),
                        ('References have inverted contrast?',   bool_to_word(self.autopick_refs_have_inverted_contrast)),
                        ('Ignore CTFs until first peak?',        bool_to_word(self.autopick_refs_ignore_ctf1stpeak)),
                        ('Use GPU acceleration?',                bool_to_word(self.autopick_do_gpu and not self.autopick_do_LoG)),
                    ]
                ]

                if self.autopick_submit_to_queue:
                    autopick_options.extend(queue_options)

                autopick_job_name, autopick_alias = (
                    ('autopick_job',  'pass 1') if ipass == 0 else
                    ('autopick2_job', 'pass 2')
                )

                autopick_job, already_had_it = addJob(
                    'AutoPick', autopick_job_name, SETUP_CHECK_FILE,
                    autopick_options, alias=autopick_alias,
                )
                runjobs.append(autopick_job)

                # Extract options
                extract_options = [
                    '{} == {}'.format(question, answer) for question, answer in [
                        ('Input coordinates:',                autopick_job + 'coords_suffix_autopick.star'),
                        ('micrograph STAR file:',             ctffind_job + 'micrographs_ctf.star'),
                        ('Diameter background circle (pix):', self.extract_bg_diameter),
                        ('Particle box size (pix):',          self.extract_boxsize),
                        ('Number of MPI procs:',              self.extract_mpi),
                    ]
                ]

                if downscale:
                    extract_options.extend([
                        '{} == {}'.format(question, answer) for question, answer in [
                            ('Rescale particles?', 'Yes'),
                            ('Re-scaled size (pixels):', (
                                self.extract_small_boxsize if ipass == 0 else 
                                self.extract2_small_boxsize
                            )),
                        ]
                    ])

                if self.extract_submit_to_queue:
                    extract_options.extend(queue_options)

                extract_job_name, extract_alias = (
                    ('extract_job',  'pass 1') if ipass == 0 else
                    ('extract2_job', 'pass 2')
                )

                extract_job, already_had_it  = addJob(
                    'Extract', extract_job_name, SETUP_CHECK_FILE,
                    extract_options, alias=extract_alias,
                )
                runjobs.append(extract_job)

                if do_classification:
                    #### Set up the Select job to split the particle STAR file into batches
                    split_options = [
                        '{} == {}'.format(question, answer) for question, answer in [
                            ('OR select from particles.star:', extract_job + 'particles.star'),
                            ('OR: split into subsets?',        'Yes'),
                            ('OR: number of subsets:',         '-1'),
                            ('Subset size:', (
                                self.batch_size if ipass == 0 else 
                                self.batch_size_pass2
                            ))
                        ]
                    ]

                    split_job_name, split_alias = (
                        ('split_job',  'into {}'.format(self.batch_size)) if ipass == 0 else
                        ('split2_job', 'into {}'.format(self.batch_size_pass2))
                    )

                    split_job, already_had_it = addJob(
                        'Select', split_job_name, SETUP_CHECK_FILE,
                        split_options, alias=split_alias
                    )

                    # Now start running stuff
                    runjobs.append(split_job)

            # Now execute the entire preprocessing pipeliner
            preprocess_schedule_name = PREPROCESS_SCHEDULE_PASS1 if ipass == 0 else PREPROCESS_SCHEDULE_PASS2
            RunJobs(runjobs, self.preprocess_repeat_times, self.preprocess_repeat_wait, preprocess_schedule_name)
            print(prefix_RELION_IT('submitted {} pipeliner with {} repeats of the preprocessing jobs'.format(
                preprocess_schedule_name, self.preprocess_repeat_times
            )))
            print(prefix_RELION_IT(' '.join((
                'this pipeliner will run in the background of your shell.',
                'You can stop it by deleting the file RUNNING_PIPELINER_{}'.format(preprocess_schedule_name)
            ))))

            # From now on, process extracted particles in batches for 2D or 3D classification,
            # only perform SGD inimodel for first batch and if no 3D reference is available

            # There is again an option to stop here...
            if do_classification:
                # If necessary, rescale the 3D reference in the second pass!
                # TODO: rescale initial reference if different from movies?
                if ipass == 1 and (self.extract_downscale or self.extract2_downscale):
                    particles_angpix = self.angpix
                    if self.images_are_movies:
                        particles_angpix *= self.motioncor_binning
                    if self.extract2_downscale:
                        particles_angpix *= self.extract_boxsize / self.extract2_small_boxsize
                        particles_boxsize = self.extract2_small_boxsize
                    else:
                        particles_boxsize = self.extract_boxsize
                    if abs(float(particles_angpix) - float(self.autopick_ref_angpix)) > 0.01:
                        # Rescale the reference for 3D classification
                        print(prefix_RELION_IT('rescaling the 3D reference from pixel size {} to {} and saving the new reference as {}'.format(
                            self.autopick_ref_angpix, particles_angpix, self.class3d_reference
                        )))
                        self.class3d_reference = self.autopick_3dreference.replace('.mrc', '_rescaled.mrc')
                        command = (
                            'relion_image_handler'
                            + ' --i ' + str(self.autopick_3dreference)
                            + ' --o ' + str(self.class3d_reference)
                            + ' --angpix ' + str(self.autopick_ref_angpix)
                            + ' --rescale_angpix ' + str(particles_angpix)
                            + ' --new_box ' + str(particles_boxsize)
                        )
                        os.system(command)

                print(prefix_RELION_IT(' '.join((
                    'now entering an infinite loop for batch-processing of particles.',
                    'You can stop this loop by deleting the file {}'.format(RUNNING_FILE)
                ))))

                # It could be that this is a restart, so check previous_batch1_size in the output directory.
                # Also check the presence of class2d_job_batch_001 in case the first job was not submitted yet.
                first_split_file = find_split_job_output(split_job + 'particles_split', 1)
                if not any(x is None for x in (
                    getJobName("class2d_job_batch_001", SETUP_CHECK_FILE), first_split_file
                )):
                    batch1 = safe_load_star(first_split_file, expected=['particles', 'rlnMicrographName'])
                    previous_batch1_size = len(batch1['particles']['rlnMicrographName'])
                else:
                    previous_batch1_size = 0

                continue_this_pass = True
                while continue_this_pass:
                    have_new_batch = False
                    nr_batches = len(glob.glob(split_job + "particles_split*.star"))
                    for ibatch in range(nr_batches):
                        iibatch = ibatch + 1
                        batch_name = find_split_job_output(split_job + "particles_split", iibatch)

                        batch = safe_load_star(batch_name, expected=['particles', 'rlnMicrographName'])
                        batch_size = len(batch['particles']['rlnMicrographName'])
                        rerun_batch1 = False
                        if ibatch == 0 and batch_size > previous_batch1_size and batch_size > self.minimum_batch_size:
                            previous_batch1_size = batch_size
                            rerun_batch1 = True

                        particles_star_file = batch_name

                        # The first batch is special: 
                        # perform 2D classification with smaller batch size 
                        # (but at least minimum_batch_size) 
                        # and keep overwriting in the same output directory
                        if rerun_batch1 or batch_size == self.batch_size:

                            # Discard particles with odd average/stddev values
                            if self.do_discard_on_image_statistics:

                                # Run a Select job to get rid of particles with outlier average/stddev values...
                                discard_options = [
                                    '{} == {}'.format(question, answer) for question, answer in [
                                        ('OR select from particles.star:',     batch_name),
                                        ('OR: select on image statistics?',    'Yes'),
                                        ('Sigma-value for discarding images:', self.discard_sigma),
                                        ('Metadata label for images:',         'rlnImageName'),
                                    ]
                                ]

                                discard_job_name = 'discard_job' if ipass == 0 else 'discard2_job'

                                if self.discard_submit_to_queue:
                                    discard_options.extend(queue_options)

                                discard_job, already_had_it = addJob(
                                    'Select', discard_job_name,
                                    SETUP_CHECK_FILE, discard_options,
                                )

                                if rerun_batch1 or not already_had_it:
                                    have_new_batch = True
                                    RunJobs([discard_job], 1, 1, 'DISCARD')
                                    print(prefix_RELION_IT("submitted job to discard based on image statistics for {} particles in {}".format(
                                        batch_size, batch_name
                                    )))

                                    # Wait until Discard job is finished.
                                    # Check every thirty seconds.
                                    WaitForJob(discard_job, 30)

                                particles_star_file = discard_job + 'particles.star'

                            # 2D classification
                            if do_2d_classification:
                                class2d_options = [
                                    '{} == {}'.format(question, answer) for question, answer in [
                                        ('Input images STAR file:',                 particles_star_file),
                                        ('Number of classes:',                      self.class2d_nr_classes),
                                        ('Mask diameter (A):',                      self.mask_diameter),
                                        ('Number of iterations:',                   self.class2d_nr_iter),
                                        ('Angular search range - psi (deg):',       self.class2d_angle_step),
                                        ('Offset search range (pix):',              self.class2d_offset_range),
                                        ('Offset search step (pix):',               self.class2d_offset_step),
                                        ('Number of pooled particles:',             self.refine_nr_pool),
                                        ('Which GPUs to use:',                      self.refine_gpu),
                                        ('Number of MPI procs:',                    self.refine_mpi),
                                        ('Number of threads:',                      self.refine_threads),
                                        ('Copy particles to scratch directory:',    self.refine_scratch_disk),
                                        ('Additional arguments:',                   self.class2d_other_args),
                                        ('Use fast subsets (for large data sets)?', bool_to_word(batch_size > self.refine_batchsize_for_fast_subsets)),
                                        ('Use GPU acceleration?',                   bool_to_word(self.refine_do_gpu)),
                                        ('Ignore CTFs until first peak?',           bool_to_word(self.class2d_ctf_ign1stpeak)),
                                        ('Pre-read all particles into RAM?',        bool_to_word(self.refine_preread_images)),
                                    ]
                                ]

                                if self.refine_submit_to_queue:
                                    class2d_options.extend(queue_options)

                                jobname, alias = (
                                    ('class2d_job_batch_{:03d}'.format(iibatch),       'pass1_batch_{:03d}'.format(iibatch)) if ipass == 0 else
                                    ('class2d_pass2_job_batch_{:03d}'.format(iibatch), 'pass2_batch_{:03d}'.format(iibatch))
                                )

                                class2d_job, already_had_it = addJob(
                                    'Class2D', jobname, SETUP_CHECK_FILE,
                                    class2d_options, alias=alias
                                )

                                if rerun_batch1 or not already_had_it:
                                    have_new_batch = True
                                    RunJobs([class2d_job], 1, 1, 'CLASS2D')
                                    print(prefix_RELION_IT("submitted 2D classification with {} particles in {}".format(batch_size, class2d_job)))

                                    # Wait until Class2D job is finished.
                                    # Check every thirty seconds.
                                    WaitForJob(class2d_job, 30)

                        # Perform 3D classification
                        if do_3d_classification:
                            # Do SGD initial model generation only in the first pass,
                            # when no reference is provided AND only for the first (complete) batch, for subsequent batches use that model
                            if not self.have_3d_reference and ipass == ibatch == 0 and batch_size == self.batch_size:

                                inimodel_options = [
                                    '{} == {}'.format(question, answer) for question, answer in [
                                        ('Input images STAR file:',                   particles_star_file),
                                        ('Symmetry:',                                 self.symmetry),
                                        ('Mask diameter (A):',                        self.mask_diameter),
                                        ('Number of classes:',                        self.inimodel_nr_classes),
                                        ('Initial angular sampling:',                 self.inimodel_angle_step),
                                        ('Offset search range (pix):',                self.inimodel_offset_range),
                                        ('Offset search step (pix):',                 self.inimodel_offset_step),
                                        ('Number of initial iterations:',             self.inimodel_nr_iter_initial),
                                        ('Number of in-between iterations:',          self.inimodel_nr_iter_inbetween),
                                        ('Number of final iterations:',               self.inimodel_nr_iter_final),
                                        ('Write-out frequency (iter):',               self.inimodel_freq_writeout),
                                        ('Initial resolution (A):',                   self.inimodel_resol_ini),
                                        ('Final resolution (A):',                     self.inimodel_resol_final),
                                        ('Initial mini-batch size:',                  self.inimodel_batchsize_ini),
                                        ('Final mini-batch size:',                    self.inimodel_batchsize_final),
                                        ('Increased noise variance half-life:',       self.inimodel_sigmafudge_halflife),
                                        ('Number of pooled particles:',               '1'),
                                        ('Which GPUs to use:',                        self.refine_gpu),
                                        ('Number of MPI procs:',                      self.refine_mpi),
                                        ('Number of threads:',                        self.refine_threads),
                                        ('Copy particles to scratch directory:',      self.refine_scratch_disk),
                                        ('Additional arguments:',                     self.inimodel_other_args),
                                        ('Flatten and enforce non-negative solvent?', bool_to_word(self.inimodel_solvent_flatten)),
                                        ('Skip padding?',                             bool_to_word(self.refine_skip_padding)),
                                        ('Use GPU acceleration?',                     bool_to_word(self.refine_do_gpu)),
                                        ('Ignore CTFs until first peak?',             bool_to_word(self.inimodel_ctf_ign1stpeak)),
                                        ('Pre-read all particles into RAM?',          bool_to_word(self.refine_preread_images)),
                                    ]
                                ]

                                if self.refine_submit_to_queue:
                                    inimodel_options.extend(queue_options)

                                inimodel_job, already_had_it = addJob(
                                    'InitialModel', 'inimodel',
                                    SETUP_CHECK_FILE, inimodel_options,
                                )

                                if not already_had_it:
                                    have_new_batch = True
                                    RunJobs([inimodel_job], 1, 1, 'INIMODEL')
                                    print(prefix_RELION_IT("submitted initial model generation with {} particles in {}".format(batch_size, inimodel_job)))

                                    # Wait until inimodel job is finished.
                                    # Check every thirty seconds.
                                    WaitForJob(inimodel_job, 30)

                                sgd_model_star = findOutputModelStar(inimodel_job)
                                if sgd_model_star is None:
                                    print(prefix_RELION_IT("Initial model generation {} does not contain expected output maps.".format(inimodel_job)))
                                    print(prefix_RELION_IT("This job should have finished, but you may continue it from the GUI."))
                                    raise Exception("ERROR!! quitting the pipeline.") # TODO: MAKE MORE ROBUST

                                # Use the model of the largest class for the 3D classification below
                                total_iter = self.inimodel_nr_iter_initial + self.inimodel_nr_iter_inbetween + self.inimodel_nr_iter_final
                                best_inimodel_class, best_inimodel_resol, best_inimodel_angpix = findBestClass(sgd_model_star, use_resol=True)
                                self.class3d_reference = best_inimodel_class
                                self.class3d_ref_is_correct_greyscale = True
                                self.class3d_ref_is_ctf_corrected = True
                                self.have_3d_reference = True


                            if self.have_3d_reference:
                                # Now perform the actual 3D classification
                                class3d_options = [
                                    '{} == {}'.format(question, answer) for question, answer in [
                                        ('Input images STAR file:',                 particles_star_file),
                                        ('Reference map:',                          self.class3d_reference),
                                        ('Initial low-pass filter (A):',            self.class3d_ini_lowpass),
                                        ('Symmetry:',                               self.symmetry),
                                        ('Regularisation parameter T:',             self.class3d_T_value),
                                        ('Reference mask (optional):',              self.class3d_reference_mask),
                                        ('Number of classes:',                      self.class3d_nr_classes),
                                        ('Mask diameter (A):',                      self.mask_diameter),
                                        ('Number of iterations:',                   self.class3d_nr_iter),
                                        ('Angular sampling interval:',              self.class3d_angle_step),
                                        ('Offset search range (pix):',              self.class3d_offset_range),
                                        ('Offset search step (pix):',               self.class3d_offset_step),
                                        ('Number of pooled particles:',             self.refine_nr_pool),
                                        ('Which GPUs to use:',                      self.refine_gpu),
                                        ('Number of MPI procs:',                    self.refine_mpi),
                                        ('Number of threads:',                      self.refine_threads),
                                        ('Copy particles to scratch directory:',    self.refine_scratch_disk),
                                        ('Additional arguments:',                   self.class3d_other_args),
                                        ('Use fast subsets (for large data sets)?', bool_to_word(batch_size > self.refine_batchsize_for_fast_subsets)),
                                        ('Ref. map is on absolute greyscale?',      bool_to_word(self.class3d_ref_is_correct_greyscale)),
                                        ('Has reference been CTF-corrected?',       bool_to_word(self.class3d_ref_is_ctf_corrected)),
                                        ('Skip padding?',                           bool_to_word(self.refine_skip_padding)),
                                        ('Use GPU acceleration?',                   bool_to_word(self.refine_do_gpu)),
                                        ('Ignore CTFs until first peak?',           bool_to_word(self.class3d_ctf_ign1stpeak)),
                                        ('Pre-read all particles into RAM?',        bool_to_word(self.refine_preread_images)),
                                    ]
                                ]

                                if self.refine_submit_to_queue:
                                    class3d_options.extend(queue_options)

                                jobname, alias = (
                                    ('class3d_job_batch_{:03d}'.format(iibatch),  'pass1_batch_{:03d}'.format(iibatch)) if ipass == 0 else
                                    ('class3d2_job_batch_{:03d}'.format(iibatch), 'pass2_batch_{:03d}'.format(iibatch))
                                )

                                class3d_job, already_had_it = addJob(
                                    'Class3D', jobname, SETUP_CHECK_FILE,
                                    class3d_options, alias=alias,
                                )

                                if rerun_batch1 or not already_had_it:
                                    have_new_batch = True
                                    RunJobs([class3d_job], 1, 1, 'CLASS3D')
                                    print(prefix_RELION_IT('submitted 3D classification with {} particles in {}'.format(batch_size, class3d_job)))

                                    # Wait until Class2D job is finished.
                                    # Check every thirty seconds.
                                    WaitForJob(class3d_job, 30)

                                class3d_model_star = findOutputModelStar(class3d_job)
                                if class3d_model_star is None:
                                    print(prefix_RELION_IT("3D Classification {} does not contain expected output maps.".format(class3d_job)))
                                    print(prefix_RELION_IT("This job should have finished, but you may continue it from the GUI."))
                                    raise Exception("ERROR!! quitting the pipeline.") # TODO: MAKE MORE ROBUST

                                best_class3d_class, best_class3d_resol, best_class3d_angpix = findBestClass(class3d_model_star, use_resol=True)

                                # Once the first batch in the first pass is completed,
                                # move on to the second pass
                                if ipass == ibatch == 0 and self.do_second_pass and best_class3d_resol < self.minimum_resolution_3dref_2ndpass:
                                    self.autopick_3dreference = best_class3d_class
                                    self.autopick_ref_angpix = best_class3d_angpix
                                    self.autopick_2dreferences = ''
                                    self.autopick_do_LoG = False
                                    self.class3d_reference = best_class3d_class
                                    self.have_3d_reference = True
                                    self.autopick_3dref_symmetry = self.symmetry

                                    # Stop the PREPROCESS pipeliner of the first pass by removing its RUNNING file
                                    filename_to_remove = 'RUNNING_PIPELINER_' + preprocess_schedule_name
                                    if os.path.isfile(filename_to_remove):
                                        print(prefix_RELION_IT('removing file {} to stop the pipeliner from the first pass'.format(filename_to_remove)))
                                        os.remove(filename_to_remove)

                                    # Generate a file to indicate we're in the second pass,
                                    # so that restarts of the python script will be smooth
                                    with open(SECONDPASS_REF3D_FILE, 'w') as writefile:
                                        writefile.write(''.join(str(x) + '\n' for x in (
                                            best_class3d_class, best_class3d_angpix
                                        )))

                                    # Move out of this ipass of the passes loop....
                                    ibatch = nr_batches + 1
                                    continue_this_pass = False
                                    print(prefix_RELION_IT('moving on to the second pass using {} for template-based autopicking'.format(self.autopick_3dreference)))
                                    # break out of the for-loop over the batches
                                    break

                    if not have_new_batch:
                        CheckForExit()
                        # Don't check the particles.star file too often
                        # This will raise a NameError,
                        # since `self.batch_repeat_time` is not defined.
                        time.sleep(60 * self.batch_repeat_time)


class RelionItGui(object):

    def __init__(self, main_window, options):
        self.main_window = main_window
        self.options = options
        self.warnings = []

        # Convenience function for making file browser buttons
        def new_browse_button(self, var_to_set, filetypes=(('MRC file', '*.mrc'), ('All files', '*'))):
            def browse_command():
                chosen_file = tkFileDialog.askopenfilename(filetypes=filetypes)
                if chosen_file is not None:
                    # Make path relative if it's in the current directory
                    if chosen_file.startswith(os.getcwd()):
                        chosen_file = os.path.relpath(chosen_file)
                    var_to_set.set(chosen_file)
            return tk.Button(self, text="Browse...", command=browse_command)

        ### Create GUI

        main_frame = tk.Frame(main_window)
        main_frame.pack(fill=tk.BOTH, expand=1)

        left_frame, right_frame = 2 * (tk.Frame(main_frame),)
        for frame, side in zip((left_frame, right_frame), (tk.LEFT, tk.RIGHT)):
            frame.pack(side=side, anchor=tk.N, fill=tk.X, expand=1)

        ###

        expt_frame = tk.LabelFrame(left_frame, text="Experimental details", padx=5, pady=5)
        expt_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        tk.Label(expt_frame, text="Voltage (kV):").grid(row=0, sticky=tk.W)
        self.voltage_entry = tk.Entry(expt_frame)
        self.voltage_entry.grid(row=0, column=1, sticky=tk.W+tk.E)
        self.voltage_entry.insert(0, str(options.voltage))

        tk.Label(expt_frame, text="Cs (mm):").grid(row=1, sticky=tk.W)
        self.cs_entry = tk.Entry(expt_frame)
        self.cs_entry.grid(row=1, column=1, sticky=tk.W+tk.E)
        self.cs_entry.insert(0, str(options.Cs))

        tk.Label(expt_frame, text="Phase plate?").grid(row=2, sticky=tk.W)
        self.phaseplate_var = tk.BooleanVar()
        phaseplate_button = tk.Checkbutton(expt_frame, var=self.phaseplate_var)
        phaseplate_button.grid(row=2, column=1, sticky=tk.W)
        if options.ctffind_do_phaseshift:
            phaseplate_button.select()

        tk.Label(expt_frame, text=u"Pixel size (\u212B):").grid(row=3, sticky=tk.W)
        self.angpix_var = tk.StringVar()  # for data binding
        self.angpix_entry = tk.Entry(expt_frame, textvariable=self.angpix_var)
        self.angpix_entry.grid(row=3, column=1, sticky=tk.W+tk.E)
        self.angpix_entry.insert(0, str(options.angpix))

        tk.Label(expt_frame, text=u"Exposure rate (e\u207B / \u212B\u00B2 / frame):").grid(row=4, sticky=tk.W)
        self.exposure_entry = tk.Entry(expt_frame)
        self.exposure_entry.grid(row=4, column=1, sticky=tk.W + tk.E)
        self.exposure_entry.insert(0, str(options.motioncor_doseperframe))

        ###

        particle_frame = tk.LabelFrame(left_frame, text="Particle details", padx=5, pady=5)
        particle_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(particle_frame, 1, weight=1)

        tk.Label(particle_frame, text=u"Longest diameter (\u212B):").grid(row=0, sticky=tk.W)
        self.particle_max_diam_var = tk.StringVar()  # for data binding
        self.particle_max_diam_entry = tk.Entry(particle_frame, textvariable=self.particle_max_diam_var)
        self.particle_max_diam_entry.grid(row=0, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_max_diam_entry.insert(0, str(options.autopick_LoG_diam_max))

        tk.Label(particle_frame, text=u"Shortest diameter (\u212B):").grid(row=1, sticky=tk.W)
        self.particle_min_diam_entry = tk.Entry(particle_frame)
        self.particle_min_diam_entry.grid(row=1, column=1, sticky=tk.W+tk.E, columnspan=2)
        self.particle_min_diam_entry.insert(0, str(options.autopick_LoG_diam_min))

        tk.Label(particle_frame, text="3D reference (optional):").grid(row=2, sticky=tk.W)
        self.ref_3d_var = tk.StringVar()  # for data binding
        self.ref_3d_entry = tk.Entry(particle_frame, textvariable=self.ref_3d_var)
        self.ref_3d_entry.grid(row=2, column=1, sticky=tk.W+tk.E)
        self.ref_3d_entry.insert(0, str(options.autopick_3dreference))

        new_browse_button(particle_frame, self.ref_3d_var).grid(row=2, column=2)

        tk.Label(particle_frame, text=u"Mask diameter (\u212B):").grid(row=3, sticky=tk.W)
        self.mask_diameter_var = tk.StringVar()  # for data binding
        self.mask_diameter_entry = tk.Entry(particle_frame, textvariable=self.mask_diameter_var)
        self.mask_diameter_entry.grid(row=3, column=1, sticky=tk.W+tk.E)
        self.mask_diameter_entry.insert(0, str(options.mask_diameter))
        self.mask_diameter_px = tk.Label(particle_frame, text="= NNN px")
        self.mask_diameter_px.grid(row=3, column=2, sticky=tk.W)

        tk.Label(particle_frame, text="Box size (px):").grid(row=4, sticky=tk.W)
        self.box_size_var = tk.StringVar()  # for data binding
        self.box_size_entry = tk.Entry(particle_frame, textvariable=self.box_size_var)
        self.box_size_entry.grid(row=4, column=1, sticky=tk.W+tk.E)
        self.box_size_entry.insert(0, str(options.extract_boxsize))
        self.box_size_in_angstrom = tk.Label(particle_frame, text=u"= NNN \u212B")
        self.box_size_in_angstrom.grid(row=4, column=2, sticky=tk.W)

        tk.Label(particle_frame, text="Down-sample to (px):").grid(row=5, sticky=tk.W)
        self.extract_small_boxsize_var = tk.StringVar()  # for data binding
        self.extract_small_boxsize_entry = tk.Entry(particle_frame, textvariable=self.extract_small_boxsize_var)
        self.extract_small_boxsize_entry.grid(row=5, column=1, sticky=tk.W+tk.E)
        self.extract_small_boxsize_entry.insert(0, str(options.extract_small_boxsize))
        self.extract_angpix = tk.Label(particle_frame, text=u"= NNN \u212B/px")
        self.extract_angpix.grid(row=5, column=2, sticky=tk.W)

        tk.Label(particle_frame, text="Calculate for me:").grid(row=6, sticky=tk.W)
        self.auto_boxsize_var = tk.BooleanVar()
        self.auto_boxsize_button = tk.Checkbutton(particle_frame, var=self.auto_boxsize_var)
        self.auto_boxsize_button.grid(row=6, column=1, sticky=tk.W)
        self.auto_boxsize_button.select()

        ###

        project_frame = tk.LabelFrame(right_frame, text="Project details", padx=5, pady=5)
        project_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(project_frame, 1, weight=1)

        tk.Label(project_frame, text="Project directory:").grid(row=0, sticky=tk.W)
        tk.Label(project_frame, text=os.getcwd(), anchor=tk.W).grid(row=0, column=1, sticky=tk.W, columnspan=2)

        tk.Label(project_frame, text="Pattern for movies:").grid(row=1, sticky=tk.W)
        self.import_images_var = tk.StringVar()  # for data binding
        self.import_images_entry = tk.Entry(project_frame, textvariable=self.import_images_var)
        self.import_images_entry.grid(row=1, column=1, sticky=tk.W+tk.E)
        self.import_images_entry.insert(0, self.options.import_images)

        import_button = new_browse_button(
            project_frame, self.import_images_var,
            filetypes=(('Image file', '{*.mrc, *.mrcs, *.tif, *.tiff}'), ('All files', '*'))
        )
        import_button.grid(row=1, column=2)

        tk.Label(project_frame, text="Gain reference (optional):").grid(row=2, sticky=tk.W)
        self.gainref_var = tk.StringVar()  # for data binding
        self.gainref_entry = tk.Entry(project_frame, textvariable=self.gainref_var)
        self.gainref_entry.grid(row=2, column=1, sticky=tk.W+tk.E)
        self.gainref_entry.insert(0, self.options.motioncor_gainreference)

        new_browse_button(project_frame, self.gainref_var).grid(row=2, column=2)

        ###

        pipeline_frame = tk.LabelFrame(right_frame, text="Pipeline control", padx=5, pady=5)
        pipeline_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)
        tk.Grid.columnconfigure(expt_frame, 1, weight=1)

        tk.Label(pipeline_frame, text="Stop after CTF estimation?").grid(row=0, sticky=tk.W)
        self.stop_after_ctf_var = tk.BooleanVar()
        stop_after_ctf_button = tk.Checkbutton(pipeline_frame, var=self.stop_after_ctf_var)
        stop_after_ctf_button.grid(row=0, column=1, sticky=tk.W)
        if options.stop_after_ctf_estimation:
            stop_after_ctf_button.select()

        tk.Label(pipeline_frame, text="Do 2D classification?").grid(row=1, sticky=tk.W)
        self.class2d_var = tk.BooleanVar()
        self.class2d_button = tk.Checkbutton(pipeline_frame, var=self.class2d_var)
        self.class2d_button.grid(row=1, column=1, sticky=tk.W)
        if options.do_class2d:
            self.class2d_button.select()

        tk.Label(pipeline_frame, text="Do 3D classification?").grid(row=2, sticky=tk.W)
        self.class3d_var = tk.BooleanVar()
        self.class3d_button = tk.Checkbutton(pipeline_frame, var=self.class3d_var)
        self.class3d_button.grid(row=2, column=1, sticky=tk.W)
        if options.do_class3d:
            self.class3d_button.select()

        tk.Label(pipeline_frame, text="Do second pass? (only if no 3D ref)").grid(row=3, sticky=tk.W)
        self.second_pass_var = tk.BooleanVar()
        self.second_pass_button = tk.Checkbutton(pipeline_frame, var=self.second_pass_var)
        self.second_pass_button.grid(row=3, column=1, sticky=tk.W)
        if options.do_second_pass:
            self.second_pass_button.select()

        tk.Label(pipeline_frame, text="Do 2D classification (2nd pass)?").grid(row=4, sticky=tk.W)
        self.class2d_pass2_var = tk.BooleanVar()
        self.class2d_pass2_button = tk.Checkbutton(pipeline_frame, var=self.class2d_pass2_var)
        self.class2d_pass2_button.grid(row=4, column=1, sticky=tk.W)
        self.class2d_pass2_button.select()
        if options.do_class2d_pass2:
            self.class2d_pass2_button.select()

        tk.Label(pipeline_frame, text="Do 3D classification (2nd pass)?").grid(row=5, sticky=tk.W)
        self.class3d_pass2_var = tk.BooleanVar()
        self.class3d_pass2_button = tk.Checkbutton(pipeline_frame, var=self.class3d_pass2_var)
        self.class3d_pass2_button.grid(row=5, column=1, sticky=tk.W)
        if options.do_class3d_pass2:
            self.class3d_pass2_button.select()

        ### Add logic to the box size boxes

        self.box_size_var.trace('w', self.update_box_size_labels)
        self.extract_small_boxsize_var.trace('w', self.update_box_size_labels)

        self.angpix_var.trace('w', self.update_box_sizes)
        self.particle_max_diam_var.trace('w', self.update_box_sizes)
        self.auto_boxsize_button.config(command=self.update_box_sizes)

        ### Add logic to the check boxes

        stop_after_ctf_button.config(command=self.update_pipeline_control_state)
        self.class3d_button.config(command=self.update_pipeline_control_state)
        self.second_pass_button.config(command=self.update_pipeline_control_state)
        self.ref_3d_var.trace('w', self.update_pipeline_control_state)

        ###

        button_frame = tk.Frame(right_frame)
        button_frame.pack(padx=5, pady=5, fill=tk.X, expand=1)

        self.run_button = tk.Button(button_frame, text="Save & run", command=self.run_pipeline)
        self.run_button.pack(padx=5, pady=5, side=tk.RIGHT)

        self.save_button = tk.Button(button_frame, text="Save options", command=self.save_options)
        self.save_button.pack(padx=5, pady=5, side=tk.RIGHT)

        # Show initial pixel sizes
        self.update_box_sizes()

    @staticmethod
    def calculate_box_size(particlesize):
        '''Return (even) side length 20% larger than particle.'''
        boxsize = int(math.ceil(1.2 * particlesize))
        return boxsize + boxsize % 2

    @staticmethod
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
        return "Box size is too large!"

    def update_box_size_labels(self, *args, **kwargs):
        try:
            angpix = float(self.angpix_entry.get())
        except ValueError:
            # Can't update any of the labels without angpix
            self.mask_diameter_px.config(text="= NNN px")
            self.box_size_in_angstrom.config(text=u"= NNN \u212B")
            self.extract_angpix.config(text=u"= NNN \u212B/px")
            return
        try:
            mask_diameter = float(self.mask_diameter_entry.get())
            mask_diameter_px = mask_diameter / angpix
            self.mask_diameter_px.config(text="= {:.1f} px".format(mask_diameter_px))
        except (ValueError, ZeroDivisionError):
            self.mask_diameter_px.config(text="= NNN px")
            # Don't return - an error here doesn't stop us calculating the other labels
        try:
            box_size = float(self.box_size_entry.get())
            box_angpix = angpix * box_size
            self.box_size_in_angstrom.config(text=u"= {:.1f} \u212B".format(box_angpix))
        except ValueError:
            # Can't update these without the box size
            self.box_size_in_angstrom.config(text=u"= NNN \u212B")
            self.extract_angpix.config(text=u"= NNN \u212B/px")
            return
        try:
            extract_small_boxsize = float(self.extract_small_boxsize_entry.get())
            small_box_angpix = box_angpix / extract_small_boxsize
            self.extract_angpix.config(text=u"= {:.3f} \u212B/px".format(small_box_angpix))
        except (ValueError, ZeroDivisionError):
            # Can't update the downscaled pixel size unless the downscaled box size is valid
            self.extract_angpix.config(text=u"= NNN \u212B/px")

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
                angpix = float(self.angpix_entry.get())
                particle_size_angstroms = float(self.particle_max_diam_entry.get())
                particle_size_pixels = particle_size_angstroms / angpix
                mask_diameter = 1.1 * particle_size_angstroms
                box_size = self.calculate_box_size(particle_size_pixels)
                small_boxsize = self.calculate_downscaled_box_size(box_size, angpix)

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
        self.class2d_button.config(state=new_state)
        self.class3d_button.config(state=new_state)
        self.particle_max_diam_entry.config(state=new_state)
        self.self.particle_min_diam_entry.config(state=new_state)
        self.ref_3d_entry.config(state=new_state)
        # Update the box size controls with care to avoid activating them when we shouldn't
        self.auto_boxsize_button.config(state=new_state)
        if new_state == tk.DISABLED:
            self.mask_diameter_entry.config(state=new_state)
            self.box_size_entry.config(state=new_state)
            self.extract_small_boxsize_entry.config(state=new_state)
        else:
            self.update_box_sizes()
        can_do_second_pass = (
            self.class3d_var.get()
            and len(self.ref_3d_var.get()) == 0
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
        opts.do_class2d                = self.class2d_var.get()
        opts.do_class3d                = self.class3d_var.get()
        opts.do_second_pass            = self.second_pass_var.get()
        opts.do_class2d_pass2          = self.class2d_pass2_var.get()
        opts.do_class3d_pass2          = self.class3d_pass2_var.get()

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

        opts.angpix = try_get(float, self.angpix_entry, "Pixel size")
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

        opts.extract_small_boxsize, opts.extract2_small_boxsize = 2 * (try_get(int, self.extract_small_boxsize_entry, "Down-sampled box size"),)
        opts.extract_downscale, opts.extract2_downscale = 2 * (True,)
        check_positive(opts.extract_small_boxsize, "Down-sampled box size")

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

        # If we have a 3D reference, do a single pass with a large batch size
        opts.autopick_do_LoG = not opts.autopick_3dreference
        if opts.autopick_3dreference:
            opts.autopick_refs_min_distance = opts.autopick_LoG_diam_max * 0.7
            opts.class3d_reference = opts.autopick_3dreference
            opts.do_second_pass = False
        else:
            # No 3D reference - do LoG autopicking in the first pass
            opts.class3d_reference = ''

        # Now set a sensible batch size (leaving batch_size_pass2 at its default 100,000)
        opts.batch_size = 10000 if opts.do_second_pass else 100000

    def save_options(self):
        '''
        Update the full set of options from the values in the GUI, and save them to a file.

        Returns:
            True if the options were valid and saved successfully, otherwise False.
        '''
        try:
            self.fetch_options_from_gui()
            if not self.warnings or tkMessageBox.askokcancel(
                "Warning", "\n".join(self.warnings), icon="warning",
                default=tkMessageBox.CANCEL
            ):
                self.calculate_full_options()
                print(prefix_RELION_IT("Writing all options to {}".format(OPTIONS_FILE)))
                if os.path.isfile(OPTIONS_FILE):
                    print(prefix_RELION_IT("File {0} already exists; renaming old copy to {0}~".format(OPTIONS_FILE)))
                    os.rename(OPTIONS_FILE, OPTIONS_FILE + '~')
                self.options.print_options(OPTIONS_FILE)
                return True
        except Exception as ex:
            tkMessageBox.showerror("Error", ex.message)
            traceback.print_exc()
        return False

    def run_pipeline(self):
        '''
        Update the full set of options from the values in the GUI, close the GUI and run the pipeline.
        '''
        if self.save_options():
            self.main_window.destroy()
            self.options.run_pipeline()


def safe_load_star(filename, max_tries=5, wait=10, expected=[]):
    for _ in range(max_tries):
        try:
            starfile = star.load(filename)
            # Ensure the expected keys are present
            # By descending through the dictionary
            entry = starfile
            for key in expected:
               entry = entry[key]
            return starfile
        except:
            print("safe_load_star is retrying to read: {}, expected key: {}".format(filename, expected))
            time.sleep(wait)
    raise Exception("Failed to read a star file: {}".format(filename))


# Don't get stuck in an infinite loop
def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print(prefix_RELION_IT("{} file no longer exists, exiting now ...".format(RUNNING_FILE)))
        exit(0)


# Allow direct progressing to the second pass
def getSecondPassReference():
    if os.path.isfile(SECONDPASS_REF3D_FILE):
        with open(SECONDPASS_REF3D_FILE, 'r') as readfile:
            filename, angpix = map(str.strip, readfile.readlines())
    else:
        filename, angpix = '', '0'
    return filename, angpix


def getJobName(name_in_script, done_file):
    jobname = None
    # See if we've done this job before,
    # i.e. whether it is in the done_file
    if os.path.isfile(done_file):
        with open(done_file, 'r') as f:
            for elems in map(str.split, f):
                if len(elems) < 3:
                    continue
                if elems[0] == name_in_script:
                    jobname = elems[2]
                    break
    return jobname


def addJob(jobtype, name_in_script, done_file, options, alias=None):
    jobname = getJobName(name_in_script, done_file)

    # If we didn't before, add it now
    if jobname is None:
        command = (
            'relion_pipeliner'
            + ' --addJob {}'.format(jobtype)
            + ' --addJobOptions "{}"'.format(''.join(opt + ';' for opt in options))
        )
        if alias is not None:
            command += ' --setJobAlias "{}"'.format(alias)
        # print("DEBUG: Running " + command)
        os.system(command)

        pipeline = safe_load_star(PIPELINE_STAR, expected=['pipeline_processes', 'rlnPipeLineProcessName'])
        jobname = pipeline['pipeline_processes']['rlnPipeLineProcessName'][-1]

        with open(done_file, 'a') as f:
            f.write('{} = {}\n'.format(name_in_script, jobname))

    # Return the name of the job in the RELION pipeline,
    # e.g. 'Import/job001/'
    return jobname, (jobname is not None)


def RunJobs(jobs, repeat, wait, schedulename):
    os.system(
        'relion_pipeliner'
        + ' --schedule {}'.format(schedulename)
        + ' --repeat {}'.format(repeat)
        + ' --min_wait {}'.format(wait)
        + ' --RunJobs "{}"'.format(' '.join(jobs))
        + ' &'
    )


def WaitForJob(job, time=30):
    '''
    `job`:  name of job (`str`)
    `time`: time to wait in seconds
    '''
    time.sleep(time)
    print(prefix_RELION_IT("waiting for job to finish in {}".format(job)))
    while True:
        pipeline = safe_load_star(
            PIPELINE_STAR, 
            expected=['pipeline_processes', 'rlnPipeLineProcessName']
        )
        processes = pipeline['pipeline_processes']
        try:
            jobnr = processes['rlnPipeLineProcessName'].index(job)
        except ValueError:
            print(prefix_ERROR("cannot find {} in {}".format(job, PIPELINE_STAR)))
            exit(1)

        if int(
            processes['rlnPipeLineProcessStatus'][jobnr]
        ) == 2:
            print(prefix_RELION_IT("job in {} has finished now".format(job)))
            break

        CheckForExit()
        time.sleep(time)


def find_split_job_output(prefix, n, max_digits=6):
    for filename in (
        prefix + str(n).rjust(i, '0') + '.star' for i in range(max_digits)
    ):
        if os.path.isfile(filename):
            return filename


def writeManualPickingGuiFile(particle_diameter):
    if not os.path.isfile('.gui_manualpickrun.job'):
        with open('.gui_manualpickrun.job', 'w') as f:
            f.write(''.join(line + '\n' for line in [
                '{} == {}'.format(question, answer) for question, answer in [
                    ('job_type',                   '3'),
                    ('Pixel size (A)',             '-1'),
                    ('Black value:',               '0'),
                    ('Blue value:',                '0'),
                    ('MetaDataLabel for color:',   'rlnParticleSelectZScore'),
                    ('Scale for CTF image:',       '1'),
                    ('Particle diameter (A):',     str(particle_diameter)),
                    ('Blue<>red color particles?', 'No'),
                    ('Highpass filter (A)',        '-1'),
                    ('Lowpass filter (A)',         '20'),
                    ('Scale for micrographs:',     '0.2'),
                    ('Red value:',                 '2'),
                    ('Sigma contrast:',            '3'),
                    ('White value:',               '0'),
                ]
            ]))


def findBestClass(model_star_file, use_resol=True):
    '''
    Identify the best class given `model_star_file` (a `str`).
    Return the index of the best class (an `int`), the best resolution (a `float`), and the pixel size.
    '''
    model_star = safe_load_star(model_star_file)
    best_class, best_size, best_resol = 0, 0, 999

    if use_resol:
        def isbetter(curr_size, curr_resol, best_size, best_resol):
            return curr_resol < best_resol or (
                curr_resol == best_resol and curr_size > best_size
            )
    else:
        def isbetter(curr_size, curr_resol, best_size, best_resol):
            return curr_size > best_size or (
                curr_size == best_size and curr_resol < best_resol
            )

    for index, size, resol in zip(
        model_star['model_classes']['rlnReferenceImage'],
        map(float, model_star['model_classes']['rlnClassDistribution']),
        map(float, model_star['model_classes']['rlnEstimatedResolution']),
    ):
        if isbetter(size, resol, best_size, best_resol):
            best_class = index
            best_size = size
            best_resol = resol

    print(prefix_RELION_IT("found best class: {} with class size of {} and resolution of {}".format(
        best_class, best_size, best_resol
    )))
    return best_class, best_resol, model_star['model_general']['rlnPixelSize']


def findOutputModelStar(job_dir):
    try:
        job_star = safe_load_star(
            str(job_dir) + 'job_pipeline.star',
            expected=['pipeline_output_edges', 'rlnPipeLineEdgeToNode']
        )
        for output_file in job_star['pipeline_output_edges']['rlnPipeLineEdgeToNode']:
            if output_file.endswith("_model.star"):
                return output_file
    except:
        return


def main():
    '''
    Run the RELION 3 pipeline.

    Options files given as command line arguments will be opened in order and
    used to update the default options.
    '''
    # Start by parsing arguments
    # (If --help is given, the program will print a usage message and exit)
    parser = argparse.ArgumentParser()
    parser.add_argument("extra_options", nargs="*", metavar="extra_options.py",
                        help="Python files containing options for relion_it.py")
    parser.add_argument("--gui", action="store_true", help="launch a simple GUI to set options")
    parser.add_argument("--continue", action="store_true", dest="continue_",
                        help="continue a previous run by loading options from ./relion_it_options.py")
    args = parser.parse_args()

    for msg in [
        '-------------------------------------------------------------------------------------------------------------------',
        'script for automated, on-the-fly single-particle analysis in RELION (>= 3.1)',
        'authors: Sjors H.W. Scheres, Takanori Nakane & Colin M. Palmer',
        '',
        'usage: ./relion_it.py [extra_options.py [extra_options2.py ....] ] [--gui] [--continue]',
        '',
        'this script will check whether processes are still running using files with names starting with RUNNING',
        '  you can restart this script after stopping previous processes by deleting all RUNNING files',
        'this script keeps track of already submitted jobs in a filed called {}'.format(SETUP_CHECK_FILE),
        '  upon a restart, jobs present in this file will be continued (for preprocessing), or ignored when already finished',
        'if you would like to re-do a specific job from scratch (e.g. because you changed its parameters)',
        '  remove that job, and those that depend on it, from the {}'.format(SETUP_CHECK_FILE),
        '-------------------------------------------------------------------------------------------------------------------',
        '',
    ]:
        print(prefix_RELION_IT(msg))

    # Exit in case another version of this script is running.
    if os.path.isfile(RUNNING_FILE):
        print(prefix_RELION_IT(' '.join((
            'ERROR: {} is already present: delete this file and make sure no other copy of this script is running.'.format(RUNNING_FILE),
            'Exiting now ...'
        ))))
        exit(0)

    # Exit in case the preprocessing pipeliners are still running.
    for checkfile in ('RUNNING_PIPELINER_' + x for x in (
        PREPROCESS_SCHEDULE_PASS1, PREPROCESS_SCHEDULE_PASS2
    )):
        if os.path.isfile(checkfile):
            print(prefix_RELION_IT(' '.join((
                'ERROR: {} is already present: delete this file and make sure no relion_pipeliner job is still running.'.format(checkfile),
                'Exiting now ...'
            ))))
            exit(0)

    if args.continue_:
        print(prefix_RELION_IT(' '.join((
            'continuing a previous run.',
            'Options will be loaded from ./relion_it_options.py'
        ))))
        args.extra_options.append(OPTIONS_FILE)

    opts = RelionItOptions()
    for user_opt_file in args.extra_options:
        print(prefix_RELION_IT('reading options from {}'.format(user_opt_file)))
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)

    if args.gui:
        print(prefix_RELION_IT('launching GUI...'))
        tk_root = tk.Tk()
        tk_root.title("relion_it.py setup")
        RelionItGui(tk_root, opts)
        tk_root.mainloop()
    else:
        opts.run_pipeline()


if __name__ == "__main__":
    main()
