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
the relevant jobs (`autopick_job_pass2' and those following) from the `RELION_IT_SUBMITTED_JOBS' file, then restarting the
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
from __future__ import division  # So that / does floating-point division

import argparse
import glob
import inspect
import os
import runpy
import time
import re
import itertools
import star
from relion_it_gui import RelionItGui, OPTIONS_FILE

# Constants
PIPELINE_STAR = 'default_pipeline.star'
RUNNING_FILE = 'RUNNING_RELION_IT'
SETUP_CHECK_FILE = 'RELION_IT_SUBMITTED_JOBS'
PREPROCESS_SCHEDULE_PASSES = ['PREPROCESS_PASS' + str(n + 1) for n in range(2)]
REF3D_FILES = ['RELION_IT_3DREF_PASS' + str(n + 1) for n in range(2)]
GUI_PROJECT_FILE = '.gui_projectdir'


def is_dunder_name(name):
    return name.startswith('__') and name.endswith('__')


def prefix_RELION_IT(msg):
    return ' RELION_IT: ' + msg


def prefix_ERROR(msg):
    return ' ERROR: ' + msg


def yesno(x):
    return 'Yes' if x else 'No'


def withRenaming(method):
    def newmethod(self, filename):
        if os.path.isfile(filename):
            print(prefix_RELION_IT("File {0} already exists. Renaming old copy to {0}~".format(filename)))
            os.rename(filename, filename + '~')
        method(self, filename)
    return newmethod


class QAList(list):
    ''' A thin wrapper around the built-in `list` type.
    '''

    def formatted_pairs(self):
        return ['{} == {}'.format(question, answer) for question, answer in self]

    def __add__(self, other):
        return QAList(list(self) + list(other))


class Command(object):

    def __init__(self, executable, options):
        self.executable = executable
        self.options = options

    def __str__(self):
        return ' '.join([self.executable] + self.formatted_pairs())

    def formatted_pairs(self):
        return ['{} {}'.format(k, v) for k, v in self.options]


class Job(object):

    jobtype = ''
    name_in_script = ''
    alias = None

    def __init__(self, options):
        self.options = options

    def add_to_pipeline(self):
        addJob(self.jobtype, self.name_in_script, self.options, self.alias)

    def findOutputModelStar(self):
        keys = ['pipeline_output_edges', 'rlnPipeLineEdgeToNode']
        job_star = star.load(self.name + 'job_pipeline.star', keys)

        for output_file in star.recursivelydescend(job_star, keys):
            if output_file.endswith("_model.star"):
                return output_file

        print(prefix_RELION_IT("{} {} does not contain expected output maps.".format(self.description.capitalize(), self.name)))
        print(prefix_RELION_IT("This job should have finished, but you may continue it from the GUI."))
        raise Exception("ERROR! Quitting the pipeline.")
        # TODO Make more robust


class ImportJob(Job):

    jobtype = 'Import'
    name_in_script = 'import_job'

    def __init__(self, options):
        self.options = QAList([
            ('Raw input files:',               options.import_images),
            ('Import raw movies/micrographs?', 'Yes'),
            ('Pixel size (Angstrom):',         options.angpix),
            ('Voltage (kV):',                  options.voltage),
            ('Spherical aberration (mm):',     options.Cs),
            ('Amplitude contrast:',            options.ampl_contrast),
            ('Are these multi-frame movies?',  yesno(options.images_are_movies)),
        ])

        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class MotionCorJob(Job):

    jobtype = 'MotionCorr'
    name_in_script = 'motioncor_job'

    def __init__(self, options, movies_star_file):
        self.options = QAList([
            ('Input movies STAR file:',    movies_star_file),
            ('MOTIONCOR2 executable:',     options.motioncor_exe),
            ('Defect file:',               options.motioncor_defectfile),
            ('Gain-reference image:',      options.motioncor_gainreference),
            ('Gain flip:',                 options.motioncor_gainflip),
            ('Gain rotation:',             options.motioncor_gainrot),
            ('Do dose-weighting?',         'Yes'),
            ('Dose per frame (e/A2):',     options.motioncor_doseperframe),
            ('Number of patches X:',       options.motioncor_patches_x),
            ('Number of patches Y:',       options.motioncor_patches_y),
            ('Bfactor:',                   options.motioncor_bfactor),
            ('Binning factor:',            options.motioncor_binning),
            ('Which GPUs to use:',         ' '.join(map(str, options.motioncor_gpu))),
            ('Other MOTIONCOR2 arguments', options.motioncor_other_args),
            ('Number of threads:',         options.motioncor_threads),
            ('Number of MPI procs:',       options.motioncor_mpi),
            ('Additional arguments:',      ' '.join('{} {}'.format(k, v) for k, v in [
                ('--eer_upsampling', options.eer_upsampling),
                ('--eer_grouping',   options.eer_grouping),
            ])),
            ('Use RELION\'s own implementation?', yesno(options.motioncor_do_own)),
        ] + ([
            ('Save sum of power spectra?', yesno(options.ctf_software == 'ctffind')),
        ] if options.motioncor_do_own else [])) + (
            options.queue_options if options.motioncor_submit_to_queue else []
        )

        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class CTFFindJob(Job):

    jobtype = 'CtfFind'
    name_in_script = 'ctffind_job'

    def __init__(self, options, micrographs_star_file):
        gctf_options = [
            ('Ignore \'Searches\' parameters?', yesno(options.gctf_ignore_search_params)),
            ('Perform equi-phase averaging?',   yesno(options.gctf_do_EPA)),
        ] if options.ctf_software == 'gctf' else []

        self.options = QAList([
            ('Amount of astigmatism (A):', options.ctffind_astigmatism),
            ('FFT box size (pix):',        options.ctffind_boxsize),
            ('Maximum defocus value (A):', options.ctffind_defocus_max),
            ('Minimum defocus value (A):', options.ctffind_defocus_min),
            ('Defocus step size (A):',     options.ctffind_defocus_step),
            ('Maximum resolution (A):',    options.ctffind_maxres),
            ('Minimum resolution (A):',    options.ctffind_minres),
            ('Gctf executable:',           options.gctf_exe),
            ('Which GPUs to use:',         ' '.join(map(str, options.gctf_gpu))),
            ('CTFFIND-4.1 executable:',    options.ctffind4_exe),
            ('Number of MPI procs:',       options.ctffind_mpi),
            ('Input micrographs STAR file:', micrographs_star_file),
            ('Use CTFFIND-4.1?',                       yesno(options.ctf_software == 'ctffind')),
            ('Use Gctf instead?',                      yesno(options.ctf_software == 'gctf')),
            ('Use power spectra from MotionCorr job?', yesno(options.ctf_software == 'ctffind')),
        ] + gctf_options + [
            ('Estimate phase shifts?', yesno(options.ctffind_do_phaseshift)),
        ]) + (options.queue_options if options.ctffind_submit_to_queue else [])

        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class AutopickJob(Job):

    jobtype = 'AutoPick'

    def __init__(self, options, passindex, micrographs_star_file):
        self.options = QAList([
            ('Input micrographs for autopick:',      micrographs_star_file),
            ('Min. diameter for LoG filter (A)',     options.autopick_LoG_diam_min),
            ('Max. diameter for LoG filter (A)',     options.autopick_LoG_diam_max),
            ('Maximum resolution to consider (A)',   options.autopick_lowpass),
            ('Adjust default threshold (stddev):',   options.autopick_LoG_adjust_threshold),
            ('Upper threshold (stddev):',            options.autopick_LoG_upper_threshold),
            ('2D references:',                       options.autopick_2dreferences),
            ('3D reference:',                        options.autopick_3dreference),
            ('Symmetry:',                            options.autopick_3dref_symmetry),
            ('Pixel size in references (A)',         options.autopick_ref_angpix),
            ('3D angular sampling:',                 options.autopick_3dref_sampling),
            ('In-plane angular sampling (deg)',      options.autopick_inplane_sampling),
            ('Picking threshold:',                   options.autopick_refs_threshold),
            ('Minimum inter-particle distance (A):', options.autopick_refs_min_distance),
            ('Mask diameter (A)',                    options.autopick_refs_mask_diam),
            ('Maximum stddev noise:',                options.autopick_stddev_noise),
            ('Minimum avg noise:',                   options.autopick_avg_noise),
            ('Shrink factor:',                       options.autopick_shrink_factor),
            ('Which GPUs to use:',                   ' '.join(map(str, options.autopick_gpus))),
            ('Additional arguments:',                options.autopick_other_args),
            ('Number of MPI procs:',                 options.autopick_mpi),
            ('OR: provide a 3D reference?',          yesno(options.autopick_3dreference)),
            ('OR: use Laplacian-of-Gaussian?',       yesno(options.autopick_do_LoG)),
            ('Are References CTF corrected?',        yesno(options.autopick_refs_are_ctf_corrected)),
            ('References have inverted contrast?',   yesno(options.autopick_refs_have_inverted_contrast)),
            ('Ignore CTFs until first peak?',        yesno(options.autopick_refs_ignore_ctf1stpeak)),
            ('Use GPU acceleration?',                yesno(options.autopick_gpus and not options.autopick_do_LoG)),
        ])

        if options.autopick_submit_to_queue:
            self.options += options.queue_options

        self.name_in_script = 'autopick_job_pass{}'.format(passindex + 1)
        self.alias='pass {}'.format(passindex + 1)

        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class ExtractJob(Job):

    jobtype = 'Extract'

    def __init__(self, options, passindex, coordinates_star_file, micrographs_star_file):
        # Extract options
        self.options = QAList([
            ('Input coordinates:',                coordinates_star_file),
            ('micrograph STAR file:',             micrographs_star_file),
            ('Diameter background circle (pix):', options.extract_bg_diameter),
            ('Particle box size (pix):',          options.extract_boxsize),
            ('Number of MPI procs:',              options.extract_mpi),
        ])

        if options.extract_downscale[passindex]:
            self.options += QAList([
                ('Rescale particles?',       'Yes'),
                ('Re-scaled size (pixels):', options.extract_small_boxsize[passindex]),
            ])

        if options.extract_submit_to_queue:
            self.options += options.queue_options

        self.name_in_script = 'extract_job_pass{}'.format(passindex + 1)
        self.alias = 'pass {}'.format(passindex + 1)
        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class SplitJob(Job):
    ''' Select job to split the particles into batches.
    '''

    jobtype = 'Select'

    def __init__(self, options, passindex, particles_star_file):
        self.options = QAList([
            ('OR select from particles.star:', particles_star_file),
            ('OR: split into subsets?',        'Yes'),
            ('OR: number of subsets:',         -1),
            ('Subset size:',                   options.batchsize[passindex])
        ])

        self.name_in_script = 'split_job_pass{}'.format(passindex + 1)
        self.alias = 'into {}'.format(options.batchsize[passindex])
        self.add_to_pipeline()
        self.name = most_recently_submitted_job()
        options.jobs.append(self.name)


class DiscardJob(Job):
    '''A Select job to get rid of particles with outlier average/stddev values...
    '''

    jobtype = 'Select'
    description = 'job to discard based on image statistics'

    def __init__(self, options, passindex, particles_star_file):
        self.options = QAList([
            ('OR select from particles.star:',     particles_star_file),
            ('OR: select on image statistics?',    'Yes'),
            ('Sigma-value for discarding images:', options.discard_sigma),
            ('Metadata label for images:',         'rlnImageName'),
        ])

        if options.discard_submit_to_queue:
            self.options += options.queue_options

        self.name_in_script = 'discard_job_pass{}'.format(passindex + 1)
        self.add_to_pipeline()
        self.name = check_whether_job_done(self.name_in_script) or most_recently_submitted_job()


class InimodelJob(Job):

    jobtype = 'InitialModel'
    name_in_script = 'inimodel'
    description = 'initial model generation'

    def __init__(self, options, particles_star_file):
        self.options = QAList([
            ('Input images STAR file:',                   particles_star_file),
            ('Symmetry:',                                 options.symmetry),
            ('Mask diameter (A):',                        options.mask_diameter),
            ('Number of classes:',                        options.inimodel_nr_classes),
            ('Initial angular sampling (deg):',           options.inimodel_angle_step),
            ('Offset search range (pix):',                options.inimodel_offset_range),
            ('Offset search step (pix):',                 options.inimodel_offset_step),
            ('Number of initial iterations:',             options.inimodel_nr_iter_initial),
            ('Number of in-between iterations:',          options.inimodel_nr_iter_inbetween),
            ('Number of final iterations:',               options.inimodel_nr_iter_final),
            ('Write-out frequency (iter):',               options.inimodel_freq_writeout),
            ('Initial resolution (A):',                   options.inimodel_resol_ini),
            ('Final resolution (A):',                     options.inimodel_resol_final),
            ('Initial mini-batch size:',                  options.inimodel_batchsize_ini),
            ('Final mini-batch size:',                    options.inimodel_batchsize_final),
            ('Increased noise variance half-life:',       options.inimodel_sigmafudge_halflife),
            ('Number of pooled particles:',               1),
            ('Which GPUs to use:',                        ' '.join(map(str, options.refine_gpus))),
            ('Number of MPI procs:',                      options.refine_mpi_procs),
            ('Number of threads:',                        options.refine_threads),
            ('Copy particles to scratch directory:',      options.refine_scratch_disk),
            ('Additional arguments:',                     options.inimodel_other_args),
            ('Flatten and enforce non-negative solvent?', yesno(options.inimodel_solvent_flatten)),
            ('Skip padding?',                             yesno(options.refine_skip_padding)),
            ('Use GPU acceleration?',                     yesno(options.refine_gpus)),
            ('Ignore CTFs until first peak?',             yesno(options.inimodel_ctf_ign1stpeak)),
            ('Pre-read all particles into RAM?',          yesno(options.refine_preread_images)),
        ])

        if options.refine_submit_to_queue:
            self.options += options.queue_options

        self.add_to_pipeline()
        self.name = check_whether_job_done(self.name_in_script) or most_recently_submitted_job()


class Class2DJob(Job):

    jobtype = 'Class2D'
    description = '2D classification'

    def __init__(self, options, passindex, batchindex, particles_star_file, batchsize):
        self.options = QAList([
            ('Input images STAR file:',                 particles_star_file),
            ('Number of classes:',                      options.class2d_nr_classes),
            ('Mask diameter (A):',                      options.mask_diameter),
            ('Number of iterations:',                   options.class2d_nr_iter),
            ('Angular search range - psi (deg):',       options.class2d_angle_step),
            ('Offset search range (pix):',              options.class2d_offset_range),
            ('Offset search step (pix):',               options.class2d_offset_step),
            ('Number of pooled particles:',             options.refine_nr_pool),
            ('Which GPUs to use:',                      ' '.join(map(str, options.refine_gpus))),
            ('Number of MPI procs:',                    options.refine_mpi_procs),
            ('Number of threads:',                      options.refine_threads),
            ('Copy particles to scratch directory:',    options.refine_scratch_disk),
            ('Additional arguments:',                   options.class2d_other_args),
            ('Use fast subsets (for large data sets)?', yesno(batchsize > options.refine_batchsize_for_fast_subsets)),
            ('Use GPU acceleration?',                   yesno(options.refine_gpus)),
            ('Ignore CTFs until first peak?',           yesno(options.class2d_ctf_ign1stpeak)),
            ('Pre-read all particles into RAM?',        yesno(options.refine_preread_images)),
        ])

        if options.refine_submit_to_queue:
            self.options += options.queue_options

        self.name_in_script = 'class2d_job_pass{}_batch_{:03d}'.format(passindex + 1, batchindex + 1)
        # e.g. 'class2d_job_pass1_batch_001'
        self.alias = 'pass{}_batch_{:03d}'.format(passindex + 1, batchindex + 1)
        self.add_to_pipeline()
        self.name = check_whether_job_done(self.name_in_script) or most_recently_submitted_job()


class Class3DJob(Job):

    jobtype = 'Class3D'
    description = '3D classification'

    def __init__(self, options, passindex, batchindex, particles_star_file, batchsize):
        self.options = QAList([
            ('Input images STAR file:',                 particles_star_file),
            ('Reference map:',                          options.class3d_reference),
            ('Initial low-pass filter (A):',            options.class3d_ini_lowpass),
            ('Symmetry:',                               options.symmetry),
            ('Regularisation parameter T:',             options.class3d_T_value),
            ('Reference mask (optional):',              options.class3d_reference_mask),
            ('Number of classes:',                      options.class3d_nr_classes),
            ('Mask diameter (A):',                      options.mask_diameter),
            ('Number of iterations:',                   options.class3d_nr_iter),
            ('Angular sampling interval:',              '{} degrees'.format(options.class3d_angle_step)),
            ('Offset search range (pix):',              options.class3d_offset_range),
            ('Offset search step (pix):',               options.class3d_offset_step),
            ('Number of pooled particles:',             options.refine_nr_pool),
            ('Which GPUs to use:',                      ' '.join(map(str, options.refine_gpus))),
            ('Number of MPI procs:',                    options.refine_mpi_procs),
            ('Number of threads:',                      options.refine_threads),
            ('Copy particles to scratch directory:',    options.refine_scratch_disk),
            ('Additional arguments:',                   options.class3d_other_args),
            ('Use fast subsets (for large data sets)?', yesno(batchsize > options.refine_batchsize_for_fast_subsets)),
            ('Ref. map is on absolute greyscale?',      yesno(options.class3d_ref_is_correct_greyscale)),
            ('Has reference been CTF-corrected?',       yesno(options.class3d_ref_is_ctf_corrected)),
            ('Skip padding?',                           yesno(options.refine_skip_padding)),
            ('Use GPU acceleration?',                   yesno(options.refine_gpus)),
            ('Ignore CTFs until first peak?',           yesno(options.class3d_ctf_ign1stpeak)),
            ('Pre-read all particles into RAM?',        yesno(options.refine_preread_images)),
        ])

        if options.refine_submit_to_queue:
            self.options += options.queue_options

        self.name_in_script = 'class3d_job_pass{}_batch_{:03d}'.format(passindex + 1, batchindex + 1)
        self.alias = 'pass{}_batch_{:03d}'.format(passindex + 1, batchindex + 1)
        self.add_to_pipeline()
        self.name = check_whether_job_done(self.name_in_script) or most_recently_submitted_job()



class Batch(object):
    
    def __init__(self, job, schedule, batchsize):
        self.job = job
        self.schedule = schedule
        self.batchsize = batchsize
    
    def run_and_wait(self):
        RunJobs([self.job.name], 1, 1, self.schedule)
        print(prefix_RELION_IT("submitted {} with {} particles in {}".format(
            self.job.description, self.batchsize, self.job.name
        )))
        WaitForJob(self.job)
        return True


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
    # Spherical aberration: Polara = 2.0; Talos/Krios = 2.7; some Cryo-ARM = 1.4
    Cs = 1.4

    ### Import images (Linux wild card; movies as *.mrc, *.mrcs, *.tiff or *.tif; single-frame micrographs as *.mrc)
    import_images = 'Movies/*.tiff'
    # Are these multi-frame movies? Set to False for single-frame micrographs (and motion-correction will be skipped)
    images_are_movies = True

    ### Motion-correction parameters
    # Dose (electrons per square Angstrom per fraction)
    motioncor_doseperframe = 1.277
    # Gain-reference image in MRC format (only necessary if input movies are not yet gain-corrected, e.g. compressed TIFFs from K2)
    motioncor_gainreference = 'Movies/gain.mrc'
    # EER upsampling (1 = 4K, 2 = 8K).
    # If you use 8K rendering, the pixel size (angpix) MUST be the half of the physical pixel size and the motioncor_binning should be 2.
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
    # To pick fewer particles, use positive values: (0,1]
    # To pick more particles, use negative values: [-1,0)
    autopick_LoG_adjust_threshold = 0.0
    autopick_LoG_upper_threshold = 999.0
    #
    # OR: 2D references for reference-based picking (when autopick_do_LoG = False)
    autopick_2dreferences = ''
    # OR: 3D references for reference-based picking (when autopick_do_LoG = False)
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
    # Down-scale particles upon extraction? (1st pass, 2nd pass)
    extract_downscale = (False, False)
    # Box size of the down-scaled particles (in pixels) (1st pass, 2nd pass)
    extract_small_boxsize = (64, 128)

    # 2D / 3D classification of extracted particles
    # Perform 2D classification? (1st pass, 2nd pass)
    do_class2d = (True, True)
    # Perform 3D classification? (1st pass, 2nd pass)
    do_class3d = (True, False)
    # Repeat 2D and/or 3D-classification for batches of this many particles  (1st pass, 2nd pass)
    batchsize = (10000, 100000)
    # Number of 2D classes to use
    class2d_nr_classes  = 50
    # Diameter of the mask used for 2D/3D classification (in Angstrom)
    mask_diameter = 190
    # Symmetry group (when using SGD for initial model generation, C1 may work best)
    symmetry = 'C1'

    ### 3D-classification parameters
    # Number of 3D classes to use
    class3d_nr_classes = 4
    # Initial reference model
    class3d_reference = ''
    # Is reference on correct greyscale?
    class3d_ref_is_correct_greyscale = False
    # Has the initial reference been CTF-corrected?
    class3d_ref_is_ctf_corrected = True
    # Initial lowpass filter on reference
    class3d_ini_lowpass = 40

    ### Use the largest 3D class from the first batch as a 3D reference for a second pass of autopicking?
    # (only when `do_class3d[0]`)
    # Only move on to template-based autopicking if the 3D references achieves this resolution (in A)
    minimum_resolution_3dref_2ndpass = 20

    do_pass = [True, True]

    ###################################################################################
    ############ The following parameters can often be kept the same for a given setup
    ###################################################################################

    ### Repeat settings for entire pipeline
    # Repeat the pre-processing runs this many times (or until RUNNING_PIPELINER_default_PREPROCESS file is deleted)
    preprocess_repeat_max_iters = 999
    # Wait at least this many minutes between each repeat cycle
    preprocess_repeat_wait = 1
    ### Stop after CTF estimation?
    # (skip autopicking, extraction, 2D/3D classification, etc)
    stop_after_ctf_estimation = False
    # How many minutes should there be between checks
    # that enough particles have been extracted for a new batch of 2D classification?
    batch_repeat_time = 1

    ### Motion-correction parameters
    # Use RELION's own (CPU-only) motion-correction instead of the UCSF implementation?
    motioncor_do_own = True
    # The number of threads (only for RELION's own implementation)
    #   Optimal when nr_movie_frames / nr_threads is an integer
    motioncor_threads = 6
    # Exectutable of UCSF MotionCor2
    motioncor_exe = '/public/EM/MOTIONCOR2/MotionCor2'
    # On which GPU(s) to execute UCSF MotionCor2
    motioncor_gpu = [0]
    # Number of MPI processes for motion correction:
    motioncor_mpi = 4
    # Local motion-estimation patches for MotionCor2
    motioncor_patches_x, motioncor_patches_y = 4, 4
    # B-factor in A^2 for downweighting of high-spatial frequencies
    motioncor_bfactor = 150
    # Use binning=2 for super-resolution movies
    motioncor_binning = 1
    # Provide a defect file for your camera if you have one
    motioncor_defectfile = ''
    # Gain reference orientation w.r.t movies
    #   (if input movies are not yet gain-corrected, e.g. TIFFs)
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
    # Which software to use for CTF correction
    # 'ctffind' for Alexis Rohou's CTFFIND4 (CPU-only)
    # 'gctf' for Kai Zhang's Gctf (GPU-accelerated)
    ctf_software = 'ctffind'
    # Estimate phase shifts (for Volta phase plate [VPP] data)
    ctffind_do_phaseshift = False
    # CTFFind options
    # CTFFIND4 executable
    ctffind4_exe = '/public/EM/ctffind/ctffind.exe'
    # How many MPI processes to use for running CTF estimation?
    ctffind_mpi = 8
    # Submit CTF estimation job to the cluster?
    ctffind_submit_to_queue = False
    # Gctf executable
    gctf_exe = '/public/EM/Gctf/bin/Gctf'
    # On which GPU(s) to execute Gctf
    gctf_gpu = [0]
    # For Gctf: ignore parameters on the 'Searches' tab?
    gctf_ignore_search_params = True
    # For Gctf: perform equi-phase averaging?
    gctf_do_EPA = True

    ### Autopick parameters
    # Which GPU(s) to use for autopicking
    # To not use GPU-acceleration for autopicking, let `autopick_gpus = []`
    autopick_gpus = [0]
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
    # 3D angular sampling for generating projections of the 3D reference for autopicking
    # (30 degrees is usually enough)
    autopick_3dref_sampling = 30
    # Pixel size in the provided 2D/3D references
    # Having this negative will cause Relion to use the pixel size of the motion-corrected movies
    autopick_ref_angpix = -1

    ### Extract parameters
    # Diameter for background normalisation (pixels)
    # Having this negative value will cause Relion to use the default box size (75%)
    extract_bg_diameter = -1
    # How many MPI processes to use for running particle extraction?
    extract_mpi = 1
    # Submit Extract job to the cluster?
    extract_submit_to_queue = False

    ## Discard particles based on average/stddev values? (this may be important for SGD initial model generation)
    do_discard_on_image_statistics = False
    # Discard images whose average pixel value is more than this many standard deviations (sigma) away from the ensemble average.
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
    # Which GPU to use (different from GPU used for pre-processing?)
    # To not use GPU-acceleration for refinement, let `refine_gpus = []`
    refine_gpus = [1]
    # How many MPI processes to use
    refine_mpi_procs = 1
    # How many threads to use
    refine_threads = 6
    # Skip padding?
    refine_skip_padding = False
    # Submit jobs to the cluster?
    refine_submit_to_queue = False
    # Use fast subsets in 2D/3D classification when batchsize is bigger than this
    refine_batchsize_for_fast_subsets = 10000

    ### 2D classification parameters
    # Wait with the first 2D classification batch until at least this many particles are extracted
    minimum_batchsize = 10000
    # Number of iterations to perform in 2D classification
    # Must be at least 20 for fast subsets
    class2d_nr_iter = 20
    # Rotational search step (degrees)
    class2d_angle_step = 6
    # Offset search range (pixels)
    class2d_offset_range = 5
    # Offset search step (pixels)
    class2d_offset_step = 1
    # Option to ignore the CTFs until their first peak
    # This option is recommended if all particles go into very few classes
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
    # Regularisation parameter T
    class3d_T_value = 4
    # Angular sampling step
    class3d_angle_step = 7.5
    # Offset search range (pixels)
    class3d_offset_range = 5
    # Offset search step (pixels)
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
    # Initial angular sampling (degrees)
    inimodel_angle_step = 15
    # Initial search range (pixels)
    inimodel_offset_range = 6
    # Initial offset search step (pixels)
    inimodel_offset_step = 2
    # Number of initial iterations
    inimodel_nr_iter_initial = 50
    # Number of in-between iterations
    inimodel_nr_iter_inbetween = 200
    # Number of final iterations
    inimodel_nr_iter_final = 50
    # Frequency to write out information
    inimodel_freq_writeout = 10
    # Initial resolution (A)
    inimodel_resol_ini = 35
    # Final resolution (A)
    inimodel_resol_final = 15
    # Initial mini-batch size
    inimodel_batchsize_ini = 100
    # Final mini-batch size
    inimodel_batchsize_final = 500
    # Increased noise variance half-life (off, i.e. -1, by default; values of ~1000 have been observed to be useful in difficult cases)
    inimodel_sigmafudge_halflife = -1
    # Additional arguments to pass to relion_refine (skip annealing to get rid of outlier particles)
    inimodel_other_args = '--sgd_skip_anneal'

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
        for key, value in [
            (key, value) for key, value in d.items() if not is_dunder_name(key)
            # exclude __name__, __builtins__, etc.
        ]:
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(prefix_RELION_IT("Unrecognised option {}".format(repr(key))))

    @withRenaming
    def write_to(self, filename):
        '''
        Write the current options to `filename`.

        This method writes the options in the same format as they are read,
        allowing options to be written to a file and re-used.

        Args:
            `filename`: A file object (optional). If supplied, options will be
                written to this file, otherwise they will be printed to
                sys.stdout.

        Raises:
            ValueError: If there is a problem printing the options.
        '''
        # NOTE The writing to stdout mentioned in the above docstring is not implemented.
        # NOTE This method requires that every assignment statement
        # between '### End of options' and '### General parameters'
        # occupy only a single line.
        with open(filename, 'w') as f:
            f.writelines(line + '\n' for line in (
                "# Options file for relion_it.py",
                "",
            ))
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

            def outlines():
                # Parse the source code for this class,
                # and write out all comments along with option lines containing new values
                sourcelines, _ = inspect.getsourcelines(self.__class__)
                for line in relevant(sourcelines):
                    if line.startswith('#') or not line:
                        # Print comments and blank lines as-is
                        yield line
                    else:
                        # Assume all other lines define an option name and value.
                        # Replace with new value.
                        m = re.match(assignmentpattern, line)
                        if m:
                            key, oldvalue = map(str.strip, m.groups())
                            try:
                                newvalue = getattr(self, key)
                                yield '{} = {}'.format(key, repr(newvalue))
                                option_names.remove(key)
                            except AttributeError:
                                # This error should not occur.
                                # If it does, there is probably a programming error.
                                raise AttributeError("Unrecognised option name '{}'".format(key))
                if option_names:
                    # This error should not occur.
                    # If it does, there is probably a programming error.
                    raise ValueError("Some options were not written to the output file: {}".format(option_names))

            for line in outlines():
                f.write(line + '\n')


    @property
    def queue_options(self):
        '''List of queue arguments'''
        # Only relevant if `self.motioncor_submit_to_queue`
        return QAList([
            ('Submit to queue?',                  'Yes'),
            ('Queue name:',                       self.queue_name),
            ('Queue submit command:',             self.queue_submit_command),
            ('Standard submission script:',       self.queue_submission_template),
            ('Minimum dedicated cores per node:', self.queue_minimum_dedicated),
        ])

    def run_pipeline(self):
        '''
        Configure and run the RELION 3 pipeline with the given options.
        '''
        # If there is no `default_pipeline.star`, make one.
        # Not really necessary.
        if not os.path.isfile(PIPELINE_STAR):
            with open(PIPELINE_STAR, 'w') as f:
                f.writelines((line + '\n' for line in (
                    'data_pipeline_general',
                    '_rlnPipeLineJobCounter 1',
                )))

        # Create `RUNNING_RELION_IT`.
        # When deleted, this script will stop.
        with open(RUNNING_FILE, 'w'):
            pass

        # Create main GUI project file, so GUI won't ask to set up a project.
        with open(GUI_PROJECT_FILE, 'w'):
            pass

        # Set up GUI file for Manualpick job to allow easy viewing of autopick results.
        if not os.path.isfile('.gui_manualpickrun.job'):
            writeManualPickingGuiFile(
                self.autopick_LoG_diam_min if self.autopick_do_LoG else
                self.autopick_refs_min_distance
            )

        # If we're only doing motioncor and CTF estimation,
        # forget about the second pass and the batch processing.
        if self.stop_after_ctf_estimation:
            self.do_class2d[0], self.do_class3d[0], self.do_pass[1] = False, False, False

        # If `REF3D_FILES[1]` exists,
        # skip the first pass and go straight to the second.
        secondpass_ref3d, secondpass_ref3d_angpix = getSecondPassReference()
        startat = int(self.do_pass[1] and bool(secondpass_ref3d))
        if self.do_pass[1] and secondpass_ref3d:
            for msg in [
                'found {} with angpix {} as a 3D reference for second pass in file {}'.format(
                    secondpass_ref3d, secondpass_ref3d_angpix, REF3D_FILES[1]
                ),
                'if the automatic reference selection turned out to be unsatisfactory,',
                'you can re-run the second pass with another reference by:',
                ' stopping the pipeline with the shell command `rm RUNNING_*`',
                ' updating the reference filename in {}'.format(REF3D_FILES[1]),
                ' deleting relevant jobs (autopick_job_pass2 and followings) in {}'.format(SETUP_CHECK_FILE),
                ' and restarting the pipeline.',
            ]:
                print(prefix_RELION_IT(msg))
            self.class3d_reference = secondpass_ref3d
            self.autopick_3dreference, self.autopick_ref_angpix = secondpass_ref3d, secondpass_ref3d_angpix
            self.autopick_2dreferences = ''
            self.autopick_do_LoG = False

        # Perform two passes through the entire pipeline (PREPROCESS and CLASS2D/3D batches).
        # In the first pass,
        # a 3D reference will be generated.
        # In the second pass,
        # this reference will be used for template-based autopicking.
        for passindex in range(startat, len(filter(bool, self.do_pass))):

            self.jobs = []
            import_job = ImportJob(self)
            movies_star_file = import_job.name + 'movies.star'
            micrographs_star_file = ((
                MotionCorJob(self, movies_star_file).name + 'corrected_'
            ) if self.images_are_movies else import_job.name) + 'micrographs.star'
            ctffind_job = CTFFindJob(self, micrographs_star_file)

            do_2d_classification = self.do_class2d[passindex]
            do_3d_classification = self.do_class3d[passindex]
            do_classification = do_2d_classification or do_3d_classification

            # There is an option to stop on-the-fly processing after CTF estimation
            if not self.stop_after_ctf_estimation:
                autopick_job = AutopickJob(
                    self, passindex,
                    ctffind_job.name + 'micrographs_ctf.star',
                )
                extract_job = ExtractJob(
                    self, passindex,
                    autopick_job.name + 'coords_suffix_autopick.star',
                    ctffind_job.name + 'micrographs_ctf.star',
                )
                if do_classification:
                    split_job = SplitJob(
                        self, passindex, 
                        extract_job.name + 'particles.star',
                    )

            # Now execute the preprocessing pipeliner
            preprocess_schedule = PREPROCESS_SCHEDULE_PASSES[passindex]
            RunJobs(self.jobs, self.preprocess_repeat_max_iters, self.preprocess_repeat_wait, preprocess_schedule)
            print(prefix_RELION_IT('submitted {} pipeliner with {} repeats of the preprocessing jobs'.format(
                preprocess_schedule, self.preprocess_repeat_max_iters
            )))
            running_pipeliner_file = 'RUNNING_PIPELINER_' + preprocess_schedule
            print(prefix_RELION_IT(' '.join((
                'this pipeliner will run in the background of your shell.',
                'You can stop it by deleting the file {}'.format(running_pipeliner_file)
            ))))

            # From now on, process extracted particles in batches for 2D or 3D classification,
            # only perform SGD inimodel for first batch and if no 3D reference is available

            # There is an option to stop here.
            if do_classification:
                # If necessary, rescale the 3D reference in the second pass!
                # TODO rescale initial reference if different from movies?
                if passindex == 1 and any(self.extract_downscale):
                    particles_angpix = self.angpix
                    if self.images_are_movies:
                        particles_angpix *= self.motioncor_binning
                    if self.extract_downscale[1]:
                        particles_angpix *= self.extract_boxsize / self.extract_small_boxsize[1]
                        particles_boxsize = self.extract_small_boxsize[1]
                    else:
                        particles_boxsize = self.extract_boxsize
                    if abs(float(particles_angpix) - float(self.autopick_ref_angpix)) > 0.01:
                        # Rescale the reference for 3D classification
                        print(prefix_RELION_IT('rescaling the 3D reference from pixel size {} to {} and saving the new reference as {}'.format(
                            self.autopick_ref_angpix, particles_angpix, self.class3d_reference
                        )))
                        self.class3d_reference = self.autopick_3dreference.replace('.mrc', '_rescaled.mrc')
                        os.system(str(Command('relion_image_handler', [
                            ('--i',              self.autopick_3dreference),
                            ('--o',              self.class3d_reference),
                            ('--angpix',         self.autopick_ref_angpix),
                            ('--rescale_angpix', particles_angpix),
                            ('--new_box',        particles_boxsize),
                        ])))

                print(prefix_RELION_IT(' '.join((
                    'entering an infinite loop for particle batch-processing.',
                    'You can stop this loop by deleting the file {}'.format(RUNNING_FILE)
                ))))

                # It could be that this is a restart,
                # so check for `previous_batch1_size` in the output directory.
                # Also check for `class2d_job_pass1_batch_001`,
                # in case the first job was not submitted yet.
                try:
                    first_split_file = sorted(glob.glob(split_job.name + 'particles_split*.star'))[0]
                except IndexError:
                    first_split_file = None
                keys = ['particles', 'rlnMicrographName']
                previous_batch1_size = len(star.recursivelydescend(
                    star.load(first_split_file, keys), keys
                )) if (
                    check_whether_job_done('class2d_job_pass1_batch_001') is not None and first_split_file is not None
                ) else 0

                move_to_next_pass = False
                while not move_to_next_pass:
                    new_batches = []
                    for batchindex, batchname in enumerate(sorted(glob.glob(
                        split_job.name + 'particle_split*.star'
                    ))):
                        batch = star.load(batchname, keys)
                        batchsize = len(star.recursivelydescend(batch), keys)
                        rerun_batch1 = batchindex == 0 and batchsize > min(previous_batch1_size, self.minimum_batchsize)
                        if rerun_batch1:
                            previous_batch1_size = batchsize
                        particles_star_file = batchname

                        # The first batch is special:
                        # Perform 2D classification with a batch size smaller than before
                        # (but no less than `minimum_batchsize`).
                        # Keep overwriting in the same output directory.
                        if batchindex == 0 or rerun_batch1:

                            # Discard particles with odd mean/stddev values
                            if self.do_discard_on_image_statistics:
                                discard_job = DiscardJob(self, passindex, batchname)
                                if check_whether_job_done(discard_job.name_in_script) is None or rerun_batch1:
                                    new_batches.append(Batch('DISCARD', discard_job, batchsize).run_and_wait())
                                particles_star_file = discard_job.name + 'particles.star'

                            # 2D classification
                            if do_2d_classification:
                                class2d_job = Class2DJob(self, passindex, batchindex, particles_star_file, batchsize)
                                if check_whether_job_done(class2d_job.name_in_script) is None or rerun_batch1:
                                    new_batches.append(Batch('CLASS2D', class2d_job, batchsize).run_and_wait())

                        # Perform 3D classification
                        if do_3d_classification:
                            # In the first pass and the first (complete) batch
                            # if no reference is provided,
                            # do SGD initial model generation.
                            # Then use this model for subsequent batches.
                            if (
                                passindex == batchindex == 0
                                and not self.class3d_reference
                            ):
                                # Calculate an initial 3D model
                                # by SGD initial model generation
                                inimodel_job = InimodelJob(self, particles_star_file)
                                if check_whether_job_done(inimodel_job.name_in_script) is None:
                                    new_batches.append(Batch('INIMODEL', inimodel_job, batchsize).run_and_wait())

                                sgd_model_star = inimodel_job.findOutputModelStar()
                                # Use the model of the largest class for 3D classification
                                self.class3d_reference, _resol, _angpix = findBestClass(sgd_model_star, use_resol=True)
                                # XXX The comment says 'largest', but the code says 'best resolved'.
                                self.class3d_ref_is_correct_greyscale = True
                                self.class3d_ref_is_ctf_corrected = True

                            if self.class3d_reference:
                                # Do 3D classification
                                class3d_job = Class3DJob(self, passindex, batchindex, particles_star_file, batchsize)
                                if check_whether_job_done(class3d_job.name_in_script) is None or rerun_batch1:
                                    new_batches.append(Batch('CLASS3D', class3d_job, batchsize).run_and_wait())

                                class3d_model_star = class3d_job.findOutputModelStar()
                                best_class3d_class, best_class3d_resol, best_class3d_angpix = findBestClass(class3d_model_star, use_resol=True)

                                # Once the first batch in the first pass is completed,
                                # progress to the second pass.
                                if (
                                    passindex == batchindex == 0
                                    and self.do_pass[1]
                                    and best_class3d_resol < self.minimum_resolution_3dref_2ndpass
                                ):
                                    self.class3d_reference = best_class3d_class
                                    self.autopick_3dreference, self.autopick_ref_angpix = best_class3d_class, best_class3d_angpix
                                    self.autopick_2dreferences = ''
                                    self.autopick_do_LoG = False
                                    self.autopick_3dref_symmetry = self.symmetry

                                    # Remove the first pass's RUNNING_PIPELINER file to stop its PREPROCESS pipeliner
                                    if os.path.isfile(running_pipeliner_file):
                                        print(prefix_RELION_IT(
                                            'removing file {} to stop the pipeliner from the first pass'.format(running_pipeliner_file)
                                        ))
                                        os.remove(running_pipeliner_file)

                                    # Generate a file to indicate we're in the second pass,
                                    # so that the script restarts smoothly
                                    with open(REF3D_FILES[1], 'w') as f:
                                        f.writelines((str(x) + '\n' for x in (
                                            best_class3d_class, best_class3d_angpix
                                        )))

                                    move_to_next_pass = True
                                    print(prefix_RELION_IT(
                                        'moving on to the second pass using {} for template-based autopicking'.format(self.autopick_3dreference)
                                    ))
                                    # break out of the for-loop over the batches
                                    break

                    if not new_batches:
                        CheckForExit()
                        # Don't check `particles.star` too often
                        time.sleep(60 * self.batch_repeat_time)


# Don't get stuck in an infinite loop
def CheckForExit():
    if not os.path.isfile(RUNNING_FILE):
        print(prefix_RELION_IT("file {} no longer exists. Exiting...".format(RUNNING_FILE)))
        exit(0)


# Allow direct progression to the second pass
def getSecondPassReference():
    if os.path.isfile(REF3D_FILES[1]):
        with open(REF3D_FILES[1], 'r') as f:
            class_, angpix = map(str.strip, f)
    else:
        class_, angpix = '', '0'
    return class_, angpix


def SETUP_CHECK_FILE_readwrite(add_job_method):
    def addJob_with_checking(jobtype, name_in_script, options, alias=None):
        '''
        Given a job's name in the script, look for whether it has already been completed.
        If not, add it to the pipeliner and add an entry in `RELION_IT_SUBMITTED_JOBS`.
        '''
        # Check whether this job has already been added to the pipeline.
        if check_whether_job_done(name_in_script) is None:
            # Only add this job if it hasn't already been added.
            add_job_method(jobtype, options, alias)
            # Add an entry in the SETUP_CHECK_FILE
            with open(SETUP_CHECK_FILE, 'a') as f:
                f.write('{} = {}\n'.format(name_in_script, most_recently_submitted_job()))
    return addJob_with_checking


@SETUP_CHECK_FILE_readwrite
def addJob(jobtype, options, alias=None):
    '''
    Given a job type, a list of options, and an optional alias:
    Pass the job to `relion_pipeliner` and return its `jobname`.
    '''
    command = Command('relion_pipeliner', [
        ('--addJob',        jobtype),
        ('--addJobOptions', '"{}"'.format(''.join(
            pair + ';' for pair in options.formatted_pairs()
        ))),
    ])
    if alias is not None:
        command.options.append(('--setJobAlias', '"{}"'.format(alias)))
    os.system(str(command))


def check_whether_job_done(name_in_script, done_file=SETUP_CHECK_FILE):
    '''
    Look for whether this job has already been completed.
    '''
    # See if we've done this job before,
    # i.e. whether it is in the SETUP_CHECK_FILE
    if os.path.isfile(done_file):
        with open(done_file, 'r') as f:
            for tokens in map(str.split, f):
                # Assuming lines follow the pattern 'name_in_script = jobname'
                if len(tokens) >= 3 and tokens[0] == name_in_script:
                    return tokens[2]  # jobname


def most_recently_submitted_job():
    # Return the name of the job in the RELION pipeline most recently submitted,
    # e.g. 'Import/job001/'
    keys = ['pipeline_processes', 'rlnPipeLineProcessName']
    pipeline = star.load(PIPELINE_STAR, keys)
    return star.recursivelydescend(pipeline, keys)[-1]


def RunJobs(jobs, repeat, min_wait, schedule):
    '''
    Pass `jobs` to `relion_pipeliner` with the `--RunJobs` option.
    '''
    os.system(str(Command('relion_pipeliner', [
        ('--schedule', schedule),
        ('--repeat',   repeat),
        ('--min_wait', min_wait),
        ('--RunJobs',  '"{}"'.format(' '.join(jobs))),
    ])) + ' &')  # Background


def WaitForJob(jobname, wait=30):
    '''
    `jobname: str` name of job
    `wait: int` time to wait (wait)
    '''
    def waiting():
        keys = ['pipeline_processes', 'rlnPipeLineProcessName']
        pipeline = star.load(PIPELINE_STAR, keys)
        processes = pipeline['pipeline_processes']
        try:
            jobnr = processes['rlnPipeLineProcessName'].index(jobname)
        except ValueError:
            print(prefix_ERROR("cannot find {} in {}".format(jobname, PIPELINE_STAR)))
            exit(1)
        status = int(processes['rlnPipeLineProcessStatus'][jobnr])
        if status == 2:
            print(prefix_RELION_IT("job in {} has finished now".format(jobname)))
            return False
        return True

    time.sleep(wait)
    print(prefix_RELION_IT("waiting for job to finish in {}".format(jobname)))
    while waiting():
        CheckForExit()
        time.sleep(wait)


def writeManualPickingGuiFile(particle_diameter):
    with open('.gui_manualpickrun.job', 'w') as f:
        f.writelines([line + '\n' for line in QAList([
            ('job_type:',                  3),
            ('Pixel size (A):',            -1),
            ('Black value:',               0),
            ('Blue value:',                0),
            ('MetaDataLabel for color:',   'rlnParticleSelectZScore'),
            ('Scale for CTF image:',       1),
            ('Particle diameter (A):',     particle_diameter),
            ('Blue<>red color particles?', 'No'),
            ('Highpass filter (A):',       -1),
            ('Lowpass filter (A):',        20),
            ('Scale for micrographs:',     0.2),
            ('Red value:',                 2),
            ('Sigma contrast:',            3),
            ('White value:',               0),
        ]).formatted_pairs()])


def findBestClass(model_star_file, use_resol=True):
    '''
    Identify the best class given `model_star_file` (a `str`).
    Return the filename of the best class, its resolution, and the pixel size.
    '''
    model_star = star.load(model_star_file)
    model_classes = model_star['model_classes']

    pair = (
        (lambda size, resol: (1/resol, size)) if use_resol else
        (lambda size, resol: (size, 1/resol))
    )

    best_class, best_size, best_resol = max(zip(
        # We are depending on these 3 iterables that we are zipping to be `list`s.
        # They might actually be `dict`s, in which case we will have to use a dictionary zip.
        model_classes['rlnReferenceImage'],                   # Class
        map(float, model_classes['rlnClassDistribution']),    # Size
        map(float, model_classes['rlnEstimatedResolution']),  # Resolution
    ), key=lambda CSRtriple: pair(CSRtriple[1], CSRtriple[2]))
    # Lexicographical comparison

    print(prefix_RELION_IT("found best class {} of size {} and resolution {}".format(
        best_class, best_size, best_resol
    )))
    return best_class, best_resol, model_star['model_general']['rlnPixelSize']


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
    rows, columns = map(int, os.popen('stty size', 'r').read().split())
    for msg in [
        (columns - 12) * '-',
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
        '  remove that job, and those that depend on it, from {}'.format(SETUP_CHECK_FILE),
        (columns - 12) * '-',
        '',
    ]:
        print(prefix_RELION_IT(msg))

    # Exit in case another version of this script is running.
    if os.path.isfile(RUNNING_FILE):
        print(prefix_RELION_IT(prefix_ERROR(' '.join((
            '{} is already present: delete this file and make sure no other copy of this script is running.'.format(RUNNING_FILE),
            'Exiting now ...'
        )))))
        exit(0)

    # Exit in case the preprocessing pipeliners are still running.
    for checkfile in (
        'RUNNING_PIPELINER_' + pass_ for pass_ in PREPROCESS_SCHEDULE_PASSES
    ):
        if os.path.isfile(checkfile):
            print(prefix_RELION_IT(prefix_ERROR(' '.join((
                '{} is already present: delete this file and make sure no relion_pipeliner job is still running.'.format(checkfile),
                'Exiting now ...'
            )))))
            exit(0)

    if args.continue_:
        print(prefix_RELION_IT(' '.join((
            'continuing a previous run.',
            'Options will be loaded from ./{}'.format(OPTIONS_FILE)
        ))))
        args.extra_options.append(OPTIONS_FILE)

    opts = RelionItOptions()
    for user_opt_file in args.extra_options:
        print(prefix_RELION_IT('reading options from {}'.format(user_opt_file)))
        user_opts = runpy.run_path(user_opt_file)
        opts.update_from(user_opts)

    if args.gui:
        RelionItGui.launch(opts)
    else:
        opts.run_pipeline()


if __name__ == "__main__":
    main()
