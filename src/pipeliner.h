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

#ifndef PIPELINER_H_
#define PIPELINER_H_
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "src/metadata_table.h"
#include "src/pipeline_jobs.h"

class Process {

    public:

    std::string name;
    std::string alias;
    int type;
    int status;
    std::vector<long int> inputNodeList;  // List of Nodes of input to this process
    std::vector<long int> outputNodeList; // List of Nodes of output from this process

    // Constructor
    Process(std::string name, int type, int status, std::string alias="None") {
        this->name = name;
        this->type = type;
        this->status = status;
        this->alias = alias;
    }

    // Destructor
    ~Process() {
        inputNodeList.clear();
        outputNodeList.clear();
    }

    inline std::string alias_or_name() {
        return alias == "None" ? name : alias;
    }

    // All the directory names of the different types of jobs defined inside the pipeline
    static constexpr const char *const IMPORT_NAME       = "Import";       // Import any file as a Node of a given type
    static constexpr const char *const MOTIONCORR_NAME   = "MotionCorr";   // Import any file as a Node of a given type
    static constexpr const char *const CTFFIND_NAME      = "CtfFind";  	   // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
    static constexpr const char *const MANUALPICK_NAME   = "ManualPick";   // Manually pick particle coordinates from micrographs
    static constexpr const char *const AUTOPICK_NAME     = "AutoPick";     // Automatically pick particle coordinates from micrographs, their CTF and 2D references
    static constexpr const char *const EXTRACT_NAME      = "Extract";      // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
    static constexpr const char *const CLASSSELECT_NAME  = "Select"; 	   // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
    static constexpr const char *const CLASS2D_NAME      = "Class2D";      // 2D classification (from input particles)
    static constexpr const char *const CLASS3D_NAME      = "Class3D";      // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
    static constexpr const char *const AUTO3D_NAME       = "Refine3D";     // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
    // static constexpr const char *const POLISH_NAME = "Polish";             // Particle-polishing (from movie-particles)
    static constexpr const char *const MASKCREATE_NAME   = "MaskCreate";   // Process to create masks from input maps
    static constexpr const char *const JOINSTAR_NAME     = "JoinStar";     // Process to create masks from input maps
    static constexpr const char *const SUBTRACT_NAME     = "Subtract";     // Process to subtract projections of parts of the reference from experimental images
    static constexpr const char *const POST_NAME         = "PostProcess";  // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
    static constexpr const char *const RESMAP_NAME       = "LocalRes";     // Local resolution estimation (from unfiltered half-maps and a 3D mask)
    // static constexpr const char *const MOVIEREFINE_NAME = "MovieRefine";   // Movie-particle extraction and refinement combined
    static constexpr const char *const INIMODEL_NAME     = "InitialModel"; // De-novo generation of 3D initial model (using SGD)
    static constexpr const char *const MULTIBODY_NAME    = "MultiBody";    // Multi-body refinement
    static constexpr const char *const MOTIONREFINE_NAME = "Polish";       // Jasenko's motion fitting program for Bayesian polishing (to replace MovieRefine?)
    static constexpr const char *const CTFREFINE_NAME    = "CtfRefine";    // Jasenko's program for defocus and beamtilt optimisation
    static constexpr const char *const EXTERNAL_NAME     = "External";     // For running non-relion programs

    enum Type {
        IMPORT, // Import any file as a Node of a given type
        MOTIONCORR, // Import any file as a Node of a given type
        CTFFIND, // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
        MANUALPICK, // Manually pick particle coordinates from micrographs
        AUTOPICK, // Automatically pick particle coordinates from micrographs, their CTF and 2D references
        EXTRACT, // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
        // SORT,  // Sort particles based on their Z-scores
        CLASSSELECT, // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
        CLASS2D,  // 2D classification (from input particles)
        CLASS3D,  // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
        AUTO3D,  // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
        // POLISH,  // Particle-polishing (from movie-particles)
        MASKCREATE,  // Process to create masks from input maps
        JOINSTAR,  // Process to create masks from input maps
        SUBTRACT,  // Process to subtract projections of parts of the reference from experimental images
        POST,  // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
        RESMAP,  // Local resolution estimation (from unfiltered half-maps and a 3D mask)
        // MOVIEREFINE,  // Movie-particle extraction and refinement combined
        INIMODEL,  // De-novo generation of 3D initial model (using SGD)
        MULTIBODY,  // Multi-body refinement
        MOTIONREFINE,  // Jasenko's motion_refine
        CTFREFINE,  // Jasenko's ctf_refine
        EXTERNAL  // External scripts
    };

    // Status a process may have
    enum ProcessStatus {
        RUNNING,           // (hopefully) running
        SCHEDULED,         // scheduled for future execution
        FINISHED_SUCCESS,  // successfully finished
        FINISHED_FAILURE,  // reported an error
        FINISHED_ABORTED   // aborted by the user
    };

};

enum { DONT_LOCK, DO_LOCK };

#define PIPELINE_HAS_CHANGED ".pipeline_has_changed"

class PipeLine {

    public:

    int job_counter;
    bool do_read_only;

    std::string name;
    std::vector<Node> nodeList; //list of all Nodes in the pipeline
    std::vector<Process> processList; //list of all Processes in the pipeline

    PipeLine() {
        name = "default";
        job_counter = 1;
        do_read_only = false;
    }

    ~PipeLine() {
        clear();
    }

    void clear() {
        nodeList.clear();
        processList.clear();
        job_counter = 1;
    }

    void setName(std::string name) {
        this->name = name;
    }

    // Add a new input Edge to the list
    // Check whether Node with that name already exists in the Node list, and if so update that one
    // The input_for_process will be added to the inputForProcessList of this Node
    //
    void addNewInputEdge(Node &node, long int input_for_process);

    // Add a new output Edge to the list
    // Check whether Node with that name already exists in the Node list, and if so update that one
    // The output_from_process will be added to the outputFromProcessList of this Node
    //
    void addNewOutputEdge(long int output_from_process, Node &node);

    // Check whether Node already exists in the nodeList. If not add and return pointer to new node, otherwise return pointer to existing node
    // Also touch entry in .Nodes directory, use touch_if_not_exist for scheduled jobs
    long int addNode(Node &node, bool touch_if_not_exist = false);

    // Add a new Process to the list (no checks are performed)
    long int addNewProcess(Process &process, bool do_overwrite = false);

    // Find nodes or process (by name or alias)
    long int findNodeByName(std::string name);
    long int findProcessByName(std::string name);
    long int findProcessByAlias(std::string name);

    // Touch each individual Node name in the temporary Nodes directory
    // Return true if Node output file exists and temporary file is written, false otherwise
    bool touchTemporaryNodeFile(Node &node, bool touch_even_if_not_exist=false);
    void touchTemporaryNodeFiles(Process &process);

    // And delete these temporary files
    void deleteTemporaryNodeFile(Node &node);
    void deleteTemporaryNodeFiles(Process &process);

    // Re-make entries of all NodeNames in the hidden .Nodes directory (for file browsing for InputNode I/O)
    void remakeNodeDirectory();

    // Check for process completion by checking for the presence of all outputNode filenames
    // Returns true if any of the running processes has completed, false otherwise
    bool checkProcessCompletion();

    // Get the command line arguments for job
    std::string getCommandLineJob(
        RelionJob &job, int current_job, bool is_main_continue,
        bool is_scheduled, bool do_makedir, bool do_overwrite_current,
        std::vector<std::string> &commands
    ) throw (std::string);

    // Adds job to the pipeline and return the id of the newprocess
    long int addJob(
        RelionJob &job, int as_status, 
        bool do_overwrite, bool do_write_minipipeline = true
    );

    // Runs a job and adds it to the pipeline
    void runJob(
        RelionJob &job, int &current_job, 
        bool only_schedule, bool is_main_continue,
        bool is_scheduled, bool do_overwrite_current
    ) throw (std::string);

    // Adds a scheduled job to the pipeline from the command line (with a name for job type)
    int addScheduledJob(std::string job_type, std::string fn_options);

    // Adds a scheduled job to the pipeline from the command line (with integer job type)
    int addScheduledJob(int job_type, std::string fn_options);

    // Add this RelionJob as scheduled to the pipeline
    int addScheduledJob(RelionJob &job, std::string fn_options="");

    void waitForJobToFinish(int current_job, bool &is_failure, bool &is_abort);

    // Runs a series of scheduled jobs, possibly in a loop, from the command line
    void runScheduledJobs(
        FileName fn_sched, FileName fn_jobids, int nr_repeat, 
        long int minutes_wait, long int minutes_wait_before = 0, long int seconds_wait_after = 10, 
        bool do_overwrite_current = false
    );

    // If I'm deleting this_job from the pipeline, which Nodes and which Processes need to be deleted?
    void deleteJobGetNodesAndProcesses(
        int this_job, bool do_recursive, 
        std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses
    );

    // Given a lists of Nodes and Processes to delete, delete them.
    void deleteNodesAndProcesses(
        std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses
    );

    // Check the presence of a file called RELION_OUTPUT_NODES.star, and add the nodes in that STAR file as output nodes for this job
    void getOutputNodesFromStarFile(int this_job);

    // Changes the status of this_job to finished in the pipeline, throws if job hasn't started yet
    void markAsFinishedJob(int this_job, bool is_failed = false) throw (std::string);

    // Set the alias for a job, return true for success, false otherwise
    void setAliasJob(int this_job, std::string alias) throw (std::string);

    // Make the flowchart for this job
    void makeFlowChart(long int current_job, bool do_display_pdf) throw (std::string);

    // Undelete a JOb from the pipeline
    void undeleteJob(FileName fn_undel);

    // Clean up intermediate files from this_job
    void cleanupJob(int this_job, bool do_harsh) throw (std::string);

    // Clean upintermediate files from all jobs in the pipeline
    void cleanupAllJobs(bool do_harsh) throw (std::string);

    void replaceFilesForImportExportOfScheduledJobs(
        FileName fn_in_dir, FileName fn_out_dir,
        std::vector<std::string> &find_pattern, std::vector<std::string> &replace_pattern
    );

    // Export all scheduled jobs
    void exportAllScheduledJobs(std::string mydir) throw (std::string);

    // Import previously exported jobs
    void importJobs(FileName fn_export);

    // Import a job into the pipeline
    // Return true if at least one process is imported correctly
    bool importPipeline(std::string _name);

    // Write out the pipeline to a STAR file
    void write(
        bool do_lock = false, FileName fn_del="", 
        std::vector<bool> deleteNode = std::vector<bool>(), 
        std::vector<bool> deleteProcess = std::vector<bool>()
    );

    // Read in the pipeline from a STAR file
    void read(bool do_lock = false, std::string lock_message = "Undefined lock message");

    struct rwlock {

        PipeLine &pipeline;

        rwlock(PipeLine& pipeline, std::string message): pipeline(pipeline) {
            pipeline.read(DO_LOCK, message);
        }

        ~rwlock() {
            pipeline.write(DO_LOCK);
        }

    };

};

class PipeLineFlowChart {

    public:

    // Use short process names, or original, full ones
    bool do_short_names;

    // Also make upwardsFlowCharts for all branches?
    bool do_branches;

    // All the processes for which a upwardFlowChart will be made
    std::vector<long int> todo_list;

    public:

    PipeLineFlowChart() {
        do_branches = true;
        do_short_names = false;
    }

    // Write how many particles or classes or whatever the node is that represents a downward arrow
    std::string getDownwardsArrowLabel(
        PipeLine &pipeline, long int lower_process, long int new_process
    );

    // The process will be added to the top
    // The function returns the parent process from which the upper_node came
    // It will return a negative value if there was no parent process
    long int addProcessToUpwardsFlowChart(
        std::ofstream &fh, PipeLine &pipeline, long int lower_process,
        long int new_process, std::vector<long int> &branched_procs
    );

    void makeOneUpwardsFlowChart(
        std::ofstream &fh, PipeLine &pipeline, long int from_node,
        std::vector<long int> &all_branches, bool is_main_flow
    );

    void makeAllUpwardsFlowCharts(
        FileName &fn_out, PipeLine &pipeline, long int from_process
    );

    // Open and close a new flowchart picture
    void openTikZPicture(std::ofstream &fh, bool is_main_flow);

    void closeTikZPicture(std::ofstream &fh, bool is_main_flow);

    void adaptNamesForTikZ(FileName &name);

    // Open and close a new output file
    void openFlowChartFile(FileName &fn_out, std::ofstream &fh);

    void closeFlowChartFile(std::ofstream &fh);

};

#endif /* PIPELINER_H_ */
