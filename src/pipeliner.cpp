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

#include "src/pipeliner.h"
#include <unistd.h>
#include <unordered_map>

// #define DEBUG

long int PipeLine::addNode(Node &node, bool touch_if_not_exist) {

    if (node.name.empty())
        REPORT_ERROR("PipeLine::addNode ERROR: Adding an empty nodename. Did you fill in all Node names correctly?");

    // Check if node has an aliased name
    FileName fn_node = node.name;
    auto search = std::find_if(
        processList.begin(), processList.end(),
        [&fn_node] (const Process &process) { return fn_node.contains(process.alias); }
    );
    // If so, revert back to the original name
    if (search != processList.end()) {
        node.name = search->name + fn_node.without(search->alias);
    }

    auto it = std::find_if(
        nodeList.begin(), nodeList.end(),
        [&node] (const Node &other_node) { return node.name == other_node.name; }
    );

    if (it == nodeList.end()) {
        nodeList.push_back(node);  // invalidates it
        it = nodeList.end() - 1;
    }

    touchTemporaryNodeFile(*it, touch_if_not_exist);

    return it - nodeList.begin();

}


void PipeLine::addNewInputEdge(Node &node, long int myProcess) {

    // 1. Check whether Node with that name already exists in the Node list
    long int oldsize = nodeList.size();
    long int myNode  = addNode(node);
    long int newsize = nodeList.size();

    // 2. Set the input_for_process in the inputForProcessList of this Node
    // But only if it doesn't exist yet
    auto &inputForProcessList = nodeList[myNode].inputForProcessList;
    if (std::find(
        inputForProcessList.begin(), inputForProcessList.end(), myProcess
    ) == inputForProcessList.end()) {
        inputForProcessList.push_back(myProcess);
        processList[myProcess].inputNodeList.push_back(myNode);
    }

    if (newsize > oldsize) {
        // This is a previously unobserved Node, being used as input in a new Process. Check whether it came from an old Process
        for (long int i = 0; i < processList.size(); i++) {
            /// TODO: check this name comparison!!!
            FileName nodename = (nodeList[myNode]).name;
            if (nodename.contains(processList[i].name)) {
                #ifdef DEBUG
                std::cerr << " nodename.find(processList[i].name) = " << nodename.find(processList[i].name) << " processList[i].name= " << processList[i].name << std::endl;
                std::cerr << " found! " <<  nodename << " as output from " << processList[i].name << std::endl;
                #endif
                processList[i].outputNodeList.push_back(myNode);
                nodeList[myNode].outputFromProcess = i;
                break;
            }
        }
    }
}

void PipeLine::addNewOutputEdge(long int myProcess, Node &node) {

    long int old_size = nodeList.size();

    // 1. Check whether Node with that name already exists in the Node list
    // Touch .Nodes entries even if they don't exist for scheduled jobs
    bool touch_if_not_exist = (processList[myProcess].status == Process::SCHEDULED);

    // 2. Set the output_from_process of this Node
    node.outputFromProcess = myProcess;
    long int myNode = addNode(node, touch_if_not_exist);

    // 3. Only for new Nodes, add this Node to the outputNodeList of myProcess
    if (myNode == old_size)
        processList[myProcess].outputNodeList.push_back(myNode);

}

long int PipeLine::addNewProcess(Process &process, bool do_overwrite) {
    // Is there already a process with the same name in the processList?
    auto it = std::find_if(
        processList.begin(), processList.end(), 
        [&process] (const Process &other_process) { return process.name == other_process.name; }
    );
    if (it == processList.end()) {
        processList.push_back(process);  // invalidates it
        it = processList.end() - 1;
        job_counter++;
    } else {
        it->status = process.status;
        if (!do_overwrite)
        REPORT_ERROR("PipeLine::addNewProcess: ERROR: trying to add existing Process to the pipeline, while overwriting is not allowed.");
    }
    return it - processList.begin();
}

long int PipeLine::findNodeByName(std::string name) {
    for (long int ipos = 0; ipos < nodeList.size(); ipos++) {
        if (nodeList[ipos].name == name) return ipos;
    }
    return -1;
}


long int PipeLine::findProcessByName(std::string name) {
    for (long int ipos = 0; ipos < processList.size(); ipos++) {
        if (processList[ipos].name == name) return ipos;
    }
    return -1;
}

long int PipeLine::findProcessByAlias(std::string name) {
    if (name == "None") return -1;

    for (long int ipos = 0; ipos < processList.size(); ipos++) {
        if (processList[ipos].alias == name) return ipos;
    }
    return -1;
}

bool PipeLine::touchTemporaryNodeFile(Node &node, bool touch_even_if_not_exist) {
    if (do_read_only) return false;

    FileName fn_dir = ".Nodes/";
    FileName fnt;

    // Check whether there is an alias for the corresponding process
    FileName fn_alias = node.outputFromProcess < 0 ? "None" : processList[node.outputFromProcess].alias;

    if (fn_alias != "None") {
        // Make sure fn_alias ends with a slash
        if (fn_alias[fn_alias.length() - 1] != '/')
            fn_alias += "/";
        FileName fn_pre, fn_jobnr, fn_post;
        if (decomposePipelineFileName(node.name, fn_pre, fn_jobnr, fn_post)) {
            fnt = fn_alias + fn_post;
        } else {
            REPORT_ERROR("PipeLine::touchTemporaryNodeFile ERROR: invalid node name: " + node.name);
        }
    } else {
        fnt = node.name;
    }

    if (exists(node.name) || touch_even_if_not_exist) {
        // Make subdirectory for each type of node
        FileName fn_type = integerToString(node.type) + "/";
        FileName mydir = fn_dir + fn_type + fnt.substr(0, fnt.rfind("/") + 1);
        FileName mynode = fn_dir + fn_type + fnt;
        std::string command;
        if (!exists(mydir))
            system(("mkdir -p " + mydir).c_str());
        touch(mynode);
        return true;
    }
    return false;
}

void PipeLine::touchTemporaryNodeFiles(Process &process) {
    if (do_read_only) return;

    bool touch_if_not_exist = process.status == Process::SCHEDULED;
    for (long int node : process.outputNodeList)
        touchTemporaryNodeFile(nodeList[node], touch_if_not_exist);
}

void PipeLine::deleteTemporaryNodeFile(Node &node) {
    if (do_read_only) return;

    FileName fn_dir = ".Nodes/";

    // Check whether there is an alias for the corresponding process
    FileName fn_alias = node.outputFromProcess < 0 ? "None" : processList[node.outputFromProcess].alias;

    FileName fnt;
    if (fn_alias != "None") {
        // Make sure fn_alias ends with a slash
        if (fn_alias[fn_alias.length()-1] != '/')
            fn_alias += "/";
        FileName fn_pre, fn_jobnr, fn_post;
        if (decomposePipelineFileName(node.name, fn_pre, fn_jobnr, fn_post)) {
            fnt = fn_alias + fn_post;
        } else {
            REPORT_ERROR("PipeLine::deleteTemporaryNodeFile ERROR: invalid node name: " + node.name);
        }
    } else {
        fnt = node.name;
    }

    FileName fn_type = integerToString(node.type) + "/";
    FileName fn = fn_dir + fn_type + fnt;
    remove(fn.c_str());

    // Also remove the directory if it is empty
    /// TODO: Check what happens if the directory is not empty yet.
    fn = fn.beforeLastOf("/");
    remove(fn.c_str());

}

void PipeLine::deleteTemporaryNodeFiles(Process &process) {
    if (do_read_only) return;

    for (long int node : process.outputNodeList)
        deleteTemporaryNodeFile(nodeList[node]);
}

void PipeLine::remakeNodeDirectory() {
    if (do_read_only) return;

    // Clear existing directory
    FileName fn_dir = ".Nodes/";
    system((" rm -rf " + fn_dir).c_str());

    for (Node &node : nodeList) {
        int process = node.outputFromProcess;
        bool touch_if_not_exist = process >= 0 && processList[process].status == Process::SCHEDULED;
        touchTemporaryNodeFile(node, touch_if_not_exist);
    }
    system(("chmod 777 -R " + fn_dir).c_str());
}


bool PipeLine::checkProcessCompletion() {
    if (do_read_only) return false;

    std::vector<long int> finished_success_processes;
    std::vector<long int> finished_failure_processes;
    std::vector<long int> finished_aborted_processes;

    for (long int i = 0; i < processList.size(); i++) {
        // Only check running processes for file existence
        if (processList[i].status == Process::RUNNING) {

            if (exists(processList[i].name + RELION_JOB_EXIT_SUCCESS)) {
                finished_success_processes.push_back(i);
            } else if (exists(processList[i].name + RELION_JOB_EXIT_FAILURE)) {
                finished_failure_processes.push_back(i);
            } else if (exists(processList[i].name + RELION_JOB_EXIT_ABORTED)) {
                finished_aborted_processes.push_back(i);
            }
        }
    }

    // Only do read/write cycle in case a process was finished, otherwise the GUI slows down too much
    if (
        finished_success_processes.empty() &&
        finished_failure_processes.empty() &&
        finished_aborted_processes.empty() &&
        !exists(PIPELINE_HAS_CHANGED)
    ) return false;

    // Read in the latest version of the pipeline, just in case anyone else made a change meanwhile...
    std::string lock_message = "";
    if (finished_success_processes.size() > 0) {
        lock_message += "checkProcessCompletion: the following jobs have successfully finished: ";
        for (long int i = 0; i < finished_success_processes.size(); i++)
            lock_message += " " + processList[finished_success_processes[i]].name;
        lock_message += "\n";
    }
    if (finished_failure_processes.size() > 0) {
        lock_message += "checkProcessCompletion: the following jobs have failed with an error: ";
        for (long int i = 0; i < finished_failure_processes.size(); i++)
            lock_message += " " + processList[finished_failure_processes[i]].name;
        lock_message += "\n";
    }
    if (finished_aborted_processes.size() > 0) {
        lock_message += "checkProcessCompletion: the following jobs have been aborted: ";
        for (long int i = 0; i < finished_aborted_processes.size(); i++)
            lock_message += " " + processList[finished_aborted_processes[i]].name;
        lock_message += "\n";
    }

    read(DO_LOCK, lock_message);

    // Set the new status of all the finished processes
    for (int i = 0; i < finished_success_processes.size(); i++) {
        int myproc = finished_success_processes[i];
        processList[myproc].status = Process::FINISHED_SUCCESS;

        // Also see whether there was an output nodes starfile
        getOutputNodesFromStarFile(myproc);

        // Also touch the outputNodes in the .Nodes directory
        for (long int j = 0; j < processList[myproc].outputNodeList.size(); j++) {
            int myNode = (processList[myproc]).outputNodeList[j];
            if (myNode < 0 || myNode >= nodeList.size())
                REPORT_ERROR("pipeline checkProcessCompletion ERROR: " + integerToString(j) + "th output node of " + processList[myproc].name + " is invalid: " + integerToString(myNode));
            if (!touchTemporaryNodeFile(nodeList[myNode]))
                std::cerr << "WARNING: output node " + nodeList[myNode].name + " does not exist, while job " + processList[myproc].name +" should have finished successfully. You can manually mark this job as failed to suppress this message." << std::endl;
        }

    }
    for (long int i : finished_failure_processes) {
        processList[i].status = Process::FINISHED_FAILURE;
    }
    for (long int i : finished_aborted_processes) {
        processList[i].status = Process::FINISHED_ABORTED;
    }

    // Always couple read/write with DO_LOCK
    // This is to make sure two different windows do not get out-of-sync
    write(DO_LOCK);
    if (exists(PIPELINE_HAS_CHANGED)) std::remove(PIPELINE_HAS_CHANGED);
    return true;

}

std::string PipeLine::getCommandLineJob(
    RelionJob &thisjob, int current_job, bool is_main_continue,
    bool is_scheduled, bool do_makedir, bool do_overwrite_current,
    std::vector<std::string> &commands
) throw (std::string) {

    if (do_overwrite_current) { is_main_continue = false; }

    // Except for continuation or scheduled jobs, all jobs get a new unique directory
    std::string my_outputname;
    if ((
        is_main_continue || is_scheduled || do_overwrite_current
    ) && current_job < processList.size()) {
        if (current_job < 0)
            REPORT_ERROR("BUG: current_job < 0");
        FileName fn_settings = processList[current_job].name;
        my_outputname = fn_settings.beforeLastOf("/") + "/";
    } else {
        my_outputname = "";
    }

    // Set is_continue flag inside the job
    thisjob.is_continue = is_main_continue;

    std::string final_command = thisjob.getCommands(my_outputname, commands, do_makedir, job_counter);
    if (commands.empty())
    throw " PipeLine::getCommandLineJob: Nothing to do...";
    return final_command;
}

// Adds thisjob to the pipeline and returns the id of the newprocess
long int PipeLine::addJob(
    RelionJob &thisjob, int as_status, bool do_overwrite, bool do_write_minipipeline
) {

    // Also write a mini-pipeline in the output directory
    PipeLine mini_pipeline;
    mini_pipeline.setName(thisjob.outputName + "job");

    // Add Process to the processList of the pipeline
    Process process(thisjob.outputName, thisjob.type, as_status);
    long int myProcess = addNewProcess(process, do_overwrite);
    mini_pipeline.addNewProcess(process);

    // Add all input nodes
    for (auto &node : thisjob.inputNodes) {
        addNewInputEdge(node, myProcess);
        mini_pipeline.addNewInputEdge(node, 0);
    }
    // Add all output nodes
    for (auto &node : thisjob.outputNodes) {
        addNewOutputEdge(myProcess, node);
        mini_pipeline.addNewOutputEdge(0, node);
    }

    if (do_write_minipipeline) {
        // Write the mini-pipeline to an updated STAR file
        mini_pipeline.write();
    }
    // Writing of the overall pipeline is done in the function calling addToPipeLine

    return myProcess;
}

void PipeLine::runJob(
    RelionJob &_job, int &current_job, bool only_schedule, bool is_main_continue,
    bool is_scheduled, bool do_overwrite_current
) throw (std::string) {

    // Remove run.out and run.err when overwriting a job
    if (do_overwrite_current) { is_main_continue = false; }

    std::vector<std::string> commands;
    std::string final_command = getCommandLineJob(
        _job, current_job, is_main_continue, is_scheduled, true, // makedir
        do_overwrite_current, commands
    );

    // Remove run.out and run.err when overwriting a job
    if (do_overwrite_current) {
        // Completely empty the output directory, NOTE that  _job.outputName+ is not defined until AFTER calling getCommandLineJob!!!
        std::string command = " rm -rf " + _job.outputName + "*";
        int res = system(command.c_str());

        // Above deletes run_submit.script too, so we have to call this again ...
        final_command = getCommandLineJob(
            _job, current_job, is_main_continue, is_scheduled, true,
            do_overwrite_current, commands
        );
    }

    // Read in the latest version of the pipeline, just in case anyone else made a change meanwhile...
    std::string lock_message = "runJob: " + _job.outputName;
    read(DO_LOCK, lock_message);

    // Save temporary hidden file with this jobs settings as default for a new job
    _job.write("");

    // Also save a copy of the GUI settings with the current output name
    _job.write(_job.outputName);

    // Make sure none of the exit or abort files from the pipeline_control system are here from before
    std::remove((_job.outputName + RELION_JOB_ABORT_NOW).c_str());
    std::remove((_job.outputName + RELION_JOB_EXIT_ABORTED).c_str());
    std::remove((_job.outputName + RELION_JOB_EXIT_SUCCESS).c_str());
    std::remove((_job.outputName + RELION_JOB_EXIT_FAILURE).c_str());


    /*
    // If this is a continuation job, check whether output files exist and move away!
    // This is to ensure that the continuation job goes OK and will show up as 'running' in the GUI
    bool do_move_output_nodes_to_old = false;
    if (!only_schedule && is_main_continue)
    {
        do_move_output_nodes_to_old = !(processList[current_job].type == Process::CLASS2D ||
                                        processList[current_job].type == Process::CLASS3D ||
                                        processList[current_job].type == Process::INIMODEL ||
                                        processList[current_job].type == Process::AUTO3D ||
                                        processList[current_job].type == Process::MULTIBODY ||
                                        processList[current_job].type == Process::MANUALPICK ||
                                        processList[current_job].type == Process::CLASSSELECT);

        // For continuation of relion_refine jobs, remove the original output nodes from the list
        if (processList[current_job].type == Process::CLASS2D ||
            processList[current_job].type == Process::CLASS3D ||
            processList[current_job].type == Process::AUTO3D ||
            processList[current_job].type == Process::MULTIBODY ||
            processList[current_job].type == Process::INIMODEL)
        {

            std::vector<bool> deleteNodes, deleteProcesses;
            deleteNodes.resize(nodeList.size(), false);
            deleteProcesses.resize(processList.size(), false);

            for (long int inode = 0; inode < (processList[current_job]).outputNodeList.size(); inode++)
            {
                long int mynode = (processList[current_job]).outputNodeList[inode];
                if (!exists(nodeList[mynode].name))
                    deleteNodes[mynode] = true;
            }

            FileName fn_del = "tmp";
            write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
            std::remove("tmpdeleted_pipeline.star");

            // Read the updated pipeline back in again
            lock_message += " part 2";
            read(DO_LOCK, lock_message);

        }
    } // end if !only_schedule && is_main_continue

    // If a job is executed with a non-continue scheduled job, then also move away any existing output node files
    if (current_job >= 0 && (is_scheduled && !is_main_continue) || do_overwrite_current)
        do_move_output_nodes_to_old = true;

    // Move away existing output nodes
    if (do_move_output_nodes_to_old)
    {

        for (int i = 0; i < processList[current_job].outputNodeList.size(); i++)
        {
            int j = processList[current_job].outputNodeList[i];
            std::string fn_node = nodeList[j].name;
            if (exists(fn_node))
            {
                std::string path2 =  fn_node + ".old";
                rename(fn_node.c_str(), path2.c_str());
            }
        }
    }
    */

    // For continuation of relion_refine jobs, remove the original output nodes from the list
    if (!only_schedule && is_main_continue) {
        if (
            processList[current_job].type == Process::CLASS2D ||
            processList[current_job].type == Process::CLASS3D ||
            processList[current_job].type == Process::AUTO3D ||
            processList[current_job].type == Process::MULTIBODY ||
            processList[current_job].type == Process::INIMODEL
        ) {

            std::vector<bool> deleteNodes, deleteProcesses;
            deleteNodes.resize(nodeList.size(), false);
            deleteProcesses.resize(processList.size(), false);

            for (long int node : processList[current_job].outputNodeList) {
                if (!exists(nodeList[node].name))
                    deleteNodes[node] = true;
            }

            FileName fn_del = "tmp";
            write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
            std::remove("tmpdeleted_pipeline.star");

            // Read the updated pipeline back in again
            lock_message += " part 2";
            read(DO_LOCK, lock_message);

        }
    }

    // Now save the job (and its status) to the PipeLine
    int mynewstatus = only_schedule ? Process::SCHEDULED : Process::RUNNING;
    bool allow_overwrite = is_main_continue || is_scheduled; // continuation and scheduled jobs always allow overwriting into the existing directory

    // Add the job to the pipeline, and set current_job to the new one
    current_job = addJob(_job, mynewstatus, allow_overwrite || do_overwrite_current);

    // Write out the new pipeline
    write(DO_LOCK);

    // Now actually execute the Job
    if (!only_schedule) {
        //std::cout << "Executing: " << final_command << std::endl;
        int res = system(final_command.c_str());

        // Also print the final_command to the note for future reference
        FileName fn_note = processList[current_job].name + "note.txt";
        std::ofstream ofs;
        ofs.open (fn_note.c_str(), std::ofstream::out | std::ofstream::app);

        // current date/time based on current system
        time_t now = time(0);
        ofs << std::endl << " ++++ Executing new job on " << ctime(&now);
        ofs <<  " ++++ with the following command(s): " << std::endl;
        for (size_t icom = 0; icom < commands.size(); icom++)
            ofs << commands[icom] << std::endl;
        ofs <<  " ++++ " << std::endl;
        ofs.close();
    }

    // Copy pipeline star file as backup to the output directory
    FileName fn_pipe = name + "_pipeline.star";
    if (exists(fn_pipe))
    copy(fn_pipe, processList[current_job].name + fn_pipe);

}

static const std::unordered_map<std::string, int> string2type_mapping {
    {{Process::IMPORT_NAME},      Process::IMPORT},
    {{Process::MOTIONCORR_NAME},  Process::MOTIONCORR},
    {{Process::CTFFIND_NAME},     Process::CTFFIND},
    {{Process::MANUALPICK_NAME},  Process::MANUALPICK},
    {{Process::AUTOPICK_NAME},    Process::AUTOPICK},
    {{Process::EXTRACT_NAME},     Process::EXTRACT},
    {{Process::CLASSSELECT_NAME}, Process::CLASSSELECT},
    {{Process::CLASS2D_NAME},     Process::CLASS2D},
    {{Process::CLASS3D_NAME},     Process::CLASS3D},
    {{Process::AUTO3D_NAME},      Process::AUTO3D},
    {{Process::MASKCREATE_NAME},  Process::MASKCREATE},
    {{Process::JOINSTAR_NAME},    Process::JOINSTAR},
    {{Process::SUBTRACT_NAME},    Process::SUBTRACT},
    {{Process::POST_NAME},        Process::POST},
    {{Process::RESMAP_NAME},      Process::RESMAP},
    {{Process::INIMODEL_NAME},    Process::INIMODEL},
    {{Process::EXTERNAL_NAME},    Process::EXTERNAL}
};

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(std::string typestring, std::string fn_options) {
    auto it = string2type_mapping.find(typestring);
    if (it == string2type_mapping.end())
        REPORT_ERROR("ERROR: unrecognised string for job type: " + typestring);

    return addScheduledJob(it->second, fn_options);
}

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(int job_type, std::string fn_options) {
    RelionJob job (job_type);
    std::vector<std::string> options = split(fn_options, ";");
    for (const auto &option : options)
        job.setOption(option);

    // Always add preprocessing jobs as continuation ones (for convenient on-the-fly processing)
    if (
        job_type == Process::MOTIONCORR || job_type == Process::CTFFIND ||
        job_type == Process::AUTOPICK   || job_type == Process::EXTRACT
    ) { job.is_continue = true; }

    try {
        int current_job = processList.size();
        runJob(job, current_job, true, job.is_continue, false, false);
        // true is only_schedule, false means !is_scheduled, 2nd false means dont overwrite current
        return current_job;
    } catch (const std::string &errmsg) {
        REPORT_ERROR(errmsg.c_str());
    }
}

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(RelionJob &job, std::string fn_options) {

    std::vector<std::string> options = split(fn_options, ";");
    for (const auto &option : options)
        job.setOption(option);

    try {
        int current_job = processList.size();
        runJob(job, current_job, true, job.is_continue, false, false);
        //         only_schedule ^        !is_scheduled ^      ^ !do_overwrite_current
        return current_job;
    } catch (const std::string &errmsg) {
        REPORT_ERROR(errmsg.c_str());
    }
}

void PipeLine::waitForJobToFinish(
    int current_job, bool &is_failure, bool &is_aborted
) {
    while (true) {
        sleep(10);
        checkProcessCompletion();
        if (
            processList[current_job].status == Process::FINISHED_SUCCESS ||
            processList[current_job].status == Process::FINISHED_ABORTED ||
            processList[current_job].status == Process::FINISHED_FAILURE
        ) {
            // Prepare a string for a more informative .lock file
            std::string lock_message = " pipeliner noticed that " + processList[current_job].name + " finished and is trying to update the pipeline";

            // Read in existing pipeline, in case some other window had changed something else
            read(DO_LOCK, lock_message);

                   if (processList[current_job].status == Process::FINISHED_FAILURE) {
                is_failure = true;
            } else if (processList[current_job].status == Process::FINISHED_ABORTED) {
                is_aborted = true;
            }

            // Write out the modified pipeline with the new status of current_job
            write(DO_LOCK);
            break;
        }
    }
}

void PipeLine::runScheduledJobs(
    FileName fn_sched, FileName fn_jobids, int nr_repeat,
    long int minutes_wait, long int minutes_wait_before, long int seconds_wait_after,
    bool do_overwrite_current
) {

    if (fn_jobids.empty())
        REPORT_ERROR("PipeLine::runScheduledJobs: Nothing to do...");

    std::vector<std::string> jobids = split(fn_jobids, " ");
    std::vector<FileName> my_scheduled_processes;
    for (const std::string &id : jobids)
        my_scheduled_processes.push_back(id);

    FileName fn_log = "pipeline_" + fn_sched + ".log";
    std::ofstream fh (fn_log.c_str(), std::ios::app);

    if (nr_repeat > 1) {
        std::cout << " PIPELINER: writing out information in logfile " << fn_log << std::endl;
    }

    // Touch the fn_check file
    FileName fn_check = "RUNNING_PIPELINER_" + fn_sched;
    bool fn_check_exists = false;
    if (nr_repeat > 1) {
        touch(fn_check);
        fn_check_exists = true;
    }

    fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    fh << " Starting a new scheduler execution called " << fn_sched << std::endl;
    fh << " The scheduled jobs are: " << std::endl;
    for (long int i = 0; i < my_scheduled_processes.size(); i++)
        fh << " - " << my_scheduled_processes[i] << std::endl;
    if (nr_repeat > 1) {
        if (minutes_wait_before > 0)
            fh << " Will wait " << minutes_wait_before << " minutes before running the first job." << std::endl;
        fh << " Will execute the scheduled jobs " << nr_repeat << " times." << std::endl;
        if (nr_repeat > 1)
            fh << " Will wait until at least " << minutes_wait << " minutes have passed between each repeat." << std::endl;
        fh << " Will be checking for existence of file " << fn_check << "; if it no longer exists, the scheduler will stop." << std::endl;
    }
    fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    // Wait this many minutes before starting the repeat cycle...
    if (minutes_wait_before > 0)
        sleep(minutes_wait_before * 60);

    bool is_failure = false;
    bool is_aborted = false;

    int repeat = 0;
    time_t now = time(0);
    for (repeat = 0 ; repeat < nr_repeat; repeat++) {
        if (nr_repeat > 1) {
            fh << " + " << ctime(&now) << " -- Starting the " << repeat + 1 << "th repeat" << std::endl;
        }

        // Get starting time of the repeat cycle
        timeval time_start, time_end;
        gettimeofday(&time_start, nullptr);

        for (long int i = 0; i < my_scheduled_processes.size(); i++) {
            int current_job = findProcessByName(my_scheduled_processes[i]);
            if (current_job < 0) {
                // Also try finding it by alias
                current_job = findProcessByAlias(my_scheduled_processes[i]);
                if (current_job < 0)
                    REPORT_ERROR("ERROR: cannot find process with name: " + my_scheduled_processes[i]);
            }
            RelionJob myjob;
            bool is_continue;
            if (!myjob.read(processList[current_job].name, is_continue, true))
                // true means also initialise the job
                REPORT_ERROR("There was an error reading job: " + processList[current_job].name);

            // Check whether the input nodes are there, before executing the job
            for (long int inode = 0; inode < processList[current_job].inputNodeList.size(); inode++) {
                long int mynode = processList[current_job].inputNodeList[inode];
                while (!exists(nodeList[mynode].name)) {
                    fh << " + -- Warning " << nodeList[mynode].name << " does not exist. Waiting 60 seconds ... " << std::endl;
                    sleep(60);
                }
            }
            now = time(0);
            fh << " + " << ctime(&now) << " ---- Executing " << processList[current_job].name  << std::endl;

            try {
                runJob(myjob, current_job, false, is_continue, true, do_overwrite_current);
                // true means is_scheduled; false=dont overwrite current
            } catch (const std::string &errmsg) {
                REPORT_ERROR(errmsg);
            }

            // Now wait until that job is done!
            while (true) {
                if (nr_repeat > 1 && !exists(fn_check)) {
                    fn_check_exists = false;
                    break;
                }

                sleep(seconds_wait_after);
                checkProcessCompletion();
                if (
                    processList[current_job].status == Process::FINISHED_SUCCESS ||
                    processList[current_job].status == Process::FINISHED_ABORTED ||
                    processList[current_job].status == Process::FINISHED_FAILURE
                ) {
                    // Prepare a string for a more informative .lock file
                    std::string lock_message =
                        " Scheduler " + fn_sched + " noticed that "
                        + processList[current_job].name + " finished and is trying to update the pipeline";

                    // Read in existing pipeline, in case some other window had changed something else
                    read(DO_LOCK, lock_message);

                    if (processList[current_job].status == Process::FINISHED_SUCCESS) {
                        // Will we go on to do another repeat?
                        if (repeat + 1 != nr_repeat) {
                            int mytype = processList[current_job].type;
                            // The following jobtypes have functionality to only do the unfinished part of the job
                            if (
                                mytype == Process::MOTIONCORR || mytype == Process::CTFFIND ||
                                mytype == Process::AUTOPICK   || mytype == Process::EXTRACT ||
                                mytype == Process::CLASSSELECT
                            ) {
                                myjob.is_continue = true;
                                // Write the job again, now with the updated is_continue status
                                myjob.write(processList[current_job].name);
                            }
                            processList[current_job].status = Process::SCHEDULED;
                        } else {
                            processList[current_job].status = Process::FINISHED_SUCCESS;
                        }
                    } else if (processList[current_job].status == Process::FINISHED_FAILURE) {
                        is_failure = true;
                    } else if (processList[current_job].status == Process::FINISHED_ABORTED) {
                        is_aborted = true;
                    }

                    // Write out the modified pipeline with the new status of current_job
                    write(DO_LOCK);
                    break;
                }
            }

            // break out of scheduled processes loop
            if (is_failure || is_aborted) break;

            if (nr_repeat > 1 && !fn_check_exists) break;

        }

        // break out of repeat loop
        if (nr_repeat > 1 && !fn_check_exists) break;

        if (is_failure || is_aborted) break;

        // Wait at least until 'minutes_wait' minutes have passed from the beginning of the repeat cycle
        gettimeofday(&time_end, nullptr);
        long int passed_minutes = (time_end.tv_sec - time_start.tv_sec) / 60;
        long int still_wait = minutes_wait - passed_minutes;
        if (still_wait > 0 && repeat + 1 != nr_repeat) {
            fh << " + -- Waiting " << still_wait << " minutes until next repeat .."<< std::endl;
            sleep(still_wait * 60);
        }
    }

    if (is_failure) {
        fh << " + Stopping pipeliner due to a job that failed with an error ..." << std::endl;
    } else if (is_aborted) {
        fh << " + Stopping pipeliner due to user abort .. " << std::endl;
    } else if (repeat == nr_repeat) {
        if (nr_repeat > 1) {
            fh << " + performed all requested repeats in scheduler " << fn_sched << ". Stopping pipeliner now ..." << std::endl;
        } else {
            fh << " + All jobs have finished, so stopping pipeliner now ..." << std::endl;
        }

        // Read in existing pipeline, in case some other window had changed it
        std::string lock_message = " Scheduler " + fn_sched + " has finished and is trying to update the pipeline";
        read(DO_LOCK, lock_message);

        // After breaking out of repeat, set status of the jobs to finished
        for (long int i = 0; i < my_scheduled_processes.size(); i++) {
            int current_job = findProcessByName(my_scheduled_processes[i]);
            processList[current_job].status = Process::FINISHED_SUCCESS;
        }

        // Write the pipeline to an updated STAR file
        write(DO_LOCK);

        // Remove the temporary file
        std::remove(fn_check.c_str());
    } else if (!fn_check_exists && nr_repeat > 1) {
        fh << " + File " << fn_check << " was removed. Stopping now .." << std::endl;
        std::cout << " PIPELINER: the " << fn_check << " file was removed. Stopping now ..." << std::endl;
    }

    fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

}

void PipeLine::deleteJobGetNodesAndProcesses(
    int this_job, bool do_recursive,
    std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses
) {

    deleteProcesses.resize(processList.size(), false);
    deleteNodes    .resize(nodeList.size(),    false);

    std::vector<long int> to_delete_processes {this_job};

    bool is_done = false;
    size_t istart = 0;
    while (!is_done) {
        size_t imax = to_delete_processes.size();
        for (long int i = istart; i < imax; i++) {
            // re-set istart for next recursive round
            istart = imax;
            long int idel = to_delete_processes[i];
            deleteProcesses[idel] = true;
            is_done = true;
            for (long int node : processList[idel].outputNodeList) {
                deleteNodes[node] = true;

                if (do_recursive) {
                    // Check whether this node is being used as input for another process, and if so, delete those as well
                    for (long int iproc : nodeList[node].inputForProcessList) {
                        // See if this process is not already in the list to be deleted
                        if (std::find(
                            to_delete_processes.begin(), to_delete_processes.end(), iproc
                        ) == to_delete_processes.end()) {
                            to_delete_processes.push_back(iproc);
                            is_done = false;
                        }
                    }
                }
            }
        }
    }
}

void PipeLine::deleteNodesAndProcesses(
    std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses
) {

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "deleteNodesAndProcesses";
    read(DO_LOCK, lock_message);

    // Write new pipeline without the deleted processes and nodes to disc and read in again
    auto it = std::find(deleteProcesses.begin(), deleteProcesses.end(), true);
    FileName fn_del = it == deleteProcesses.end() ? "" :
        processList[it - deleteProcesses.begin()].name;

    write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);

    // Delete the output directories for all selected processes from the hard disk
    // Do this after pipeline.write to get the deleted_pipeline.star still in the correct directory
    for (int i = 0; i < deleteProcesses.size(); i++) {
        if (deleteProcesses[i]) {
            FileName alldirs = ((FileName) processList[i].name).beforeLastOf("/");
            // Move entire output directory (with subdirectory structure) to the Trash folder
            FileName firstdirs = alldirs.beforeLastOf("/");
            FileName fn_tree = "Trash/" + firstdirs;
            mktree(fn_tree);
            std::string command = "mv -f " + alldirs + " " + "Trash/" + firstdirs + "/.";
            system(command.c_str());
            // Also remove the symlink if it exists
            FileName fn_alias = processList[i].alias;
            if (fn_alias != "None") {
                unlink(fn_alias.beforeLastOf("/").c_str());
            }

            deleteTemporaryNodeFiles(processList[i]);
        }
    }

    // Read new pipeline back in again
    read(DO_LOCK, lock_message + " part 2");
    write(DO_LOCK);

}

void PipeLine::getOutputNodesFromStarFile(int this_job) {

    // See if a STAR file called RELION_OUTPUT_NODES.star exists, and if so, read in which output nodes were created
    FileName outnodes = processList[this_job].name + "RELION_OUTPUT_NODES.star";
    if (exists(outnodes)) {

        MetaDataTable MDnodes;
        MDnodes.read(outnodes, "output_nodes");

        for (long int i : MDnodes) {
            FileName nodename = MDnodes.getValue<std::string>(EMDL::PIPELINE_NODE_NAME, i);
            int      nodetype = MDnodes.getValue<int>     (EMDL::PIPELINE_NODE_TYPE, i);

            // if this node does not exist yet, then add it to the pipeline
            if (findNodeByName(nodename) < 0) {
                Node node(nodename, nodetype);
                addNewOutputEdge(this_job, node);
            }
        }
    }
}

void PipeLine::markAsFinishedJob(int this_job, bool is_failed) throw (std::string) {

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "markAsFinishedJob";
    read(DO_LOCK, lock_message);

    processList[this_job].status = is_failed ? Process::FINISHED_FAILURE : Process::FINISHED_SUCCESS;

    // For relion_refine jobs, add last iteration optimiser.star, data.star, model.star and class???.mrc to the pipeline
    if (
        processList[this_job].type == Process::CLASS2D ||
        processList[this_job].type == Process::CLASS3D ||
        processList[this_job].type == Process::AUTO3D  ||
        processList[this_job].type == Process::INIMODEL
    ) {
        // Get the last iteration optimiser file
        FileName fn_opt;
        FileName fn_root1 = processList[this_job].alias_or_name();

        std::vector<FileName> fn_opts;
        fn_opt = fn_root1 + "run_it*optimiser.star";
        fn_opt.globFiles(fn_opts);
        // It could also be a continuation
        fn_opt = fn_root1 + "run_ct?_it???_optimiser.star";
        fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
        // It could also be a continuation
        fn_opt = fn_root1 + "run_ct??_it???_optimiser.star";
        fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
        if (fn_opts.size() > 0) {

            fn_opt = fn_opts[fn_opts.size() - 1]; // the last one

            // Also get data.star
            FileName fn_data = fn_opt.without("_optimiser.star") + "_data.star";
            Node node2(fn_data, Node::PART_DATA);
            addNewOutputEdge(this_job, node2);

            FileName fn_root = fn_opt.without("_optimiser.star");
            if (processList[this_job].type == Process::AUTO3D)
                fn_root += "_half1";

            FileName fn_model = fn_root + "_model.star";
            Node node3(fn_model, Node::MODEL);
            addNewOutputEdge(this_job, node3);

            FileName fn_map = fn_root + "_class???.mrc";
            std::vector<FileName> fn_maps;
            fn_map.globFiles(fn_maps);
            for (int i = 0; i < fn_maps.size(); i++) {
                Node node4(fn_maps[i], Node::REF3D);
                addNewOutputEdge(this_job, node4);
            }
        } else {
            processList[this_job].status = Process::RUNNING;
            write(DO_LOCK);
            throw
            "You are trying to mark a relion_refine job as finished that hasn't even started. \n"
            " This will be ignored. Perhaps you wanted to delete it instead?";
        }
    }

    // Also see whether there is an output nodes star file...
    getOutputNodesFromStarFile(this_job);

    // Remove any of the expected output nodes from the pipeline if the corresponding file doesn't already exist
    std::vector<bool> deleteNodes, deleteProcesses;
    deleteNodes.resize(nodeList.size(), false);
    deleteProcesses.resize(processList.size(), false);

    for (long int inode = 0; inode < (processList[this_job]).outputNodeList.size(); inode++) {
        long int mynode = (processList[this_job]).outputNodeList[inode];
        if (!exists(nodeList[mynode].name))
            deleteNodes[mynode] = true;
    }
    FileName fn_del = "tmp";
    write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
    std::remove("tmpdeleted_pipeline.star");

    // Read the updated pipeline back in and write it out again
    // With the locking mechanism, each pipeline.read(bool, DO_LOCK) needs to be followed soon by a pipeline.write(DO_LOCK)!
    lock_message += " part 2";
    read(DO_LOCK, lock_message);
    write(DO_LOCK);

}

// Set the alias for a job
void PipeLine::setAliasJob(int this_job, std::string alias) throw (std::string) {

    Process &process = processList[this_job];

    FileName fn_pre, fn_jobnr, fn_post;
    if (!decomposePipelineFileName(
        process.name, fn_pre, fn_jobnr, fn_post
    ))
        REPORT_ERROR("PipeLine::setAlias ERROR: invalid pipeline process name: " + process.name);

    const std::vector<char> forbidden_chars {
        '(', ')', '{', '}', '<', '>', '\"', '/', '\\', 
        '&', '|', '*', '?', '#', '%', '$'
    };

    if (alias.empty()) {
        alias = "None";
    } else if (alias.length() < 2) {
        throw "Alias cannot be less than 2 characters, please provide another one";
    } else if (
        alias.substr(0, 3) == "job"
    ) {
        throw "Alias cannot start with 'job', please provide another one";
    } else if (std::find_if(
        forbidden_chars.begin(), forbidden_chars.end(),
        [&alias] (char c) { return alias.find(c) != std::string::npos; }
    ) != forbidden_chars.end()) {
        throw "Alias cannot contain following symbols: " 
                + join(forbidden_chars, ", ") + '.';
    } else {

        // Turn spaces in alias into underscores
        for (auto &c : alias) {
            if (c == ' ') { c = '_'; }
        }

        // Make sure the alias ends with a slash
        if (alias[alias.length() - 1] != '/')
            alias += "/";

        if (alias.length() < 1 || std::find_if(
            processList.begin(), processList.end(),
            // Check uniqueness of the alias
            [&fn_pre, &alias] (const Process &process) {
                return process.alias == fn_pre + alias && alias != "None";  // alias will never be "None" though
            }
        ) != processList.end()) {
            throw "Alias is not unique, please provide another one";
        }
    }

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "setAliasJob";
    read(DO_LOCK, lock_message);

    // Remove the original .Nodes entry
    deleteTemporaryNodeFiles(process);

    FileName fn_old_alias = process.alias;
    if (fn_old_alias != "None")
        fn_old_alias = fn_old_alias.beforeLastOf("/");

    // No alias if the alias contains a unique jobnr string because of continuation of relion_refine jobs
    if (alias == "None" ) {
        process.alias = "None";
    } else {
        // If this was already an alias: remove the old symbolic link
        if (fn_old_alias != "None")
            unlink(fn_old_alias.c_str());

        // Set the alias in the pipeline
        process.alias = fn_pre + alias;

        //Make the new symbolic link
        FileName path1 = "../" + process.name;
        FileName path2 = process.alias;
        symlink(path1.c_str(), path2.beforeLastOf("/").c_str());

    }

    // Remake the new .Nodes entry
    touchTemporaryNodeFiles(process);

    // Write new pipeline to disc
    write(DO_LOCK);

}

void PipeLine::makeFlowChart(long int current_job, bool do_display_pdf) throw (std::string) {

    if (current_job < 0)
    throw " You can only make flowcharts for existing jobs ... ";

    const char *pdf_viewer_exe = getenv("RELION_PDFVIEWER_EXECUTABLE");
    if (!pdf_viewer_exe) { pdf_viewer_exe = DEFAULT::PDFVIEWER; }
    std::string myviewer(pdf_viewer_exe);

    PipeLineFlowChart flowchart;
    FileName fn_dir = processList[current_job].name;
    FileName fn_out = "flowchart.tex";
    flowchart.makeAllUpwardsFlowCharts(fn_out, *this, current_job);
    std::string command = "latex flowchart.tex > flowchart.log && dvipdf flowchart >>flowchart.log && mv flowchart* " + fn_dir;
    std:: cout << " Executing: " << command << std::endl;
    int res = std::system(command.c_str());
    if (do_display_pdf) {
        command = myviewer + " " + fn_dir + "flowchart.pdf &";
        res = std::system(command.c_str());
    }

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "makeFlowChart";
    read(DO_LOCK, lock_message);

    // Add the PDF file as a logfile to the outputnodes of this job, so it can be visualised from the Display button
    Node node(fn_dir + "flowchart.pdf", Node::PDF_LOGFILE);
    addNewOutputEdge(current_job, node);

    write(DO_LOCK);

}

// Undelete a Job from the pipeline, move back from Trash and insert back into the graph
void PipeLine::undeleteJob(FileName fn_undel) {

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "undeleteJob";
    read(DO_LOCK, lock_message);

    importPipeline(fn_undel.beforeLastOf("_pipeline.star"));

    // Copy all processes in the STAR file back into the ProjectDirectory
    MetaDataTable MDproc;
    MDproc.read(fn_undel, "pipeline_processes");
    std::cout <<"  Undeleting from Trash ... " << std::endl;
    for (long int i : MDproc) {
        FileName fn_proc = MDproc.getValue<std::string>(EMDL::PIPELINE_PROCESS_NAME, i);

        // Copy the job back from the Trash folder
        FileName fn_dest = fn_proc.beforeLastOf("/"); //gets rid of ending "/"
        FileName fn_dir_dest = fn_dest.beforeLastOf("/"); // Now only get the job-type directory
        if (!exists(fn_dir_dest)) {
            mktree(fn_dir_dest);
        }
        std::string command = "mv Trash/" + fn_dest + " " + fn_dest;
        std::cout << command << std::endl;
        int res = system(command.c_str());

        // Also re-make all entries in the .Nodes directory
        long int myproc = findProcessByName(fn_proc);
        touchTemporaryNodeFiles(processList[myproc]);
    }
    std::cout << " Done undeleting! " << std::endl;

    // Write the new pipeline to disk and reread it back in again
    write(DO_LOCK);
}

void PipeLine::cleanupJob(int this_job, bool do_harsh) throw (std::string) {
    const auto this_process = processList[this_job];
    std::cout << "Cleaning up " << this_process.name << std::endl;
    if (this_job < 0 || this_process.status != Process::FINISHED_SUCCESS) {
        throw " You can only clean up finished jobs ... ";
    }

    // These job types do not have cleanup:
    if (
        this_process.type == Process::IMPORT ||
        this_process.type == Process::MANUALPICK ||
        this_process.type == Process::CLASSSELECT ||
        this_process.type == Process::MASKCREATE ||
        this_process.type == Process::JOINSTAR ||
        this_process.type == Process::RESMAP
    ) return;

    // Find any subdirectories
    std::vector<FileName> fns_subdir;
    FileName fn_curr_dir = "";
    // Recursively find all subdirectories
    for (int idir = -1, istop = 0; idir < istop; istop = fns_subdir.size(), idir++) {
        FileName fn_curr_dir = idir == -1 ? this_process.name : this_process.name + fns_subdir[idir];
        DIR *dir = opendir(fn_curr_dir.c_str());
        for (dirent *entry = readdir(dir); entry; entry = readdir(dir)) {
            // Only want directories, and not '.' or '..'
            if (entry->d_type == DT_DIR && entry->d_name[0] != '.') {
                FileName fnt = idir == -1 ? entry->d_name : fns_subdir[idir] + entry->d_name;
                fns_subdir.push_back(fnt + "/");
            }
        }
        closedir(dir);
    }

    std::vector<FileName> fns_del;
    FileName fn_pattern;

    // In all jobs cleanup the .old files (from continuation runs)
    //fn_pattern = this_process.name + "*.old";
    //fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del

    ////////// Now see which jobs needs cleaning up
    if (this_process.type == Process::MOTIONCORR) {
        for (const FileName &fn_subdir : fns_subdir) {
            if (do_harsh) {
                // Remove entire movies directory
                fns_del.push_back(this_process.name + fn_subdir);
            } else {
                fn_pattern = this_process.name + fn_subdir + "*.com";
                fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
                fn_pattern = this_process.name + fn_subdir + "*.err";
                fn_pattern.globFiles(fns_del, false);
                fn_pattern = this_process.name + fn_subdir + "*.out";
                fn_pattern.globFiles(fns_del, false);
                fn_pattern = this_process.name + fn_subdir + "*.log";
                fn_pattern.globFiles(fns_del, false);
            }
        }
    } else if (this_process.type == Process::CTFFIND) {
        fn_pattern = this_process.name + "gctf*.out";
        fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
        fn_pattern = this_process.name + "gctf*.err";
        fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
        for (const FileName &fn_subdir : fns_subdir) {
            // Remove entire Micrographs directory structure
            fns_del.push_back(this_process.name + fn_subdir);
        }

    } else if (this_process.type == Process::AUTOPICK) {
        for (const FileName &fn_subdir : fns_subdir) {
            // Remove the Spider files with the FOM maps
            fn_pattern = this_process.name + fn_subdir + "*.spi";
            fn_pattern.globFiles(fns_del, false);
        }
    } else if (this_process.type == Process::EXTRACT) {
        for (const FileName &fn_subdir : fns_subdir) {
            if (do_harsh) {
                // Remove entire directory (STAR files and particle stacks!)
                fns_del.push_back(this_process.name + fn_subdir);
            } else {
                // Only remove the STAR files with the metadata (this will only give moderate file savings)
                fn_pattern = this_process.name + fn_subdir + "*_extract.star";
                fn_pattern.globFiles(fns_del, false);
            }
        }
    } else if (
        this_process.type == Process::CLASS2D  ||
        this_process.type == Process::CLASS3D  ||
        this_process.type == Process::AUTO3D   ||
        this_process.type == Process::INIMODEL ||
        this_process.type == Process::MULTIBODY
    ) {
        // First find the _data.star from each iteration
        std::vector<FileName> fns_iter;
        fn_pattern = this_process.name + "run_it[0-9][0-9][0-9]_data.star";
        fn_pattern.globFiles(fns_iter);
        fn_pattern = this_process.name + "run_ct[0-9]_it[0-9][0-9][0-9]_data.star";
        fn_pattern.globFiles(fns_iter, false);
        fn_pattern = this_process.name + "run_ct[0-9][0-9]_it[0-9][0-9][0-9]_data.star";
        fn_pattern.globFiles(fns_iter, false);
        fn_pattern = this_process.name + "run_ct[0-9][0-9][0-9]_it[0-9][0-9][0-9]_data.star";
        fn_pattern.globFiles(fns_iter, false);
        // Keep everything for the last iteration, such thatone could for example still do a multibody refinement after gentle cleaning
        // Loop over ifile (i.e. the _data.star files from all iterations)
        for (int ifile = 0; ifile < (signed int)(fns_iter.size())-1; ifile++) {
            FileName fn_file = fns_iter[ifile].without("_data.star");
            // Find the iterations to keep: i.e. those that are part of the pipeline
            bool is_in_pipeline = false;
            for (const auto &node : nodeList) {
                if (((FileName) node.name).contains(fn_file)) {
                    is_in_pipeline = true;
                    break;
                }
            }
            // Delete all files from this iteration
            if (!is_in_pipeline) {
                fn_pattern = fn_file + "*";
                fn_pattern.globFiles(fns_del, false);
            }

            // Also clean up maps for PCA movies when doing harsh cleaning
            if (do_harsh && this_process.type == Process::MULTIBODY) {
                fn_pattern = this_process.name + "analyse_component???_bin???.mrc";
                fn_pattern.globFiles(fns_del, false);
            }
        }
    } else if (this_process.type == Process::CTFREFINE) {
        for (const FileName &fn_subdir : fns_subdir) {
            // remove the temporary output files
            fn_pattern = this_process.name + fn_subdir + "*_wAcc_optics-group*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_xyAcc_optics-group*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_aberr-Axx_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_aberr-Axy_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_aberr-Ayy_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_aberr-bx_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_aberr-by_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_mag_optics-group_*.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_fit.star";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_fit.eps";
            fn_pattern.globFiles(fns_del, false);
        }
    } else if (this_process.type == Process::MOTIONREFINE) {
        for (const FileName &fn_subdir : fns_subdir) {
            // remove the temporary output files
            fn_pattern = this_process.name + fn_subdir + "*_FCC_cc.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_FCC_w0.mrc";
            fn_pattern.globFiles(fns_del, false);
            fn_pattern = this_process.name + fn_subdir + "*_FCC_w1.mrc";
            fn_pattern.globFiles(fns_del, false);

            if (do_harsh) {
                fn_pattern = this_process.name + fn_subdir + "*_shiny.mrcs";
                fn_pattern.globFiles(fns_del, false);
                fn_pattern = this_process.name + fn_subdir + "*_shiny.star";
                fn_pattern.globFiles(fns_del, false);
            }
        }
    } else if (this_process.type == Process::SUBTRACT) {
        if (do_harsh) {
            fn_pattern = this_process.name + "subtracted.*";
            fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
        }
    } else if (this_process.type == Process::POST) {
        fn_pattern = this_process.name + "*masked.mrc";
        fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
    }

    // Now actually move all the files
    FileName fn_old_dir = "";
    // Loop over all files to be deleted
    for (const FileName &fn_del : fns_del) {
        FileName fn_dest = "Trash/" + fn_del;
        FileName fn_dir = fn_dest.beforeLastOf("/");
        if (fn_dir != fn_old_dir && !exists(fn_dir))
            int res = mktree(fn_dir);
        // by removing entire directories, it could be the file is gone already
        if (exists(fn_del)) {
            std::string command = "mv -f " + fn_del + " "+ fn_dir;
            int res = system(command.c_str());
        }
    }
}

// Clean upintermediate files from all jobs in the pipeline
void PipeLine::cleanupAllJobs(bool do_harsh) throw (std::string) {
    for (int myjob = 0; myjob < processList.size(); myjob++) {
        if (processList[myjob].status == Process::FINISHED_SUCCESS) {
            if (do_harsh && exists(processList[myjob].name + "NO_HARSH_CLEAN"))
                continue;
            cleanupJob(myjob, do_harsh);
        }
    }
}

void PipeLine::replaceFilesForImportExportOfScheduledJobs(
    FileName fn_in_dir, FileName fn_out_dir,
    std::vector<std::string> &find_pattern, std::vector<std::string> &replace_pattern
) {
    int res;
    std::string command;
    std::vector<std::string> fns {"run.job", "note.txt", "job_pipeline.star"};

    // Copy the run.job, the note.txt and the job_pipeline.star
    // Replace all instances of all find_pattern's with the replace_pattern's
    for (const std::string &fn : fns) {
    for (int ipatt = 0; ipatt < find_pattern.size(); ipatt++) {
        FileName outfile = fn_out_dir + fn;
        FileName tmpfile = fn_out_dir + "tmp";
        FileName infile = ipatt == 0 ? fn_in_dir + fn : tmpfile;
        // Create directory first time round
        if (ipatt == 0) {
            FileName dirs = outfile.beforeLastOf("/");
            system(("mkdir -p " + dirs).c_str());
        }
        command = "sed 's|" + find_pattern[ipatt] + "|" + replace_pattern[ipatt] + "|g' < " + infile + " > " + outfile;
        // std::cerr << " Executing: " << command<<std::endl;
        res = system(command.c_str());
        if (ipatt + 1 < find_pattern.size()) {
            std::rename(outfile.c_str(), tmpfile.c_str());
            // std::cerr << " Excuting: mv " << outfile<<" "<<tmpfile<<std::endl;
        }
    }
    }
}

void PipeLine::exportAllScheduledJobs(std::string mydir) throw (std::string) {
    mydir.append("/");  // Make sure the directory name ends with a slash
    system(("mkdir -p ExportJobs/" + mydir).c_str());

    MetaDataTable MDexported;

    // Loop through all the Scheduled jobs and export them one-by-one
    int iexp = 0;
    std::vector<std::string> find_pattern, replace_pattern;
    for (const auto &process : processList) {
        if (process.status == Process::SCHEDULED) {
            iexp++;
            if (process.alias != "None")
            throw "ERROR: aliases are not allowed on Scheduled jobs that are to be exported! Make sure all scheduled jobs are made with unaliases names.";

            // A general name for the exported job:
            const FileName expname = FileName(process.name).beforeFirstOf("/") + "/exp" + integerToString(iexp, 3) + "/";
            find_pattern.push_back(process.name);
            replace_pattern.push_back(expname);

            MDexported.setValue(EMDL::PIPELINE_PROCESS_NAME, expname, MDexported.addObject());

            // Copy the run.job, the note.txt and the job_pipeline.star and replace patterns
            replaceFilesForImportExportOfScheduledJobs(process.name, "ExportJobs/" + mydir + expname, find_pattern, replace_pattern);
        }
    }

    MDexported.write("ExportJobs/" + mydir + "exported.star");
}

void PipeLine::importJobs(FileName fn_export) {

    FileName fn_export_dir = fn_export.beforeLastOf("/") + "/";

    //FileName fn_dir_export = fn_export.beforeLastOf("/")+"/";
    MetaDataTable MDexported;
    MDexported.read(fn_export);

    // Read in existing pipeline, in case some other window had changed it
    std::string lock_message = "importJobs";
    read(DO_LOCK, lock_message);

    std::vector<std::string> find_pattern, replace_pattern;
    for (long int i : MDexported) {
        FileName expname = MDexported.getValue<std::string>(EMDL::PIPELINE_PROCESS_NAME, i);
        find_pattern.push_back(expname);
        // Make a new name for this job
        FileName newname = expname.beforeFirstOf("/") + "/job" + integerToString(job_counter, 3) + "/";
        //std::cerr << " expname= " << expname << " newname= " << newname << std::endl;
        replace_pattern.push_back(newname);
        replaceFilesForImportExportOfScheduledJobs(fn_export_dir + expname, newname, find_pattern, replace_pattern);
        // Import the job into the pipeline
        importPipeline(newname + "job");
        job_counter++;
    }

    // Write the new pipeline to disk
    write(DO_LOCK);
}

// Import a job into the pipeline
bool PipeLine::importPipeline(std::string _name) {
    if (_name == name) {
        std::cerr << " importPipeline WARNING: ignoring request to import myself! "<<std::endl;
        return false;
    }

    PipeLine mini_pipeline;
    mini_pipeline.setName(_name);
    mini_pipeline.read();

    // TODO: vectors that map imported process and nodes numbers to the new ones!!
    std::vector<long int> imported_process_nr;
    std::vector<long int> imported_node_nr;
    long int ori_nr_proc = processList.size();
    long int ori_nr_node = nodeList.size();

    bool imported = false;
    for (auto &process : mini_pipeline.processList) {
        // Check that the new processes all have unique names
        if (findProcessByName(process.name) >= 0)
            REPORT_ERROR("importPipeline ERROR: cannot import pipeline with non-unique job name: " + process.name);

        if (findProcessByAlias(process.alias) >= 0) {
            std::cerr << "importPipeline WARNING: resetting non-unique imported alias: " << process.alias << std::endl;
            process.alias = "None";
        }

        imported_process_nr.push_back(processList.size());
        processList.push_back(process);
        imported = true;

    }

    if (imported) {
        for (const auto &node : mini_pipeline.nodeList) {
            // Only push_back nodes that weren't present yet
            int mynode = findNodeByName(node.name);
            if (mynode < 0) {
                imported_node_nr.push_back(nodeList.size());
                nodeList.push_back(node);
            } else {
                //std::cerr << "mynode= "<<mynode << " name=" << nodeList[mynode].name << std::endl;
                imported_node_nr.push_back(mynode);
            }
        }

        // Now fix the numbers of the lists in the imported processes and nodes
        for (auto &process : processList) {
            for (long int &node : process.inputNodeList) {
                int ori_node = node;
                //std::cerr << " ori_node=" << ori_node << " name= "<<nodeList[ori_node].name << std::endl;
                //std::cerr << " imported_node_nr[ori_node]= " << imported_node_nr[ori_node] << " name= "<<nodeList[imported_node_nr[ori_node]].name<<std::endl;
                node = imported_node_nr[ori_node];
            }
            for (long int &node : process.outputNodeList) {
                int ori_node = node;
                node = imported_node_nr[ori_node];
            }
        }
        for (auto &node : nodeList) {
            for (long int &process : node.inputForProcessList) {
                int ori_proc = process;
                process = imported_process_nr[ori_proc];
            }
            int ori_proc2 = node.outputFromProcess;
            node.outputFromProcess = imported_process_nr[ori_proc2];
        }
    }

    return imported;
}

// Read pipeline from STAR file
void PipeLine::read(bool do_lock, std::string lock_message) {

    #ifdef DEBUG_LOCK
    std::cerr << "entering read lock_message=" << lock_message << std::endl;
    #endif
    FileName name_wo_dir = name;
    FileName dir_lock=".relion_lock", fn_lock=".relion_lock/lock_" + name_wo_dir.afterLastOf("/") + "_pipeline.star";;
    if (do_lock && !do_read_only) {
        int iwait = 0;
        int status = mkdir(dir_lock.c_str(), S_IRWXU);

        #ifdef DEBUG_LOCK
        std::cerr <<  " A status= " << status << std::endl;
        #endif
        while (status != 0) {
            if (errno == EACCES) /** NOTE: interestingly, not EACCESS! */
                REPORT_ERROR("ERROR: PipeLine::read cannot create a lock directory " + dir_lock + ". You don't have write permission to this project. If you want to look at other's project directory (but run nothing there), please start RELION with --readonly.");

            // If the lock exists: wait 3 seconds and try again
            // First time round, print a warning message
            if (iwait == 0) {
                std::cout << " WARNING: trying to read pipeline.star, but directory " << dir_lock << " exists (which protects against simultaneous writing by multiple instances of the GUI)" << std::endl;
            }
            sleep(3);
            status =  mkdir(dir_lock.c_str(), S_IRWXU);
            #ifdef DEBUG_LOCK
            std::cerr <<  " B status= " << status << std::endl;
            #endif

            iwait++;
            if (iwait > 40) {

                REPORT_ERROR("ERROR: PipeLine::read has waited for 2 minutes for lock directory to disappear. You may want to manually remove the file: " + fn_lock);
            }

        }
        // Generate the lock file
        std::ofstream  fh;
        fh.open(fn_lock.c_str(), std::ios::out);
        if (!fh) REPORT_ERROR((std::string) "ERROR: Cannot open file: " + fn_lock);
        fh << lock_message << std::endl;
        fh.close();
    }

    // Start from scratch
    clear();

    FileName fn = name + "_pipeline.star";
    std::ifstream in(fn.c_str(), std::ios_base::in);

    if (in.fail()) REPORT_ERROR((std::string) "PipeLine::read: File " + fn + " cannot be read.");

    MetaDataTable MDgen, MDnode, MDproc, MDedge1, MDedge2;

    // This if allows for older version of the pipeline without the jobcounter
    // TODO: remove after alpha-testing
    if (MDgen.readStar(in, "pipeline_general")) {
        int jobcounter = MDgen.getValue<int>(EMDL::PIPELINE_JOB_COUNTER, MDgen.index());
        if (job_counter < 0) REPORT_ERROR("PipeLine::read: rlnPipeLineJobCounter must not be negative!");
    }

    MDnode.readStar(in, "pipeline_nodes");
    for (long int i : MDnode) try {
        std::string name = MDnode.getValue<std::string>(EMDL::PIPELINE_NODE_NAME, i);
        int         type = MDnode.getValue<int>        (EMDL::PIPELINE_NODE_TYPE, i);
        Node newNode(name, type);
        nodeList.push_back(newNode);
    } catch (const char *errmsg) {
        REPORT_ERROR("PipeLine::read: cannot find name or type in pipeline_nodes table");
    }

    MDproc.readStar(in, "pipeline_processes");
    for (long int i : MDproc) {
        try {
            std::string name   = MDproc.getValue<std::string>(EMDL::PIPELINE_PROCESS_NAME, i);
            std::string alias  = MDproc.getValue<std::string>(EMDL::PIPELINE_PROCESS_ALIAS, i);
            int         type   = MDproc.getValue<int>        (EMDL::PIPELINE_PROCESS_TYPE, i);
            int         status = MDproc.getValue<int>        (EMDL::PIPELINE_PROCESS_STATUS, i);
            processList.push_back(Process(name, type, status, alias));

            // Make a symbolic link to the alias if it isn't there...
            if (alias != "None") {
                // Also make a symbolic link for the output directory!
                // Make sure it doesn't end in a slash
                FileName fn_alias = alias;
                if (fn_alias[fn_alias.length() - 1] == '/')
                    fn_alias = fn_alias.beforeLastOf("/");

                // Only make the alias if it doesn't exist yet, otherwise you end up with recursive ones.
                if (!exists(fn_alias)) {
                    // ln -s ../{name} {fn_alias}
                    symlink(("../" + name).c_str(), fn_alias.c_str());
                }
            }

        } catch (const char *errmsg) {
            REPORT_ERROR("PipeLine::read: cannot find name or type in pipeline_processes table");
        }
    }

    // Read in all input (Node->Process) edges
    MDedge1.readStar(in, "pipeline_input_edges");
    for (long int i : MDedge1) {
        std::string fromnodename, procname;
        try {
            procname     = MDedge1.getValue<std::string>(EMDL::PIPELINE_EDGE_PROCESS, i);
            fromnodename = MDedge1.getValue<std::string>(EMDL::PIPELINE_EDGE_FROM, i);
        } catch (const char *errmsg) {
            REPORT_ERROR("PipeLine::read: cannot find procname or fromnodename in pipeline_edges table");
        }

        // Now fill in all To and FromEdgeLists of all Nodes
        long int myProcess = findProcessByName(procname);
        bool found_both = true;
        if (myProcess < 0 || myProcess >= processList.size()) {
            std::cerr << "PipeLine WARNING: cannot find child process with name: " << procname << std::endl;
            found_both = false;
            //REPORT_ERROR("PipeLine::read ERROR: cannot find to-process with name: " + procname);
        }
        long int fromNode = findNodeByName(fromnodename);
        if (fromNode < 0 || fromNode >= nodeList.size()) {
            std::cerr << "PipeLine WARNING: cannot find parent node with name: " << fromnodename << std::endl;
            found_both = false;
            //REPORT_ERROR("PipeLine::read ERROR: cannot find from-node with name: " + fromnodename);
        }
        if (found_both) {
            processList[myProcess].inputNodeList.push_back(fromNode);
            nodeList[fromNode].inputForProcessList.push_back(myProcess);
        }
    }

    // Read in all output (Process->Node) edges
    MDedge2.readStar(in, "pipeline_output_edges");
    for (long int i : MDedge2) {
        std::string tonodename, procname;
        try {
            tonodename = MDedge2.getValue<std::string>(EMDL::PIPELINE_EDGE_TO, i);
            procname   = MDedge2.getValue<std::string>(EMDL::PIPELINE_EDGE_PROCESS, i);
        } catch (const char *errmsg) {
            REPORT_ERROR("PipeLine::read: cannot find procname or tonodename in pipeline_edges table");
        }

        // Now fill in all To and FromEdgeLists of all Nodes
        long int myProcess = findProcessByName(procname);
        bool found_both = true;
        if (myProcess < 0 || myProcess >= processList.size()) {
            std::cerr << "PipeLine WARNING: cannot find parent process with name: " << procname << std::endl;
            found_both = false;
            //REPORT_ERROR("PipeLine::read ERROR: cannot find from-process with name: " + procname);
        }
        long int toNode = findNodeByName(tonodename);
        if (toNode < 0 || toNode >= nodeList.size()) {
            std::cerr << "PipeLine WARNING: cannot find child node with name: " << tonodename << std::endl;
            found_both = false;
            //REPORT_ERROR("PipeLine::read ERROR: cannot find to-node with name: " + tonodename);
        }
        if (found_both) {
            processList[myProcess].outputNodeList.push_back(toNode);
            nodeList[toNode].outputFromProcess = myProcess;
        }
    }
}

void PipeLine::write(
    bool do_lock, FileName fn_del, std::vector<bool> deleteNode, std::vector<bool> deleteProcess
) {
    if (do_read_only) return;

    FileName name_wo_dir = name;
    FileName dir_lock = ".relion_lock", fn_lock = ".relion_lock/lock_" + name_wo_dir.afterLastOf("/") + "_pipeline.star";
    if (do_lock) {

        #ifdef DEBUG_LOCK
        if (exists(fn_lock)) {
            std::cerr << "writing pipeline: " << fn_lock << " exists as expected" << std::endl;
        }
        #endif

        int iwait = 0;
        while (!exists(fn_lock)) {
            // If the lock exists: wait 3 seconds and try again
            // First time round, print a warning message
            if (iwait == 0) {
                std::cerr << " WARNING: was expecting a file called "+fn_lock+ " but it isn't there. Will wait for 1 minute to see whether it appears" << std::endl;
            }
            sleep(3);
            iwait++;
            if (iwait > 40) {
                REPORT_ERROR("ERROR: PipeLine::read has waited for 2 minutes for lock file to appear, but it doesn't. This should not happen. Is something wrong with the disk access?");
            }
        }
    }

    FileName fn = name + "_pipeline.star";
    std::ofstream fh (fn.c_str(), std::ios::out);
    if (fh.fail()) REPORT_ERROR("ERROR: cannot write to pipeline file: " + fn);

    std::ofstream fh_del;
    if (!fn_del.empty()) {
        FileName fnt = fn_del + "deleted_pipeline.star";
        fh_del.open(fnt.c_str(), std::ios::out);
        if (deleteNode.size() != nodeList.size())
            REPORT_ERROR("PipeLine::write BUG: not enough entries in deleteNode vector!");
        if (deleteProcess.size() != processList.size())
            REPORT_ERROR("PipeLine::write BUG: not enough entries in deleteProcess vector!");
    }

    MetaDataTable MDgen, MDnode, MDproc, MDedge1, MDedge2;
    MetaDataTable MDgen_del, MDnode_del, MDproc_del, MDedge1_del, MDedge2_del;

    #ifdef DEBUG
    std::cerr << " writing pipeline as " << fn << std::endl;
    #endif

    MDgen.name = "pipeline_general";
    MDgen.isList = true;
    MDgen.setValue(EMDL::PIPELINE_JOB_COUNTER, job_counter, MDgen.addObject());
    MDgen.write(fh);

    if (!fn_del.empty()) {
        MDgen_del.name = "pipeline_general";
        MDgen_del.isList = true;
        MDgen_del.setValue(EMDL::PIPELINE_JOB_COUNTER, job_counter, MDgen_del.addObject());
        MDgen_del.write(fh_del);
    }

    MDproc.name = MDproc_del.name = "pipeline_processes";
    for (long int i = 0 ; i < processList.size(); i++) {
        const auto &process = processList[i];
        auto &mdt = fn_del.empty() || !deleteProcess[i] ? MDproc : MDproc_del;
        const long int index = mdt.addObject();
        mdt.setValue(EMDL::PIPELINE_PROCESS_NAME,   process.name,   index);
        mdt.setValue(EMDL::PIPELINE_PROCESS_ALIAS,  process.alias,  index);
        mdt.setValue(EMDL::PIPELINE_PROCESS_TYPE,   process.type,   index);
        mdt.setValue(EMDL::PIPELINE_PROCESS_STATUS, process.status, index);
    }
    #ifdef DEBUG
    MDproc.write(std::cerr);
    #endif
    MDproc.write(fh);
    if (!fn_del.empty())
        MDproc_del.write(fh_del);

    MDnode.name = MDnode_del.name = "pipeline_nodes";
    for (long int i = 0; i < nodeList.size(); i++) {
        const auto &node = nodeList[i];
        auto &mdt = fn_del.empty() || !deleteNode[i] ? MDnode : MDnode_del;
        const long int index = mdt.addObject();
        mdt.setValue(EMDL::PIPELINE_NODE_NAME, node.name, index);
        mdt.setValue(EMDL::PIPELINE_NODE_TYPE, node.type, index);
    }
    #ifdef DEBUG
    MDnode.write(std::cerr);
    #endif
    MDnode.write(fh);
    if (!fn_del.empty())
        MDnode_del.write(fh_del);

    // Also write all (Node->Process) edges to a single table
    MDedge1.name = MDedge1_del.name = "pipeline_input_edges";
    for (long int i = 0; i < processList.size(); i++)
    for (long int input_node : processList[i].inputNodeList) {
        auto &mdt = !fn_del.empty() || !deleteProcess[i] && !deleteNode[input_node] ? MDedge1 : MDedge1_del;
        const long int index = mdt.addObject();
        mdt.setValue(EMDL::PIPELINE_EDGE_FROM, nodeList[input_node].name, index);
        mdt.setValue(EMDL::PIPELINE_EDGE_PROCESS, processList[i].name, index);
    }
    #ifdef DEBUG
    MDedge1.write(std::cerr);
    #endif
    MDedge1.write(fh);
    if (!fn_del.empty())
        MDedge1_del.write(fh_del);

    // Also write all (Process->Node) edges to a single table
    MDedge2.name = MDedge2_del.name = "pipeline_output_edges";
    for (long int i = 0; i < processList.size(); i++)
    for (long int output_node : processList[i].outputNodeList) {
        auto &mdt = fn_del.empty() || !deleteProcess[i] && !deleteNode[output_node] ?  MDedge2 : MDedge2_del;
        const long int index = mdt.addObject();
        mdt.setValue(EMDL::PIPELINE_EDGE_PROCESS,  processList[i].name, index);
        mdt.setValue(EMDL::PIPELINE_EDGE_TO, nodeList[output_node].name, index);
    }
    MDedge2.write(fh);
    if (!fn_del.empty()) MDedge2_del.write(fh_del);

    #ifdef DEBUG
    MDedge2.write(std::cerr);
    #endif

    if (!fn_del.empty()) fh_del.close();

    if (do_lock) {

        #ifdef DEBUG_LOCK
        std::cerr << " write pipeline: now deleting " << fn_lock << std::endl;
        #endif

        if (!exists(fn_lock))
            REPORT_ERROR("ERROR: PipeLine::write was expecting a file called "+fn_lock+ " but it is no longer there.");
        if (std::remove(fn_lock.c_str()))
            REPORT_ERROR("ERROR: PipeLine::write reported error in removing file "+fn_lock);
        if (rmdir(dir_lock.c_str()))
            REPORT_ERROR("ERROR: PipeLine::write reported error in removing directory "+dir_lock);
    }

    // Touch a file to indicate to the GUI that the pipeline has just changed
    touch(PIPELINE_HAS_CHANGED);
}

std::string PipeLineFlowChart::getDownwardsArrowLabel(
    PipeLine &pipeline, long int lower_process, long int upper_process
) {
    // What is the type of the node between upper_process and lower_process?
    bool is_found = false;
    long int mynode = -1;
    for (
        long int i = 0; i < pipeline.processList[lower_process].inputNodeList.size(); i++
    ) {
        long int inode= pipeline.processList[lower_process].inputNodeList[i];
        // Find this one in the outputNodeList of the upper_process
        if (pipeline.nodeList[inode].outputFromProcess == upper_process) {
            is_found = true;
            mynode = inode;
            break;
        }
    }

    if (!is_found)
        REPORT_ERROR("PipeLineFlowChart::getDownwardsArrowLabel ERROR: cannot find node connecting " + pipeline.processList[upper_process].name + " -> " + pipeline.processList[lower_process].name);

    MetaDataTable MD;
    long int nr_obj;

    switch (pipeline.nodeList[mynode].type) {

        case Node::MOVIES:
        nr_obj = MD.read(pipeline.nodeList[mynode].name, "", true); // true means: only count nr entries;
        return integerToString(nr_obj) + " movies";

        case Node::MICS:
        nr_obj = MD.read(pipeline.nodeList[mynode].name, "", true); // true means: only count nr entries;
        return integerToString(nr_obj) + " micrographs";

        case Node::PART_DATA:
        nr_obj = MD.read(pipeline.nodeList[mynode].name, "", true); // true means: only count nr entries;
        return integerToString(nr_obj) + " particles";

        case Node::REFS2D:
        return "2Drefs";

        case Node::REF3D:
        return "3D ref";

        case Node::MASK:
        return "mask";

        case Node::MODEL:
        nr_obj = MD.read(pipeline.nodeList[mynode].name, "model_classes", true); // true means: only count nr entries;
        return integerToString(nr_obj) + " classes";

        case Node::OPTIMISER:
        return "continue";

        case Node::HALFMAP:
        return "half-map";

        case Node::FINALMAP:
        return "final map";

        case Node::RESMAP:
        return "local-res map";

        default:
        return "";
    }
}

void PipeLineFlowChart::adaptNamesForTikZ(FileName &name) {
    name.replaceAllSubstrings((std::string) "_", (std::string) "\\_");
    name.replaceAllSubstrings((std::string) ".", (std::string) "-");
    name.replaceAllSubstrings((std::string) ",", (std::string) "-");
    name.replaceAllSubstrings((std::string) "^", (std::string) "\\textasciicircum ");
    name.replaceAllSubstrings((std::string) "~", (std::string) "\\textasciitilde ");
}

long int PipeLineFlowChart::addProcessToUpwardsFlowChart(
    std::ofstream &fh, PipeLine &pipeline,
    long int lower_process, long int new_process, std::vector<long int> &branched_procs
) {
    branched_procs.clear();
    FileName procname = pipeline.processList[new_process].alias_or_name();

    if (do_short_names) {
        procname = procname.beforeFirstOf("/");
    } else {
        FileName longname = procname.afterFirstOf("/").beforeLastOf("/");
        adaptNamesForTikZ(longname);
        procname = procname.beforeFirstOf("/") + "\\\\" + longname;
    }

    FileName new_nodename= pipeline.processList[new_process].name;
    adaptNamesForTikZ(new_nodename);
    FileName lower_nodename;

    // First put the box of the process
    // If this is the lowest process, don't use "above-of" statement, and don't draw an arrow
    if (lower_process < 0) {
        fh << "\\node [block] (" << new_nodename << ") {" << procname << "};" << std::endl;
    } else {
        lower_nodename = pipeline.processList[lower_process].name;
        adaptNamesForTikZ(lower_nodename);

        fh << "\\node [block, above of="<< lower_nodename <<"] (" << new_nodename << ") {" << procname << "};" << std::endl;
        std::string mylabel = getDownwardsArrowLabel(pipeline, lower_process, new_process);
        // Make an arrow from the box to the node it came from
        fh << "\\path [line] ("<< new_nodename <<") -- node[right] {" << mylabel << "} ("<< lower_nodename <<");" << std::endl;
    }

    // See if there are any branchings side-wards, e.g. masks, 2D/3D references, coords, model, optimiser, etc
    long int result = -1;
    if (pipeline.processList[new_process].inputNodeList.empty()) {
        // Reached the top of the tree!
        return -1;
    }
    if (pipeline.processList[new_process].inputNodeList.size() > 1) {

        std::string rightname, leftname;
        for (int inode = 0; inode < pipeline.processList[new_process].inputNodeList.size(); inode++) {
            bool is_left = false;
            bool is_right = false;
            bool is_upper_left = false;
            bool is_upper_right = false;
            std::string right_label="", left_label="";

            long int inputnode = pipeline.processList[new_process].inputNodeList[inode];
            int mynodetype = pipeline.nodeList[inputnode].type;

            if (pipeline.processList[new_process].type == Process::AUTOPICK) {
                is_right = mynodetype == Node::REFS2D;
                right_label = "2D refs";
            } else if (pipeline.processList[new_process].type == Process::EXTRACT) {

                // If the coordinates come from Node::MIC_COORDS, then straight up is the CTF info
                // If the coordinates come from Node::PART_DATA, then that should be straight up
                // therefore first check whether this node has Node::PART_DATA input
                bool has_part_data = false;
                for (int inode2 = 0; inode2 < pipeline.processList[new_process].inputNodeList.size(); inode2++) {
                    long int inputnode2 = pipeline.processList[new_process].inputNodeList[inode2];
                    if (pipeline.nodeList[inputnode2].type == Node::PART_DATA) {
                        has_part_data = true;
                        break;
                    }
                }
                if (has_part_data) {
                    is_right = mynodetype == Node::MICS;
                    right_label = "mics";
                } else {
                    is_right = mynodetype == Node::MIC_COORDS;
                    right_label = "coords";
                }
            } else if (pipeline.processList[new_process].type == Process::CLASS3D) {
                is_right = mynodetype == Node::REF3D;
                right_label = "3D ref";
                is_left = mynodetype == Node::MASK;
                left_label = "mask";
            } else if (pipeline.processList[new_process].type == Process::AUTO3D) {
                is_right = mynodetype == Node::REF3D;
                right_label = "3D ref";
                is_left = mynodetype == Node::MASK;
                left_label = "mask";
            } else if (pipeline.processList[new_process].type == Process::JOINSTAR) {
                // For joinstar: there will be no parent process that returns a postive value!
                // Thereby, joinstar will always end in the 2-4 input processes, each of for which a new flowchart will be made on a new tikZpicture
                if (mynodetype == Node::MOVIES) {
                    right_label = left_label = "mics";
                } else if (mynodetype == Node::PART_DATA) {
                    right_label = left_label = "parts";
                }
                is_right       = inode == 0;
                is_left        = inode == 1;
                is_upper_right = inode == 2;
                is_upper_left  = inode == 3;
            } else if (pipeline.processList[new_process].type == Process::SUBTRACT) {
                is_right = mynodetype == Node::REF3D;
                right_label = "3D ref";
                is_left = mynodetype == Node::MASK;
                left_label = "mask";
            } else if (pipeline.processList[new_process].type == Process::POST) {
                is_left = mynodetype == Node::MASK;
                left_label = "mask";
            } else if (pipeline.processList[new_process].type == Process::RESMAP) {
                is_left = mynodetype == Node::MASK;
                left_label = "mask";
            }

            if (is_right || is_left || is_upper_right || is_upper_left) {
                FileName hyperrefname;
                FileName parent_nodename, newprocname;
                long int parent_process = pipeline.nodeList[inputnode].outputFromProcess;
                if (parent_process < 0) {
                    std::cout << " WARNING: cannot get parent of node: " << pipeline.nodeList[inputnode].name << std::endl;
                    parent_nodename = (is_right || is_upper_right) ? new_nodename + "_rigth" : new_nodename + "_left";
                    newprocname = "unknown";
                } else {
                    // Keep track of all the side-wards branches
                    branched_procs.push_back(parent_process);
                    newprocname = pipeline.processList[parent_process].alias_or_name();
                    if (do_short_names) {
                        newprocname = newprocname.beforeFirstOf("/");
                    } else {
                        FileName longname2 = newprocname.afterFirstOf("/").beforeLastOf("/");
                        adaptNamesForTikZ(longname2);
                        hyperrefname = "sec:" + newprocname.beforeFirstOf("/") + "/" + longname2;
                        if (pipeline.processList[parent_process].type == Process::IMPORT) {
                            newprocname = newprocname.beforeFirstOf("/") + "\\\\" + longname2;
                        } else {
                            newprocname = " \\hyperlink{" + hyperrefname + "}{" + newprocname.beforeFirstOf("/") + "}\\\\" + longname2;
                        }
                    }

                    parent_nodename = pipeline.processList[parent_process].name;
                    adaptNamesForTikZ(parent_nodename);
                    std::string labelpos;
                    if (is_right) { rightname = parent_nodename; }
                    if (is_left)  { leftname  = parent_nodename; }
                }

                if (is_right || is_left) {
                    std::string pos = is_right ? "right" : "left";
                    fh << "\\node [block2, "<< pos << " of=" << new_nodename << "] (" << parent_nodename << ") {" << newprocname << "};" << std::endl;
                } else if (is_upper_right || is_upper_left) {
                    std::string abovename = is_upper_right ? rightname : leftname;
                    fh << "\\node [block2b, above of="<< abovename <<"] (" << parent_nodename << ") {" << newprocname << "};" << std::endl;
                }

                // Make an arrow from the box to the process it came from
                std::string arrowlabel = (is_right || is_upper_right) ? right_label : left_label;
                fh << "\\path [line] ("<< parent_nodename <<") -- node[above] {" << arrowlabel << "} ("<< new_nodename <<");" << std::endl;
            } else {
                result = pipeline.nodeList[inputnode].outputFromProcess;
            }
        }

        return result;
    } else {
        // Only a single input node: return the process that one came from
        long int inputnode = pipeline.processList[new_process].inputNodeList[0];
        return pipeline.nodeList[inputnode].outputFromProcess;
    }
}

void PipeLineFlowChart::makeOneUpwardsFlowChart(
    std::ofstream &fh, PipeLine &pipeline, long int from_process,
    std::vector<long int> &all_branches, bool is_main_flow
) {
    openTikZPicture(fh, is_main_flow);
    long int prev_process = -1;
    long int current_process = from_process;
    bool do_stop = false;
    int counter = 0;
    while (!do_stop) {

        std::vector<long int> branched_procs;
        long int next_process = addProcessToUpwardsFlowChart(
            fh, pipeline, prev_process, current_process, branched_procs
        );

        if (counter > 10) {
            closeTikZPicture(fh, false);
            counter = 0;
            next_process = current_process;
            current_process = prev_process;
            prev_process = -1;
            openTikZPicture(fh, false);
        }

        if (next_process < 0) {
            do_stop = true;
        } else {
            prev_process = current_process;
            current_process = next_process;
        }

        // See if there are any new branches, and if so add them to the all_branches vector
        if (do_branches) {
            for (long int branch :  branched_procs) {
                if (
                    // If the process of this branch is not an import process
                    pipeline.processList[branch].type != Process::IMPORT &&
                    // and the branch is not to be found,
                    std::find(
                        all_branches.begin(), all_branches.end(), branch
                    ) == all_branches.end()
                ) {
                    // add it.
                    all_branches.push_back(branch);
                }
            }
        }

        counter++;

    }
    closeTikZPicture(fh, is_main_flow);
}

void PipeLineFlowChart::makeAllUpwardsFlowCharts(
    FileName &fn_out, PipeLine &pipeline, long int from_process
) {
    std::ofstream fh;
    openFlowChartFile(fn_out, fh);

    // At the beginning of the flowchart file, first make an overview flowchart with short names
    do_short_names = true;
    do_branches = false;
    FileName myorititle = pipeline.processList[from_process].alias_or_name();
    myorititle = myorititle.beforeLastOf("/");
    adaptNamesForTikZ(myorititle);
    fh << "\\section*{Overview flowchart for " << myorititle << "}" << std::endl;
    std::vector<long int> dummy;
    makeOneUpwardsFlowChart(fh, pipeline, from_process, dummy, true);

    // Then, make fully branched flowcharts below
    do_short_names = false;
    do_branches = true;
    std::vector<long int> all_branches {from_process};
    for (int i = 0; i < all_branches.size(); ++i) {
        FileName mytitle = pipeline.processList[all_branches[i]].alias_or_name();
        mytitle = mytitle.beforeLastOf("/");
        adaptNamesForTikZ(mytitle);
        if (i == 0) {
            std::cout << " Making main branched flowchart ... " << std::endl;
            fh << "\\section*{Branched flowchart for " << mytitle << "}" << std::endl;
        } else {
            std::cout << " Making flowchart for branch: " << integerToString(i) << " ... " << std::endl;
            std::string hypertarget = "sec:" + mytitle;
            fh << "\\subsection*{Flowchart for branch " << integerToString(i) << ": " << mytitle << "\\hypertarget{" << hypertarget << "}{}}" << std::endl;
        }

        makeOneUpwardsFlowChart(fh, pipeline, all_branches[i], all_branches, i == 0);
    }

    closeFlowChartFile(fh);
}

void PipeLineFlowChart::openTikZPicture(std::ofstream &fh, bool is_main_flow) {
    if (is_main_flow) {
        fh << "% For large flowcharts: try reducing the fraction on the next line." << std::endl;
        fh << "\\resizebox{!}{0.75\\textheight}{" << std::endl;
    }
    fh << "\\begin{tikzpicture}[scale=1, auto]" << std::endl;
    // Override the long-name styles with the shorter ones
    if (do_short_names) {
        fh << "\\tikzstyle{block} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 1.6cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
        fh << "\\tikzstyle{block2} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 4cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
        fh << "\\tikzstyle{block2b} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 1.6cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
    }
}

void PipeLineFlowChart::closeTikZPicture(std::ofstream &fh, bool is_main_flow) {
    fh << "\\end{tikzpicture}" << std::endl;
    if (is_main_flow) {
        fh << "% For large flowcharts: close resizebox here..." << std::endl;
        fh << "}" << std::endl; // closes resizebox
    }
}

void PipeLineFlowChart::openFlowChartFile(FileName &fn_out, std::ofstream &fh) {

    fh.open(fn_out.c_str(), std::ios::out);
    if (!fh) REPORT_ERROR((std::string) "PipeLineFlowChart ERROR: Cannot write to file: " + fn_out);

    // Set up the LaTex header
    fh << "\\documentclass{article}" << std::endl;
    fh << "\\usepackage{tikz,hyperref}" << std::endl;
    fh << "\\usetikzlibrary{shapes,arrows}" << std::endl;
    fh << "\\begin{document}" << std::endl;
    // These are the styles for the long names!
    fh << "\\tikzstyle{block} = [rectangle, draw, fill=white,text width=3.5cm, node distance = 1.8cm, text centered, rounded corners]" << std::endl;
    fh << "\\tikzstyle{block2} = [rectangle, draw, fill=blue!20,text width=3.5cm, node distance = 5cm, text centered, rounded corners]" << std::endl;
    fh << "\\tikzstyle{block2b} = [rectangle, draw, fill=blue!20,text width=3.5cm, node distance = 1.8cm, text centered, rounded corners]" << std::endl;


    fh << "\\tikzstyle{line} = [draw, very thick, color=black!50, -latex']" << std::endl << std::endl;
}

void PipeLineFlowChart::closeFlowChartFile(std::ofstream &fh) {
    fh << "\\end{document}" << std::endl;
    fh.close();
}
