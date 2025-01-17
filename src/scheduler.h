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

#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "src/time.h"
#include "src/pipeliner.h"
#include "src/jaz/obs_model.h"
#define SCHEDULE_HAS_CHANGED ".schedule_has_changed";

class SchedulerFloatVariable {

    public:

    RFLOAT value, original_value;

    SchedulerFloatVariable() {};

    SchedulerFloatVariable(RFLOAT _value, RFLOAT _original_value) {
        value = _value;
        original_value = _original_value;
    }
};

class SchedulerBooleanVariable {

    public:
    bool value, original_value;

    SchedulerBooleanVariable() {};

    SchedulerBooleanVariable(bool _value, bool _original_value) {
        value = _value;
        original_value = _original_value;
    }

};

class SchedulerStringVariable {

    public:
    FileName value, original_value;

    SchedulerStringVariable() {};

    SchedulerStringVariable(FileName _value, FileName _original_value) {
        value = _value;
        original_value = _original_value;
    }

};

bool isBooleanVariable(std::string name);
bool isFloatVariable(std::string name);
bool isStringVariable(std::string name);
bool isScheduleOperator(std::string name);

// A class that performs operations on variables
class SchedulerOperator {

    public:

    std::string type, input1, input2, output;

    public:

    SchedulerOperator() {};

    SchedulerOperator(std::string _type, std::string _input1="undefined", std::string _input2="undefined", std::string _output="undefined");

    void initialise(
        const std::string &type, 
        const std::string &input1="undefined",
        const std::string &input2="undefined",
        const std::string &output="undefined"
    );

    // Generate a meaningful current_name for the operator
    std::string getName();

    // Read a specific value from a STAR file
    void readFromStarFile() const;

    bool performOperation() const;

};

class SchedulerJob {

    public:

    std::string current_name, mode;

    bool job_has_started;

    public:

    SchedulerJob() {};

    SchedulerJob(std::string _name, std::string _mode, bool _has_started = false) {
        current_name = _name;
        mode = _mode;
        job_has_started = _has_started;
    }

    // If not a JOB, perform operation and return true.
    // Otherwise, just return false.
    bool performOperation();

};

// Send an email
void schedulerSendEmail(std::string message, std::string subject = "Scheduler");

// A class that defines the edges between a graph that defines execution order, where the nodes are individual JOB instances
// An edge can also be a fork, where the output is controlled through a boolean variable
class SchedulerEdge {

    public:

    std::string inputNode, outputNode, outputNodeTrue;
    std::string myBooleanVariable;
    bool is_fork;

    std::string getOutputNode() const;

    SchedulerEdge(
        std::string _input, std::string _output, bool _is_fork, 
        std::string _mybool, std::string _output_if_true
    ) {
        inputNode = _input;
        outputNode= _output;
        is_fork = _is_fork;
        outputNodeTrue = _output_if_true;
        myBooleanVariable = _mybool;
    }

    SchedulerEdge(std::string _input, std::string _output) {
        inputNode  = _input;
        outputNode = _output;
        is_fork = false;
        outputNodeTrue    = "undefined";
        myBooleanVariable = "undefined";
    }

};

class Schedule {

    public:

    std::string name, current_node;
    bool do_read_only;
    int verb;

    std::map<std::string, SchedulerJob> jobs;
    std::vector<SchedulerEdge> edges;

    PipeLine schedule_pipeline;

    // Operators that return a bool
    static constexpr const char *BOOLEAN_OPERATOR_AND         = "bool=and";
    static constexpr const char *BOOLEAN_OPERATOR_OR          = "bool=or";
    static constexpr const char *BOOLEAN_OPERATOR_NOT         = "bool=not";
    static constexpr const char *BOOLEAN_OPERATOR_GT          = "bool=gt";
    static constexpr const char *BOOLEAN_OPERATOR_LT          = "bool=lt";
    static constexpr const char *BOOLEAN_OPERATOR_GE          = "bool=ge";
    static constexpr const char *BOOLEAN_OPERATOR_LE          = "bool=le";
    static constexpr const char *BOOLEAN_OPERATOR_EQ          = "bool=eq";
    static constexpr const char *BOOLEAN_OPERATOR_FILE_EXISTS = "bool=file_exists";
    static constexpr const char *BOOLEAN_OPERATOR_READ_STAR   = "bool=read_star";

    // Operators that return a float
    static constexpr const char *FLOAT_OPERATOR_SET                      = "float=set";
    static constexpr const char *FLOAT_OPERATOR_PLUS                     = "float=plus";
    static constexpr const char *FLOAT_OPERATOR_MINUS                    = "float=minus";
    static constexpr const char *FLOAT_OPERATOR_MULT                     = "float=mult";
    static constexpr const char *FLOAT_OPERATOR_DIVIDE                   = "float=divide";
    static constexpr const char *FLOAT_OPERATOR_ROUND                    = "float=round";
    static constexpr const char *FLOAT_OPERATOR_COUNT_IMAGES             = "float=count_images";
    static constexpr const char *FLOAT_OPERATOR_COUNT_WORDS              = "float=count_words";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR                = "float=read_star";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_MAX      = "float=star_table_max";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_MIN      = "float=star_table_min";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_AVG      = "float=star_table_avg";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_MAX_IDX  = "float=star_table_max_idx";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_MIN_IDX  = "float=star_table_min_idx";
    static constexpr const char *FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX = "float=star_table_sort_idx";

    // Operators that return a string
    static constexpr const char *STRING_OPERATOR_JOIN         = "string=join";
    static constexpr const char *STRING_OPERATOR_BEFORE_FIRST = "string=before_first";
    static constexpr const char *STRING_OPERATOR_AFTER_FIRST  = "string=after_first";
    static constexpr const char *STRING_OPERATOR_BEFORE_LAST  = "string=before_last";
    static constexpr const char *STRING_OPERATOR_AFTER_LAST   = "string=after_last";
    static constexpr const char *STRING_OPERATOR_READ_STAR    = "string=read_star";
    static constexpr const char *STRING_OPERATOR_GLOB         = "string=glob";
    static constexpr const char *STRING_OPERATOR_NTH_WORD     = "string=nth_word";

    // I/O
    static constexpr const char *OPERATOR_TOUCH_FILE           = "touch_file";
    static constexpr const char *OPERATOR_COPY_FILE            = "copy_file";
    static constexpr const char *OPERATOR_MOVE_FILE            = "move_file";
    static constexpr const char *OPERATOR_DELETE_FILE          = "delete_file";
    static constexpr const char *WAIT_OPERATOR_SINCE_LAST_TIME = "wait";
    static constexpr const char *EMAIL_OPERATOR                = "email";
    static constexpr const char *EXIT_OPERATOR                 = "exit";

    class NodeJobMode {

        public:

        static constexpr const char *NEW       = "new";
        static constexpr const char *CONTINUE  = "continue";
        static constexpr const char *OVERWRITE = "overwrite";

    };

    public:

    Schedule() {
        clear();
    }

    void clear();

    void setName(std::string _name) {
        name = _name;
        schedule_pipeline.setName(_name + "schedule");
    }

    void read(bool do_lock = false, FileName fn = "");

    bool isWriteLocked();
    void write(bool do_lock = false, FileName fn = "");

    void reset();

    void setCurrentNode(std::string _name);
    void setOriginalStartNode(std::string _name);

    bool isNode(std::string name);
    bool isJob(std::string name);
    bool isOperator(std::string name);

    std::string findJobByCurrentName(std::string name);

    // Get/set Variables and Operators(scheduler_floats is only visible in this file!)
    float getFloatVariableValue(std::string name);
    float getFloatOriginalVariableValue(std::string name);
    void setFloatVariableValue(std::string name, RFLOAT val);
    void setFloatOriginalVariableValue(std::string name, RFLOAT val);

    bool getBooleanVariableValue(std::string name);
    bool getBooleanOriginalVariableValue(std::string name);
    void setBooleanVariableValue(std::string name, bool val);
    void setBooleanOriginalVariableValue(std::string name, bool val);

    std::string getStringVariableValue(std::string name);
    std::string getStringOriginalVariableValue(std::string name);
    void setStringVariableValue(std::string name, std::string val);
    void setStringOriginalVariableValue(std::string name, std::string val);

    std::string getVariableValueAsString(std::string name);

    std::string getOperatorName(std::string type, std::string input1, std::string input2, std::string output) {
        SchedulerOperator op(type, input1, input2, output);
        return op.getName();
    }
    void setOperatorParameters(std::string name, std::string type, std::string input1, std::string input2, std::string output);
    void getOperatorParameters(std::string name, std::string &type, std::string &input1, std::string &input2, std::string &output);

    // Get vectors with current Variables / Operators
    std::map<std::string, SchedulerFloatVariable> getCurrentFloatVariables();
    std::map<std::string, SchedulerBooleanVariable> getCurrentBooleanVariables();
    std::map<std::string, SchedulerStringVariable> getCurrentStringVariables();
    std::map<std::string, SchedulerOperator> getCurrentOperators();

    // Get/set operators

    // Add variables
    void setVariable(std::string name, FileName value); // (Add new one if exists, otherwise set value)
    void setOriginalVariable(std::string name, FileName value); // (Add new one if exists, otherwise set original_value)
    void addFloatVariable(std::string name, RFLOAT value);
    void addBooleanVariable(std::string name, bool value);
    void addStringVariable(std::string name, FileName value);

    // Add operators (of any kind), also adds its corresponding node
    SchedulerOperator initialiseOperator(
        std::string type, std::string input_name, std::string input2_name,
        std::string output_name
    );
    void addOperator(SchedulerOperator op);

    // Add a new job, also adds its corresponding node
    void addJob(RelionJob &myjob, std::string jobname, std::string mode);

    void addExitNode();

    // Remove variables/operators/jobs
    void removeVariable(std::string name);
    void removeEdgesWithThisInputOutputOrBoolean(const std::string &name);
    void removeOperator(std::string name);
    void removeOperatorsWithThisInputOrOutput(std::string name);
    void removeJob(std::string name);
    void removeEdge(int idx);

    // Rename this scheduler into a new directory
    void copy(FileName newname);

    // Add edges and forks in between the nodes
    void addEdge(std::string inputnode_name, std::string outputnode_name);
    void addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false );

    // Test integrity of the Schedule. Warn for unused variables, nodes, etc.
    bool isValid();

    std::string getNextNode();
    std::string getPreviousNode();

    bool gotoNextNode();
    bool gotoNextJob();

    // Modify a job to set variables and input nodes from the Scheduler
    void setVariablesInJob(RelionJob &job, FileName original_job_name, bool &needs_a_restart);

    // Run the Schedule
    void run(PipeLine &pipeline);

    // Remove the lock file for read/write protection
    void unlock();

    // Abort a running schedule
    void abort();

    struct rwlock {

        Schedule &schedule;

        rwlock(Schedule& schedule): schedule(schedule) {
            schedule.read(DO_LOCK);
        }

        ~rwlock() {
            schedule.write(DO_LOCK);
        }

    };

};

#endif /* SCHEDULER_H_ */
