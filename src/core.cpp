#include "core.h"

#include <iostream>
#include <cstdlib>

#include "build.h"
#include "analyze.h"

using std::cerr;
using std::endl;

int MRFCmdProcessor::run() {
    if (cmd_line->is_valid()) {
        return cmd_line->process_command(this);
    } else {
        cmd_line->show_error();
        exit(EXIT_FAILURE);
    }
}

void assert_file_handle(std::ifstream& ifs, const std::string& filename) {
    if (!ifs) {
        cerr << "can't open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

void assert_file_handle(std::ofstream& ofs, const std::string& filename) {
    if (!ofs) {
        cerr << "can't open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

/**
   @class MRFMainProcessor
 */

MRFMainProcessor::MRFMainProcessor(int argc, char** argv) {
    cmd_line = new MRFMainCommandLine(argc, argv);
    this->argc = argc;
    this->argv = argv;
}

MRFMainProcessor::~MRFMainProcessor() {
    delete cmd_line;
}

int MRFMainProcessor::run_mrf_cmd(const SubCommand& cmd) {
    argv++;
    argc--;
    switch (cmd) {
    case NONE:
    case HELP:
        cmd_line->show_help();
        exit(0);
    case BUILD:
        return MRFBuildProcessor(argc, argv).run();
    case INFER:
        return MRFInferProcessor(argc, argv).run();
    case STAT:
        return MRFStatProcessor(argc, argv).run();
    }
    return 1;
}

/**
   @class MRFBuildProcessor
 */

MRFBuildProcessor::MRFBuildProcessor(int argc, char** argv) {
    cmd_line = new MRFBuildCommandLine(argc, argv);
}

MRFBuildProcessor::~MRFBuildProcessor() {
    delete cmd_line;
}

int MRFBuildProcessor::build(const string& msa_filename, const string& out_filename) {
    MRFModelBuilder builder(AA);
    builder.opt = ((MRFBuildCommandLine*)cmd_line)->opt;
    int ret = builder.build(msa_filename, out_filename);
    return ret;
}

/**
   @class MRFStatProcessor
 */

MRFStatProcessor::MRFStatProcessor(int argc, char** argv) {
    cmd_line = new MRFStatCommandLine(argc, argv);
}

MRFStatProcessor::~MRFStatProcessor() {
    delete cmd_line;
}

int MRFStatProcessor::stat(const string& mrf_filename) {
    MRFModelAnalyzer analyzer(AA);
    analyzer.opt = ((MRFStatCommandLine*)cmd_line)->opt;
    int ret = analyzer.stat(mrf_filename);
    return ret;
}

/**
   @class MRFInferProcessor
 */

MRFInferProcessor::MRFInferProcessor(int argc, char** argv) {
    cmd_line = new MRFInferCommandLine(argc, argv);
}

MRFInferProcessor::~MRFInferProcessor() {
    delete cmd_line;
}

int MRFInferProcessor::infer(const string& mrf_filename, const string& seq_filename) {
    MRFModelAnalyzer analyzer(AA);
    int ret = analyzer.infer(mrf_filename, seq_filename);
    return ret;
}
