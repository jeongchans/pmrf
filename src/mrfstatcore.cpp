#include "mrfstatcore.h"

#include "analyze.h"

MRFStatProcessor::MRFStatProcessor(int argc, char** argv) {
    cmd_line = new MRFStatCommandLine(argc, argv);
}

MRFStatProcessor::~MRFStatProcessor() {
    delete cmd_line;
}

int MRFStatProcessor::stat(const string& mrf_filename) {
    MRFModelAnalyzer analyzer(AA);
    analyzer.stat_opt = ((MRFStatCommandLine*)cmd_line)->opt.stat_opt;
    int ret = analyzer.stat(mrf_filename);
    return ret;
}
