#include "mrfinfercore.h"

#include "analyze.h"

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
