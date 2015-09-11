#include "mrfbuildcore.h"

#include "build.h"

MRFBuildProcessor::MRFBuildProcessor(int argc, char** argv) {
    cmd_line = new MRFBuildCommandLine(argc, argv);
}

MRFBuildProcessor::~MRFBuildProcessor() {
    delete cmd_line;
}

int MRFBuildProcessor::build(const string& msa_filename, const string& out_filename) {
//    MRFModelBuilder builder(AMINO_ACID);
    MRFModelBuilder builder(AA_GAP);   // TODO: AA_GAP needs to be replaced with AMINO_ACID
    builder.opt = ((MRFBuildCommandLine*)cmd_line)->opt.build_opt;
    int ret = builder.build(msa_filename, out_filename);
    return ret;
}
