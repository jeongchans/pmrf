#include "mrfmaincore.h"

#include <cstdlib>

#include "mrfbuildcore.h"
#include "mrfinfercore.h"
#include "mrfstatcore.h"

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
