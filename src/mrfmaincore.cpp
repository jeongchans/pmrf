#include "mrfmaincore.h"

#include <cstdlib>

#include "mrfbuildcore.h"

MRFMainProcessor::MRFMainProcessor(int argc, char** argv) {
    cmd_line = new MRFMainCommandLine(argc, argv);
    this->argc = argc;
    this->argv = argv;
}

MRFMainProcessor::~MRFMainProcessor() {
    delete cmd_line;
}

int MRFMainProcessor::run_mrf_cmd(const SubCommand& cmd) {
    switch (cmd) {
    case BUILD:
        MRFBuildProcessor processor(argc, argv);
        return processor.run();
    }
}
