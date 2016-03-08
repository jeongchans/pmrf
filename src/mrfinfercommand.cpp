#include "mrfinfercommand.h"
#include "mrfinfercore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string option_message =
    " -h, --help            Help\n";

MRFInferCommandLine::MRFInferCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFInferCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFInferProcessor*)processor)->infer(opt.mrf_filename, opt.seq_filename);
}

void MRFInferCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " infer <file> [options]" << endl
         << endl
         << option_message
         << endl;
}

bool MRFInferCommandLine::parse_command_line(int argc, char** argv) {
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        switch (c) {
        case 0:
            switch (opt_idx) {
            case 0:
                show_help();
                exit(0);
            }
            break;
        case 'h':
            show_help();
            exit(0);
        default:
            error_message = "Unknown option\n";
            return false;
        }
    }
    optind++;
    if (optind < argc) {
        opt.mrf_filename = argv[optind++];
        opt.seq_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}
