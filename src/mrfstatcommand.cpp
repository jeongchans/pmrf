#include "mrfstatcommand.h"
#include "mrfstatcore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string option_message =
//    "Options:\n"
//    " --mode                Calculation mode of evolutionary constraints\n"
//    "                       pos: positional mode\n"
//    "                       pair: pairwise mode\n"
//    "\n"
    " -h, --help            Help\n";

MRFStatCommandLine::MRFStatCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFStatCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFStatProcessor*)processor)->stat(opt.mrf_filename);
}

void MRFStatCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " stat <mrf_file> [options]" << endl
         << endl
         << option_message
         << endl;
}

bool MRFStatCommandLine::parse_command_line(int argc, char** argv) {
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"mode", required_argument, 0, 0},
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
            case 1:
                //TODO
                return false;
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
        opt.mrf_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}
