#include "mrfmaincommand.h"
#include "mrfmaincore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string command_list_message =
    "generate MRF model\n"
    "  build        Build MRF model\n"
    "\n"
    "examine evolutionary information\n"
    "  infer        Estimate sequence distribution with MRF model\n"
    "  stat         Estimate evolutionary constraints for MRF model\n"
    "\n"
    "Use `pmrf --version' to print the PMRF version\n"
    "Use `pmrf -h' or `pmrf --help' to print this help message\n"
    "Use `pmrf <command> -h' to print the help message on a specific command\n";

MRFMainCommandLine::MRFMainCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    subcmd = NONE;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFMainCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFMainProcessor*)processor)->run_mrf_cmd(subcmd);
}

void MRFMainCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " <command> [<args>]" << endl
         << endl
         << "The following commands are available for MRF modeling and varisous applications:" << endl
         << endl
         << command_list_message
         << endl;
}

void MRFMainCommandLine::show_version() {
    cout << PROGNAME << " version " << VERSION << endl;
}

bool MRFMainCommandLine::parse_command_line(int argc, char** argv) {
    if (argc < 2) {
        show_help();
        exit(0);
    }
    if (argv[1][0] != '-') {
        string cmd = argv[1];
        if (cmd == "help") subcmd = HELP;
        else if (cmd == "build") subcmd = BUILD;
        else if (cmd == "infer") subcmd = INFER;
        else if (cmd == "stat") subcmd = STAT;
        else {
            error_message = string("Unknown command: ") + cmd + "\n";
            return false;
        }
        return true;
    }
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"version", 0, 0, 0},
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
                show_version();
                exit(0);
            }
            break;
        case 'h':
            show_help();
            exit(0);
        default:
            return false;
        }
    }
    return true;
}
