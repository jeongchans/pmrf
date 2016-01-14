#include "mrfmaincommand.h"
#include "mrfmaincore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string command_list_message =
    "generate MRF model\n"
    "  build        Build an MRF form an input alignment\n"
    "\n"
    "examine evolutionary information\n"
    "  infer        Calculate a sequence distribution\n";

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
         << "The available commands are:" << endl
         << endl
         << command_list_message
         << endl;
    cout << "Use '" << PROGNAME << " <command> -h' for help on a specific command" << endl;
}

bool MRFMainCommandLine::parse_command_line(int argc, char** argv) {
    optind = 1;
    opterr = 0;     // disable getopt_long error message
    bool help = false;
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
                help = true;
                break;
            }
            break;
        case 'h':
            help = true;
            break;
        }
    }
    if (optind < argc) {
        string s = argv[optind];
        if (s == "help") subcmd = HELP;
        else if (s == "build") subcmd = BUILD;
        else if (s == "infer") subcmd = INFER;
        else {
            error_message = "Unknown subcommand\n";
            return false;
        }
    } else {
        if (help) subcmd = HELP;
        else {
            error_message = "Subcommand is required\n";
            return false;
        }
    }
    return true;
}
