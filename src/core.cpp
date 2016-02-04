#include "core.h"

#include <iostream>
#include <cstdlib>

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
