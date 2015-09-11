#include "core.h"

#include <iostream>
#include <cstdlib>

using std::cerr;
using std::endl;

AminoAcidGap::AminoAcidGap() : Alphabet("ACDEFGHIKLMNPQRSTVWY-",
                                        "",
                                        "BJZOU",
                                        "X",
                                        "*",
                                        "~") {
    set_degeneracy('B', 'N');
    set_degeneracy('B', 'D');
    set_degeneracy('J', 'I');
    set_degeneracy('J', 'L');
    set_degeneracy('Z', 'Q');
    set_degeneracy('Z', 'E');
    set_degeneracy('O', 'K');
    set_degeneracy('U', 'C');
}

bool AminoAcidGap::is_canonical_aa(const char& x) const {
    return is_canonical(x) && x != '-';
}

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
