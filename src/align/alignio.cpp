#include "alignio.h"

#include <climits>

using std::endl;

void FastaAlignmentExporter::export_alignment(const Alignment& aln, std::ostream& os) {
    os << ">" << aln.name << " ";
    export_pos_idx(aln.pos_idx, os);
    os << endl;
    export_seq(aln.seq, os);
}

void FastaAlignmentExporter::export_pos_idx(const std::vector<int>& idxs, std::ostream& os) {
    int beg_idx, end_idx;
    for (std::vector<int>::const_iterator pos = idxs.begin(); pos != idxs.end(); ++pos) {
        if (*pos != GAP) {
            beg_idx = *pos + 1;
            break;
        }
    }
    for (std::vector<int>::const_reverse_iterator pos = idxs.rbegin(); pos != idxs.rend(); ++pos) {
        if (*pos != GAP) {
            end_idx = *pos + 1;
            break;
        }
    }
    os << beg_idx << "-" << end_idx;
}

void FastaAlignmentExporter::export_seq(const std::string seq, std::ostream& os) {
    std::string::size_type beg = 0;
    while (beg < seq.size()) {
        os << seq.substr(beg, width) << endl;
        beg += width;
    }
}
