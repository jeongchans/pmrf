#ifndef _ALIGNIO_H_
#define _ALIGNIO_H_

#include <istream>
#include <gtest/gtest_prod.h>

#include "alignment.h"

class AlignmentExporter {
  public:
    virtual void export_alignment(const Alignment& aln, std::ostream& os) = 0;
};

class FastaAlignmentExporter : public AlignmentExporter {
  public:
    FastaAlignmentExporter(const size_t& width=80)
    : width(width) {};

    virtual void export_alignment(const Alignment& aln, std::ostream& os);

  private:
    size_t width;

    FRIEND_TEST(FastaAlignmentExporterTest, test_export_pos_idx);

    void export_pos_idx(const std::vector<int>& pos_idx, std::ostream& os);
    void export_seq(const std::string seq, std::ostream& os);
};

#endif
