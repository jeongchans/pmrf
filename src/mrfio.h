#ifndef _MRFIO_H_
#define _MRFIO_H_

#include <istream>

#include "mrf.h"

using std::string;
using std::istream;
using std::ostream;

class EdgeIndexImporter {
  public:
    EdgeIndexVector import(std::istream& is);
};

class MRFExporter {
  public:
    void export_model(const MRF& model, ostream& os);

  protected:
    void export_header(const MRF& model, ostream& os);
    void export_body(const MRF& model, ostream& os);

    void export_seq(const string& seq, const size_t& width, ostream& os);
    void export_node_symbol(const string& sym, ostream& os);
    void export_psfm(const Float1dArray& w, ostream& os);
    void export_node_weight(const Float1dArray& w, ostream& os);
    void export_edge_symbol(const string& sym, ostream& os);
    void export_edge_weight(const Float2dArray& w, ostream& os);

    FRIEND_TEST(MRFExporter_Test, test_export_seq);
    FRIEND_TEST(MRFExporter_Test, test_export_node_symbol);
    FRIEND_TEST(MRFExporter_Test, test_export_node_weight);
    FRIEND_TEST(MRFExporter_Test, test_export_edge_symbol);
    FRIEND_TEST(MRFExporter_Test, test_export_edge_weight);

  private:
    void export_elem(const double& x, ostream& os, bool pre_sep=false);
    void export_elem(const char& x, ostream& os, bool pre_sep=false);
    void export_elem(const string& x, ostream& os, bool pre_sep=false);
};

class MRFImporter {
  public:
    MRF import_model(istream& is);
};

#endif
