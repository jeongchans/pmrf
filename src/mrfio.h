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
};

class MRFImporter {
  public:
    MRF import_model(istream& is, const Alphabet& abc);
};

ostream& operator<<(ostream& os, const MRF& model);

#endif
