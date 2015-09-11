#ifndef _BUILD_H_
#define _BUILD_H_

#include <string>
#include <memory>

#include "align/trace.h"
#include "msaanalyze.h"
#include "mrf.h"
#include "parameterize.h"

using std::string;
using std::shared_ptr;

class MRFModelBuilder {
  public:
    MRFModelBuilder(const Alphabet& abc);

    int build(const string& msa_filename, const string& out_filename);

    struct Option {
        MSAFormat msa_fmt;
        MSAAnalyzer::Option msa_analyzer_opt;
        string eidx_filename;
        MRFParameterizer::ObjectiveFunction::Option parameterizer_opt;
        MRFParameterizer::Parameter::Option optim_opt;
    } opt;

  private:
    const Alphabet& abc;

    shared_ptr<MRFParameterizer> mrf_parameterizer;

    TraceVector read_traces(const string& filename);
    MRF build_mrf(const TraceVector& traces);
    void export_mrf(const MRF& model, const string& out_filename);
};

#endif
