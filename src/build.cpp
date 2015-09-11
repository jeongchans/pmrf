#include "build.h"

#include <fstream>
#include <iostream>

#include "mrfio.h"
#include "core.h"

MRFModelBuilder::MRFModelBuilder(const Alphabet& abc)
: abc(abc) {
}

int MRFModelBuilder::build(const string& filename, const string& out_filename) {
    TraceVector traces = read_traces(filename);
    MRF model = build_mrf(traces);
    export_mrf(model, out_filename);
    return 0;
}

TraceVector MRFModelBuilder::read_traces(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return TraceImporter(abc).import(ifs, opt.msa_fmt).second;
}

MRF MRFModelBuilder::build_mrf(const TraceVector& traces) {
    EdgeIndexVector eidxs;
    EdgeIndexVector* p=NULL;
    if (!opt.eidx_filename.empty()) {
        std::ifstream ifs(opt.eidx_filename.c_str());
        assert_file_handle(ifs, opt.eidx_filename);
        eidxs = EdgeIndexImporter().import(ifs);
        p = &eidxs;
    }
    MSAAnalyzer msa_analyzer(opt.msa_analyzer_opt, abc);
    MRF model(traces[0].get_seq(), abc, p);
    MRFParameterizer mrf_parameterizer(msa_analyzer);
    mrf_parameterizer.opt = opt.parameterizer_opt;
    mrf_parameterizer.optim_opt = opt.optim_opt;
    mrf_parameterizer.parameterize(model, traces);
    return model;
}

void MRFModelBuilder::export_mrf(const MRF& model, const string& out_filename) {
    MRFExporter mrf_exporter;
    if (out_filename.empty()) {
        mrf_exporter.export_model(model, std::cout);
    } else {
        std::ofstream ofs(out_filename.c_str());
        assert_file_handle(ofs, out_filename);
        mrf_exporter.export_model(model, ofs);
    }
}
