#include "analyze.h"

#include <fstream>

#include "mrfio.h"
#include "core.h"

int MRFModelAnalyzer::infer(const string& mrf_filename, const string& seq_filename) {
    MRF model = read_mrf(mrf_filename);
    TraceVector traces = read_traces(seq_filename);
    for (TraceVector::const_iterator pos = traces.begin(); pos != traces.end(); ++pos) {
        double pll = calc_pll(model, pos->get_matched_aseq());
        std::cout << pos->get_matched_aseq() << "\t" << pll << std::endl;
    }
    return 0;
}

MRF MRFModelAnalyzer::read_mrf(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return MRFImporter().import_model(ifs, abc);
}

TraceVector MRFModelAnalyzer::read_traces(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return TraceImporter(abc).import(ifs, AFASTA).second;
}

double MRFModelAnalyzer::calc_pll(const MRF& model, const string& aseq) {
    double pll = 0.;
    size_t n = model.get_length();
    for (size_t i = 0; i < n; ++i) {
        int idx1 = model.get_idx(aseq[i]);
        pll += model.get_node(i).get_weight()(idx1);
        for (size_t j = i + 1; j < n; ++j) {
            int idx2 = model.get_idx(aseq[j]);
            pll += model.get_edge(i, j).get_weight()(idx1, idx2);
        }
    }
    return pll;
}
