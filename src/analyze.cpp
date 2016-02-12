#include "analyze.h"

#include <fstream>

#include "mrfio.h"
#include "core.h"

using std::setw;
using std::setprecision;
using std::fixed;
using std::left;
using std::right;

int MRFModelAnalyzer::infer(const string& mrf_filename, const string& seq_filename) {
    MRF model = read_mrf(mrf_filename);
    double wtscore = calc_pll(model, model.get_seq());
    TraceVector traces = read_traces(seq_filename);
    string delim = " ";
    std::cout << setw(4) << right << "#" << delim
              << setw(10) << "Score" << delim
              << setw(10) << "Diff" << delim
              << setw(30) << left << "Description" << std::endl;
    std::cout << "----" << delim
              << "----------" << delim
              << "----------" << delim
              << "--------------------" << std::endl;
    size_t idx = 0;
    for (TraceVector::const_iterator pos = traces.begin(); pos != traces.end(); ++pos) {
        double mtscore = calc_pll(model, pos->get_matched_aseq());
        std::cout << setw(4) << right << ++idx << delim
                  << setw(10) << fixed << setprecision(2) << mtscore << delim
                  << setw(10) << fixed << setprecision(2) << mtscore - wtscore << delim
                  << setw(30) << left << pos->get_desc().substr(0, 30) << std::endl;
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
        int idx1 = model.get_var_idx(aseq[i]);
        pll += model.get_node(i).get_weight()(idx1);
        for (size_t j = i + 1; j < n; ++j) {
            int idx2 = model.get_var_idx(aseq[j]);
            if (model.has_edge(i, j)) pll += model.get_edge(i, j).get_weight()(idx1, idx2);
        }
    }
    return pll;
}
