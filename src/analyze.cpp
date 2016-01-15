#include "analyze.h"

#include <fstream>

#include "mrfio.h"
#include "core.h"

int MRFModelAnalyzer::infer(const string& mrf_filename, const string& seq_filename) {
    std::ifstream ifs(mrf_filename.c_str());
    assert_file_handle(ifs, mrf_filename);
    MRF model = MRFImporter().import_model(ifs, abc);

    //TODO: read sequence
    //TODO: calcuate pseudolikelihood
    return 0;
}
