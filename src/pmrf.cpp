#include "mrfmaincore.h"

int main(int argc, char** argv) {
    MRFMainProcessor processor(argc, argv);
    return processor.run();
}
