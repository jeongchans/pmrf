#include "command.h"

bool MRFCommandLine::parse_bool(char* optarg, bool& arg) {
    int n = atoi(optarg);
    if (n == 0) arg = false;
    else if (n == 1) arg = true;
    else {
        error_message = "Not acceptable: " + string(optarg);
        return false;
    }
    return true;
}

bool MRFCommandLine::parse_int(char* optarg, int& arg) {
    arg = atoi(optarg);
    return true;
}

bool MRFCommandLine::parse_float(char* optarg, double& arg) {
    arg = atof(optarg);
    return true;
}

bool MRFCommandLine::parse_str(char* optarg, string& arg) {
    arg = string(optarg);
    return true;
}
