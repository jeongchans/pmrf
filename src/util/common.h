#ifndef _PROTBINFO_COMMON_H_
#define _PROTBINFO_COMMON_H_

#include <utility>
#include <vector>
#include <map>
#include <string>
#include <deque>
#include <exception>
#include <iostream>

using std::pair;
using std::vector;
using std::deque;
using std::map;
using std::string;
using std::make_pair;
using std::endl;

class ConvergeException : public std::exception {
  public:
    virtual const char* what() const throw() {
        return "Failed to converge";
    }
};

#endif
