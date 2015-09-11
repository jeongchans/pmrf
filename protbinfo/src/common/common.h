#ifndef _PROTBINFO_COMMON_H_
#define _PROTBINFO_COMMON_H_

#include <utility>
#include <vector>
#include <map>
#include <string>
#include <deque>
#include <exception>

using std::pair;
using std::vector;
using std::deque;
using std::map;
using std::string;
using std::make_pair;

template <typename T1, typename T2>
static inline const pair<T1, T2>& max_pair(const pair<T1, T2>& x, const pair<T1, T2>& y) {
    return (x.first > y.first) ? x : y;
}

template <typename T1, typename T2>
static inline const pair<T1, T2>& max3_pair(const pair<T1, T2>& x, const pair<T1, T2>& y, const pair<T1, T2>& z) {
    return max_pair(max_pair(x, y), z);
}

template <typename T1, typename T2>
static inline const pair<T1, T2>& max4_pair(const pair<T1, T2>& x, const pair<T1, T2>& y, const pair<T1, T2>& z, const pair<T1, T2>& w) {
    return max_pair(max_pair(max_pair(x, y), z), w);
}

template <typename T1, typename T2>
static inline const pair<T1, T2>& max5_pair(const pair<T1, T2>& x, const pair<T1, T2>& y, const pair<T1, T2>& z, const pair<T1, T2>& w, const pair<T1, T2>& v) {
    return max_pair(max_pair(max_pair(max_pair(x, y), z), w), v);
}

template <typename T1, typename T2>
static inline const pair<T1, T2>& max6_pair(const pair<T1, T2>& x, const pair<T1, T2>& y, const pair<T1, T2>& z, const pair<T1, T2>& w, const pair<T1, T2>& v, const pair<T1, T2>& u) {
    return max_pair(max_pair(max_pair(max_pair(max_pair(x, y), z), w), v), u);
}

class ConvergeException : public std::exception {
  public:
    virtual const char* what() const throw() {
        return "Failed to converge";
    }
};

#endif
