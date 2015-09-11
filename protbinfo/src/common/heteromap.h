#ifndef _HETEROMAP_H_
#define _HETEROMAP_H_

#include <map>
#include <string>
#include <memory>

#include "numeric.h"

using std::string;
using std::shared_ptr;

class Handle {
  public:
    Handle() {};
    Handle(const FloatType& x) {
        ptr = shared_ptr<FloatType>(new FloatType(x));
    }
    Handle(const Float1dArray& x) {
        ptr = shared_ptr<Float1dArray>(new Float1dArray(x.copy()));
    }
    Handle(const Float2dArray& x) {
        ptr = shared_ptr<Float2dArray>(new Float2dArray(x.copy()));
    }

    shared_ptr<void> ptr;
};

class HeteroMap {
  public:
    typedef string key_type;

    void insert(const key_type& key, const FloatType& val) {
        data.insert(make_pair(key, Handle(val)));
    }
    void insert(const key_type& key, const Float1dArray& val) {
        data.insert(make_pair(key, Handle(val)));
    }
    void insert(const key_type& key, const Float2dArray& val) {
        data.insert(make_pair(key, Handle(val)));
    }
    Handle& operator[](const key_type& x) { return data[x]; }
    const Handle& get(const key_type& x) const { return data.find(x)->second; }

  private:

    std::map<key_type, Handle> data;
};

#endif
