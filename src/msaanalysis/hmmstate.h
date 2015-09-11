#ifndef _PROTBINFO_HMMSTATE_H_
#define _PROTBINFO_HMMSTATE_H_

#include "common/numeric.h"
#include "align/trace.h"
#include "common/heteromap.h"

class HMMStateVisitor;

class HMMState {
  public:
    HMMState();
    HMMState(const size_t& num_emit);

    // emission probabilties
    const Float1dArray& get_emit() const { return emit; }
    const FloatType *get_emit_data() const { return emit.data(); }
    void set_emit(const Float1dArray& v) { emit = v; }

    // transition probabilities
    const Float1dArray& get_transit() const { return transit; }
    const FloatType *get_transit_data() const { return transit.data(); }
    FloatType get_transit_to(const StateType& t) const { return transit(t); }
    void set_transit(const Float1dArray& v) { transit = v; }

    // effective number of sequences
    double get_eff_num() const { return eff_num; }
    void set_eff_num(const double& neff) { eff_num = neff; }

    HMMState& operator=(const HMMState& rhs);

    virtual void accept(HMMStateVisitor* visitor, const size_t& idx) = 0;

    // property getter and setter
    void *get_property(const string& key)
        { return property[key].ptr.get(); }
    const void *get_property(const string& key) const
        { return property.get(key).ptr.get(); }
    void set_property(const string& key, const FloatType& value)
        { property.insert(key, value); }
    void set_property(const string& key, const Float1dArray& value)
        { property.insert(key, value); }
    void set_property(const string& key, const Float2dArray& value)
        { property.insert(key, value); }

  private:
    Float1dArray transit;
    Float1dArray emit;
    double eff_num;
    HeteroMap property;
};

class HMMMatchState : public HMMState {
  public:
    HMMMatchState() : HMMState() {};
    HMMMatchState(const size_t& num_emit) : HMMState(num_emit) {};

    virtual void accept(HMMStateVisitor* visitor, const size_t& idx);
};

class HMMDeleteState : public HMMState {
  public:
    HMMDeleteState() : HMMState(0) {};

    virtual void accept(HMMStateVisitor* visitor, const size_t& idx);
};

class HMMInsertState : public HMMState {
  public:
    HMMInsertState() : HMMState() {};
    HMMInsertState(const size_t& num_emit) : HMMState(num_emit) {};

    virtual void accept(HMMStateVisitor* visitor, const size_t& idx);
};

class HMMStateVisitor {
  public:
    virtual void visit_match(HMMMatchState* state, const size_t& idx) = 0;
    virtual void visit_delete(HMMDeleteState* state, const size_t& idx) = 0;
    virtual void visit_insert(HMMInsertState* state, const size_t& idx) = 0;
};

#endif
