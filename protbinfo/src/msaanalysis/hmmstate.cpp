#include "hmmstate.h"

inline void resize_and_fill(Float1dArray& v, const size_t& size, const double& value) {
    v.resize(size);
    v = value;
}

HMMState::HMMState() : eff_num(0) {
    resize_and_fill(transit, NUM_STATE_TYPE, 0);
}

HMMState::HMMState(const size_t& num_emit) : eff_num(0) {
    resize_and_fill(transit, NUM_STATE_TYPE, 0);
    if (num_emit > 0) resize_and_fill(emit, num_emit, 0);
}

HMMState& HMMState::operator=(const HMMState& rhs) {
    emit.resize(rhs.emit.shape());
    set_emit(rhs.emit);
    transit.resize(rhs.transit.shape());
    set_transit(rhs.transit);
    eff_num = rhs.eff_num;
    // TODO: copy State::property
    return *this;
}

void HMMMatchState::accept(HMMStateVisitor* visitor, const size_t& idx) {
    visitor->visit_match(this, idx);
}

void HMMDeleteState::accept(HMMStateVisitor* visitor, const size_t& idx) {
    visitor->visit_delete(this, idx);
}

void HMMInsertState::accept(HMMStateVisitor* visitor, const size_t& idx) {
    visitor->visit_insert(this, idx);
}
