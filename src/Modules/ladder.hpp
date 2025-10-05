#pragma once
#include "Coulomb/QkTable.hpp"
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace IO {
class InputBlock;
}
namespace Angular {
class SixJTable;
}

namespace Module {

//! Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

void check_L_symmetry(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &valence,
                      const Coulomb::QkTable &qk, bool include_L4,
                      const Angular::SixJTable &sj,
                      const Coulomb::LkTable *const lk = nullptr);

} // namespace Module
