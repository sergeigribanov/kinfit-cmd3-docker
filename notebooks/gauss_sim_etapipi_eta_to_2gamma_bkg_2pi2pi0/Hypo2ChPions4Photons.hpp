#ifndef _Hypo2ChPions4Photons_HPP_
#define _Hypo2ChPions4Photons_HPP_

#include <kfcmd/core/Hypothesis.hpp>

namespace kfcmd {
namespace hypos {
/**
 * Implementation of (pi+, pi-, gamma, gamma) hypothesis
 */
class Hypo2ChPions4Photons : public kfcmd::core::Hypothesis {
public:
  //! A constructor
  /*!
   * @param energy (center-of-mass energy)
   *
   * @param magneticField (magnetic field)
   *
   * @param nIter (maximum number of iterations)
   *
   * @param tolerance (optimization tolerance)
   */
  Hypo2ChPions4Photons(double, double, long = 20, double = 1.e-4);
  //! A destructor
  virtual ~Hypo2ChPions4Photons();
};
} // namespace hypos
} // namespace kfcmd

#endif
