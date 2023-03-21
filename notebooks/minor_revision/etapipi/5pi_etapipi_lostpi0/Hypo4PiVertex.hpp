#ifndef _KFCmd_Hypo4PiVertex_HPP_
#define _KFCmd_Hypo4PiVertex_HPP_
#include "kfcmd/core/Hypothesis.hpp"

namespace kfcmd {
  namespace hypos {
    /**
     * Implementation of K-pi+pi-pi+ hypothesis
     */
    class Hypo4PiVertex : public kfcmd::core::Hypothesis {
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
      Hypo4PiVertex(double, double, long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Hypo4PiVertex();
    };
  } // namespace hypos
}  // namespace kfcmd

#endif
