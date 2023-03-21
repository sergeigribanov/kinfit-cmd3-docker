#ifndef _KFCmd_Hypo2PiVertex_HPP_
#define _KFCmd_Hypo2PiVertex_HPP_
#include <kfcmd/core/Hypothesis.hpp>

namespace kfcmd {
  namespace hypos {
    /**
     * Implementation of K-pi+pi-pi+ hypothesis
     */
    class Hypo2PiVertex : public kfcmd::core::Hypothesis {
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
      Hypo2PiVertex(double, double, long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Hypo2PiVertex();
    };
  } // namespace hypos
}  // namespace kfcmd

#endif
