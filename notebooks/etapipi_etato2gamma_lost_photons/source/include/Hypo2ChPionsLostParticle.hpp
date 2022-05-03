#ifndef _KFCmd_Hypo2ChPionsLostParticle_HPP_
#define _KFCmd_Hypo2ChPionsLostParticle_HPP_
#include "kfcmd/core/Hypothesis.hpp"

class Hypo2ChPionsLostParticle : public kfcmd::core::Hypothesis {
  public:
  Hypo2ChPionsLostParticle(double, double, long = 20, double = 1.e-5);
  virtual ~Hypo2ChPionsLostParticle();
};

#endif
