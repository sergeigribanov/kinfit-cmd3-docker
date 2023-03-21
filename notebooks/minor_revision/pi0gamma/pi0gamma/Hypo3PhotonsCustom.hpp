#ifndef _Hypo3PhotonsCustom_HPP_
#define _Hypo3PhotonsCustom_HPP_

#include "kfcmd/core/Hypothesis.hpp"

class Hypo3PhotonsCustom : public kfcmd::core::Hypothesis {
public:
  Hypo3PhotonsCustom(double, double=0., long = 20, double = 1.e-5);
  virtual ~Hypo3PhotonsCustom();
};

#endif
