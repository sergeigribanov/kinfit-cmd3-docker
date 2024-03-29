/*
 * KFCmd library
 * See LICENSE file at the top of the source tree.
 *
 * This product includes software developed by the
 * CMD-3 collaboration (https://cmd.inp.nsk.su/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 */

/**
 * @file /HypoKsKPlusPiMinus_NoKsMass_ISR.cpp
 *
 * @brief Implementation of HypoKsKPlusPiMinus_NoKsMass_ISR methods
 *
 * @ingroup KFCmd
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "HypoKsKPlusPiMinus_NoKsMass_ISR.hpp"

using namespace kfcmd::hypos;

HypoKsKPlusPiMinus_NoKsMass_ISR::HypoKsKPlusPiMinus_NoKsMass_ISR(double energy, double magnetField,
                                            long nIter, double tolerance)
    : kfcmd::core::Hypothesis(energy, magnetField, nIter, tolerance) {
  addVertexXYZ("vtx0");
  addVertexXYZ("vtx1");
  auto pipl1 = new kfcmd::core::PiPlusMeson("pi+_1");
  addChargedParticle(pipl1);
  auto pimi1 = new kfcmd::core::PiMinusMeson("pi-_1");
  addChargedParticle(pimi1);
  auto pimi0 = new kfcmd::core::PiMinusMeson("pi-_0");
  addChargedParticle(pimi0);
  auto kpl = new kfcmd::core::KPlusMeson("k+");
  addChargedParticle(kpl);
  addIntermediateNeutralParticle("ks", TDatabasePDG::Instance()->GetParticle(310)->Mass(), "vtx0");
  addAltPhoton("isr");
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addMomentumConstraints("momentum-vtx0", {getParticle("origin")},
                         {pimi0, kpl, getParticle("ks"), getParticle("isr")});
  addMomentumConstraints("momentum-vtx1", {getParticle("ks")},
                         {pipl1, pimi1});
  addEnergyConstraint("energy-constraint", {getParticle("origin")},
                      {pimi0, kpl, pipl1, pimi1, getParticle("isr")});
  addOutputVertexConstraintsXYZ("pi-_0", "vtx0");
  addOutputVertexConstraintsXYZ("k+", "vtx0");
  addOutputVertexConstraintsXYZ("pi+_1", "vtx1");
  addOutputVertexConstraintsXYZ("pi-_1", "vtx1");
  addInputVertexConstraintsXYZ("ks", "vtx1");
}

HypoKsKPlusPiMinus_NoKsMass_ISR::~HypoKsKPlusPiMinus_NoKsMass_ISR() {}
