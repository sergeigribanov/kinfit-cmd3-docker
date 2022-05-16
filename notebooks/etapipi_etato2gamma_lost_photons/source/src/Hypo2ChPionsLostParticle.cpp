#include "Hypo2ChPionsLostParticle.hpp"
#include "kfcmd/core/PiMinusMeson.hpp"
#include "kfcmd/core/PiPlusMeson.hpp"

Hypo2ChPionsLostParticle::Hypo2ChPionsLostParticle(double energy,
                                                   double magneticField,
                                                   long nIter,
                                                   double tolerance)
    : kfcmd::core::Hypothesis(energy, magneticField, nIter, tolerance) {
  addVertexXYZ("vtx0");
  auto pip = new kfcmd::core::PiPlusMeson("pi+");
  addChargedParticle(pip);
  auto pim = new kfcmd::core::PiMinusMeson("pi-");
  addChargedParticle(pim);
  // Since energy constraint is not used, the mass of particle X dosen't play any role
  addParticlePxPyPz("X", 1.);
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addMomentumConstraints("em-vtx0", {getParticle("origin")},
                         {pip, pim, getParticle("X")});
  addOutputVertexConstraintsXYZ("pi+", "vtx0");
  addOutputVertexConstraintsXYZ("pi-", "vtx0");
}

Hypo2ChPionsLostParticle::~Hypo2ChPionsLostParticle() {}
