#include "Hypo2ChPions4Photons.hpp"

using namespace kfcmd::hypos;

Hypo2ChPions4Photons::Hypo2ChPions4Photons(double energy,
                                           double magneticField,
                                           long nIter,
                                           double tolerance)
    : kfcmd::core::Hypothesis(energy, magneticField, nIter, tolerance) {
  addVertexXYZ("vtx0");
  auto pip = new kfcmd::core::PiPlusMeson("pi+");
  addChargedParticle(pip);
  auto pim = new kfcmd::core::PiMinusMeson("pi-");
  addChargedParticle(pim);
  addPhoton("g0", "vtx0");
  addPhoton("g1", "vtx0");
  addPhoton("g2", "vtx0");
  addPhoton("g3", "vtx0");
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addEnergyMomentumConstraints("em-vtx0", {getParticle("origin")},
                               {pip, pim, getParticle("g0"), getParticle("g1"),
                                getParticle("g2"), getParticle("g3")});
  addOutputVertexConstraintsXYZ("pi+", "vtx0");
  addOutputVertexConstraintsXYZ("pi-", "vtx0");
}

Hypo2ChPions4Photons::~Hypo2ChPions4Photons() {}
