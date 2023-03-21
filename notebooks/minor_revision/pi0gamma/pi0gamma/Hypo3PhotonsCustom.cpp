#include "Hypo3PhotonsCustom.hpp"

Hypo3PhotonsCustom::Hypo3PhotonsCustom(double energy, double magneticField,
                                       long nIter, double tolerance)
    : kfcmd::core::Hypothesis(energy, magneticField, nIter, tolerance) {
  addVertexXYZ("vtx0");
  addPhoton("g0", "vtx0");
  addPhoton("g1", "vtx0");
  addPhoton("g2", "vtx0");
  addConstantMomentumParticle("origin", energy, Eigen::Vector3d::Zero());
  addEnergyMomentumConstraints("em-constraint", {getParticle("origin")},
                               {getParticle("g0"),
                                getParticle("g1"),
                                getParticle("g2")});
}

Hypo3PhotonsCustom::~Hypo3PhotonsCustom() {}
