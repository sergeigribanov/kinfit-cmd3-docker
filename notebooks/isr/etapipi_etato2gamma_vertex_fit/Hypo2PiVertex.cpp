#include "Hypo2PiVertex.hpp"

using namespace kfcmd::hypos;

Hypo2PiVertex::Hypo2PiVertex(double energy, double magnetField,
                             long nIter, double tolerance)
    : kfcmd::core::Hypothesis(energy, magnetField, nIter, tolerance) {
  addVertexXYZ("vtx1");
  auto pipl1 = new kfcmd::core::PiPlusMeson("pi+_1");
  addChargedParticle(pipl1);
  auto pimi1 = new kfcmd::core::PiMinusMeson("pi-_1");
  addChargedParticle(pimi1);
  addOutputVertexConstraintsXYZ("pi+_1", "vtx1");
  addOutputVertexConstraintsXYZ("pi-_1", "vtx1");
}

Hypo2PiVertex::~Hypo2PiVertex() {}