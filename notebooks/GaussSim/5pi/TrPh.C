#include <TH1F.h>
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <gaussgen.hpp>
#include "TrPh.h"

TrPh::TrPh(TTree *tree)
    : kfcore::TrPh(tree), entry_(0), nevents_(1000000) {}

TrPh::~TrPh() {}

bool TrPh::cutTracks_() {
  trackIndices_.clear();
  for (int i = 0; i < nt; i++) {
    bool point =
      (std::fabs(tz[i]) < 12.0) &&
      (std::fabs(trho[i]) < 1.0);
    bool dedx =
      (tdedx[i] > 0.0) &&
      (tdedx[i] < 15000.0);
    bool ptot =
      (tptot[i] > 5.0) &&
      (tptot[i] < 1000.0);
    if (point && dedx && ptot) trackIndices_.push_back(i);
  }
  if (trackIndices_.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons_() {
  photonIndices_.clear();
  for (int i = 0; i < nph; i++)
    if ((phen[i] > 20.0) && (phen[i] < 1000.0))
      photonIndices_.push_back(i);
  if (photonIndices_.size() >= 2) return true;
  return false;
}

bool TrPh::cut_() {
  if (!cutTracks_()) return false;
  if (!cutPhotons_()) return false;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
	    [&charges](int i, int j) { return charges[i] < charges[j]; });
  return true;
}

bool TrPh::fitOnce_(kfhypos::Hypo4ChPions2Photons* hypo) {
  bool result = false;
  if (!hypo->fillTrack("pi-_0", trackIndices_[0], *this)) return false;
  if (!hypo->fillTrack("pi-_1", trackIndices_[1], *this)) return false;
  if (!hypo->fillTrack("pi+_0", trackIndices_[2], *this)) return false;
  if (!hypo->fillTrack("pi+_1", trackIndices_[3], *this)) return false;
  double tchi2 = std::numeric_limits<double>::infinity();
  for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
    if (!hypo->fillPhoton("g0", photonIndices_[iph], *this)) continue;
    for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
      if (!hypo->fillPhoton("g1", photonIndices_[jph], *this)) continue;
      hypo->optimize();
      if (hypo->getErrorCode() != 0) continue;
      if (hypo->getChiSquare() > tchi2) continue;
      tchi2 = hypo->getChiSquare();
      fillParams_(*hypo);
      result = true;
    }
  }
  return result;
}

void TrPh::fillParams_(const kfhypos::Hypo4ChPions2Photons &hypo) {
  params_.parPiMi0 = hypo.getFinalParameters("pi-_0");
  params_.parPiMi1 = hypo.getFinalParameters("pi-_1");
  params_.parPiPl0 = hypo.getFinalParameters("pi+_0");
  params_.parPiPl1 = hypo.getFinalParameters("pi+_1");
  params_.parGamma0 = hypo.getFinalParameters("g0");
  params_.parGamma1 = hypo.getFinalParameters("g1");
  params_.invCovPiMi0 = hypo.getInverseErrorMatrix("pi-_0");
  params_.invCovPiMi1 = hypo.getInverseErrorMatrix("pi-_1");
  params_.invCovPiPl0 = hypo.getInverseErrorMatrix("pi+_0");
  params_.invCovPiPl1 = hypo.getInverseErrorMatrix("pi+_1");
  params_.invCovGamma0 = hypo.getInverseErrorMatrix("g0");
  params_.invCovGamma1 = hypo.getInverseErrorMatrix("g1");
}

void TrPh::refit_(kfhypos::Hypo4ChPions2Photons* hypo, TH1F* chi2Hist) {
  hypo->setInitialParticleParams("pi-_0", gaussgen::gaussgen(params_.parPiMi0, params_.invCovPiMi0, &rnd_));
  hypo->setInitialParticleParams("pi-_1", gaussgen::gaussgen(params_.parPiMi1, params_.invCovPiMi1, &rnd_));
  hypo->setInitialParticleParams("pi+_0", gaussgen::gaussgen(params_.parPiPl0, params_.invCovPiPl0, &rnd_));
  hypo->setInitialParticleParams("pi+_1", gaussgen::gaussgen(params_.parPiPl1, params_.invCovPiPl1, &rnd_));
  hypo->setInitialParticleParams("g0", gaussgen::gaussgen(params_.parGamma0, params_.invCovGamma0, &rnd_));
  hypo->setInitialParticleParams("g1", gaussgen::gaussgen(params_.parGamma1, params_.invCovGamma1, &rnd_));
  hypo->optimize();
  if(hypo->getErrorCode() != 0) return;
  chi2Hist->Fill(hypo->getChiSquare());
}

void TrPh::setEntry(int entry) { entry_ = entry; }

void TrPh::setNEvents(int nevents) { nevents_ = nevents; }

void TrPh::Loop(const std::string &outpath, double mfield) {
  if
    (fChain == 0) return;
  TFile* outfl = TFile::Open(outpath.c_str(), "recreate");
  TH1F chi2Hist("kf_chi2", "", 128, 0, 128);
  fChain->GetEntry(0);
  kfhypos::Hypo4ChPions2Photons hypo(2 * emeas, mfield);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = fChain->GetEntry(entry_);
  if (!cut_()) {
    std::cout << (boost::format("[!] Entry %1% didn't pass selection criteria") % entry_).str() << std::endl;
    outfl->Close();
    delete outfl;
    return;
  }
  hypo.setBeamXY(xbeam, ybeam);
  hypo.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
  hypo.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
  if (!fitOnce_(&hypo)) {
    std::cout << "[!] Fit isn't converged" << std::endl;
    outfl->Close();
    delete outfl;
    return;
  }
  for (int event = 0; event < nevents_; ++event) {
    refit_(&hypo, &chi2Hist);
  }
  outfl->cd();
  chi2Hist.Write();
  outfl->Close();
  delete outfl;
}
