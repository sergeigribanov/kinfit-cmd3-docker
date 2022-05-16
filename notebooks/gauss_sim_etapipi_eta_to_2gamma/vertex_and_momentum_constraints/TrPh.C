#include <cmath>
#include <iostream>
#include <gaussgen.hpp>
#include <kfcmd/hypos/Hypo2ChPions2Photons.hpp>
#include "TrPh.h"

const double TrPh::dZ_ = 20;
const double TrPh::dRho_ = 1;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const double TrPh::minPhEn_ = 20;
const double TrPh::maxPhEn_ = 1000;

TrPh::TrPh(TTree *tree)
    : kfcmd::core::TrPh(tree),
      entry_(0),
      nevents_(100000)  {}

TrPh::~TrPh() {}

bool TrPh::cutTracks_() {
  trackIndices_.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < dZ_) && (std::fabs(trho[i]) < dRho_);
    bool dedx = (tdedx[i] > mindEdX_) && (tdedx[i] < maxdEdX_);
    bool ptot = (tptot[i] > minTPtot_) && (tptot[i] < maxTPtot_);
    if (point && dedx && ptot) trackIndices_.push_back(i);
  }
  if (trackIndices_.size() == 2) {
    int totalCharge = 0;
    for (int i = 0; i < 2; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons_() {
  photonIndices_.clear();
  for (int i = 0; i < nph; i++)
    if ((phen[i] > minPhEn_) && (phen[i] < maxPhEn_))
      photonIndices_.push_back(i);
  if (photonIndices_.size() > 1) return true;
  return false;
}

Int_t TrPh::Cut(Long64_t) {
  bool bvar;
 bvar = (nt > 1) && (nph > 1);
  if (!bvar) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
     [&charges](int i, int j) { return charges[i] < charges[j]; });
  if (!cutPhotons_()) return -1;
  return 1;
}

bool TrPh::fitOnce_(kfhypos::Hypo2ChPions2Photons *hypo) {
  bool result = false;
  if (!hypo->fillTrack("pi-", trackIndices_[0], *this)) return false;
  if (!hypo->fillTrack("pi+", trackIndices_[1], *this)) return false;
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

void TrPh::fillParams_(const kfhypos::Hypo2ChPions2Photons &hypo) {
  params_.parPiMi = hypo.getParticleFinalParams("pi-");
  params_.parPiPl = hypo.getParticleFinalParams("pi+");
  params_.parGamma0 = hypo.getParticleFinalParams("g0");
  params_.parGamma1 = hypo.getParticleFinalParams("g1");
  params_.invCovPiMi = hypo.getParticleInverseCovarianceMatrix("pi-");
  params_.invCovPiPl = hypo.getParticleInverseCovarianceMatrix("pi+");
  params_.invCovGamma0 = hypo.getParticleInverseCovarianceMatrix("g0");
  params_.invCovGamma1 = hypo.getParticleInverseCovarianceMatrix("g1");
}

void TrPh::refit_(kfhypos::Hypo2ChPions2Photons *hypo, TH1F *chi2Hist) {
  hypo->setInitialParticleParams("pi-", gaussgen::gaussgen_nfirst(5, params_.parPiMi, params_.invCovPiMi, &rnd_));
  hypo->setInitialParticleParams("pi+", gaussgen::gaussgen_nfirst(5, params_.parPiPl, params_.invCovPiPl, &rnd_));
  hypo->setInitialParticleParams("g0", gaussgen::gaussgen(params_.parGamma0, params_.invCovGamma0, &rnd_));
  hypo->setInitialParticleParams("g1", gaussgen::gaussgen(params_.parGamma1, params_.invCovGamma1, &rnd_));
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  chi2Hist->Fill(hypo->getChiSquare());
}

void TrPh::setEntry(int entry) { entry_ = entry; }

void TrPh::setNEvents(int nevents) { nevents_ = nevents; }

void TrPh::Loop(const std::string &outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  fChain->GetEntry(0);
  TH1F chi2Hist("kf_chi2", "", 512, 0, 128);
  kfcmd::hypos::Hypo2ChPions2Photons hypo(2.e-3 * emeas, magneticField, 100, 1.e-8);
  hypo.disableConstraint("#momentum-constraint-em-vtx0-pe");
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = fChain->GetEntry(entry_);
  if (Cut(entry_) < 0) {
    std::cout << "[!] Entry " << entry_ << " didn't pass selection criteria" << std::endl;
    outfl->Close();
    delete outfl;
    return;
  }
  hypo.setBeamXY(xbeam, ybeam);
  hypo.fixVertexParameter("vtx0", 0, xbeam);
  hypo.fixVertexParameter("vtx0", 0, ybeam);
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
