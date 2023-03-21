#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <set>
#include <gaussgen.hpp>
#include "TrPh.h"

const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;

void saveMatrix(const std::string& path, const Eigen::MatrixXd& mx) {
  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file(path);
    if (file.is_open())
    {
        file << mx.format(CSVFormat);
        file.close();
    }
}

Histograms::Histograms() {
    
  // chi2Hist = new TH1F("kf_chi2", "", 128, 0, 70);
 // qHist = new TH1F("qhist", "", 1024, -20, 60);
  /// chi2_vs_q = new TH2F("chi2_vs_q", "", 1024, 0, 70, 1024, -40, 40);
    
    chi2Hist = new TH1F("kf_chi2", "", 128, 0, 400);
  qHist = new TH1F("qhist", "", 1024, -100, 400);
  chi2_vs_q = new TH2F("chi2_vs_q", "", 1024, 0, 400, 1024, -100, 400);
}

Histograms::~Histograms() {
  delete chi2Hist;
  delete qHist;
  delete chi2_vs_q;
}

void Histograms::write(const std::string &path) {
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  chi2Hist->Write();
  qHist->Write();
  chi2_vs_q->Write();
  fl->Close();
   delete fl;
}

TrPh::TrPh(TTree *tree)
    : kfcmd::core::TrPh(tree), entry_(0),
      nevents_(100)  {}

TrPh::~TrPh() {}

bool TrPh::cutTracks_() {
  trackIndices_.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < dZ_) && (std::fabs(trho[i]) < dRho_);
    bool dedx = (tdedx[i] > mindEdX_) && (tdedx[i] < maxdEdX_);
    bool ptot = (tptot[i] > minTPtot_) && (tptot[i] < maxTPtot_);
    if (point && dedx && ptot) trackIndices_.push_back(i);
  }
  if (trackIndices_.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
  return false;
}

Int_t TrPh::Cut(Long64_t) {
  if (nt < 4) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

bool TrPh::detectSimHypo_() {
  sim_hypo_ = -1;
  for (int i = 0; i < nsim; ++i) {
    if (simorig[i] != 0) continue;
    switch (simtype[i]) {
    case 321:
      sim_hypo_ = 0;
      return true;
      break;
    case -321:
      sim_hypo_ = 1;
      return true;
      break;
    default:
      break;
    }
  }
  return false;
}

bool TrPh::fitOnceKPlus_(kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass *hypo) {
  bool result = false;
  double tchi2 = std::numeric_limits<double>::infinity();
  for (int im = 0; im < 2; ++im) {
    if (!hypo->fillTrack("pi-_1", trackIndices_[im], *this))
      continue;
    if (!hypo->fillTrack("pi-_0", trackIndices_[1 - im], *this))
      continue;
    for (int ip = 2; ip < 4; ++ip) {
      if (!hypo->fillTrack("pi+_1", trackIndices_[ip], *this))
        continue;
      if (!hypo->fillTrack("k+", trackIndices_[5 - ip], *this))
        continue;
      hypo->updateInitialParams();
      auto ks_im = hypo->getInitialMomentum("pi-_1") +
                   hypo->getInitialMomentum("pi+_1");
      Eigen::VectorXd tmpv(4);
      tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
      hypo->setInitialParticleParams("ks", tmpv);
      hypo->optimize();
      if (hypo->getErrorCode() != 0) continue;
      if (hypo->getChiSquare() > tchi2) continue;
      tchi2 = hypo->getChiSquare();
      result = true;
      params_.chMinus0_track_index = trackIndices_[1 - im];
      params_.chPlus0_track_index = trackIndices_[5 - ip];
      params_.piMinus1_track_index = trackIndices_[im];
      params_.piPlus1_track_index = trackIndices_[ip];
      fillParamsKPlus_(*hypo);
    }
  }
  return result;
}

void TrPh::fillParamsKPlus_(const kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass& hypo) {
  params_.parChMinus0 = hypo.getParticleFinalParams("pi-_0");
  params_.parChPlus0 = hypo.getParticleFinalParams("k+");
  params_.parPiMinus1 = hypo.getParticleFinalParams("pi-_1");
  params_.parPiPlus1 = hypo.getParticleFinalParams("pi+_1");
  params_.invCovChMinus0 = hypo.getParticleInverseCovarianceMatrix("pi-_0");
  params_.invCovChPlus0 = hypo.getParticleInverseCovarianceMatrix("k+");
  params_.invCovPiMinus1 = hypo.getParticleInverseCovarianceMatrix("pi-_1");
  params_.invCovPiPlus1 = hypo.getParticleInverseCovarianceMatrix("pi+_1");
}

bool TrPh::fitOnceKMinus_(kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass *hypo) {
  bool result = false;
  double tchi2 = std::numeric_limits<double>::infinity();
  for (int im = 0; im < 2; ++im) {
    if (!hypo->fillTrack("pi-_1", trackIndices_[im], *this))
      continue;
    if (!hypo->fillTrack("k-", trackIndices_[1 - im], *this))
      continue;
    for (int ip = 2; ip < 4; ++ip) {
      if (!hypo->fillTrack("pi+_1", trackIndices_[ip], *this))
        continue;
      if (!hypo->fillTrack("pi+_0", trackIndices_[5 - ip], *this))
        continue;
      hypo->updateInitialParams();
      auto ks_im = hypo->getInitialMomentum("pi-_1") +
                   hypo->getInitialMomentum("pi+_1");
      Eigen::VectorXd tmpv(4);
      tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
      hypo->setInitialParticleParams("ks", tmpv);
      hypo->optimize();
      if (hypo->getErrorCode() != 0) continue;
      if (hypo->getChiSquare() > tchi2) continue;
      tchi2 = hypo->getChiSquare();
      result = true;
      params_.chMinus0_track_index = trackIndices_[1 - im];
      params_.chPlus0_track_index = trackIndices_[5 - ip];
      params_.piMinus1_track_index = trackIndices_[im];
      params_.piPlus1_track_index = trackIndices_[ip];
      fillParamsKMinus_(*hypo);
    }
  }
  return result;
}

void TrPh::fillParamsKMinus_(
    const kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass& hypo) {
  params_.parChMinus0 = hypo.getParticleFinalParams("k-");
  params_.parChPlus0 = hypo.getParticleFinalParams("pi+_0");
  params_.parPiMinus1 = hypo.getParticleFinalParams("pi-_1");
  params_.parPiPlus1 = hypo.getParticleFinalParams("pi+_1");
  params_.invCovChMinus0 = hypo.getParticleInverseCovarianceMatrix("k-");
  params_.invCovChPlus0 = hypo.getParticleInverseCovarianceMatrix("pi+_0");
  params_.invCovPiMinus1 = hypo.getParticleInverseCovarianceMatrix("pi-_1");
  params_.invCovPiPlus1 = hypo.getParticleInverseCovarianceMatrix("pi+_1");
}

void TrPh::refitKPlus_(kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass *hypo,
                       Histograms *hists) {
  hypo->setInitialParticleParams("pi-_1",
                                 gaussgen::gaussgen_nfirst(5, params_.parPiMinus1,
                                                           params_.invCovPiMinus1,
                                                           &rnd_));
  hypo->setInitialParticleParams("pi-_0",
                                 gaussgen::gaussgen_nfirst(5, params_.parChMinus0,
                                                           params_.invCovChMinus0,
                                                           &rnd_));
  hypo->setInitialParticleParams("pi+_1",
                                 gaussgen::gaussgen_nfirst(5, params_.parPiPlus1,
                                                           params_.invCovPiPlus1,
                                                           &rnd_));
  hypo->setInitialParticleParams("k+",
                                 gaussgen::gaussgen_nfirst(5, params_.parChPlus0,
                                                           params_.invCovChPlus0,
                                                           &rnd_));
  auto ks_im = hypo->getInitialMomentum("pi-_1") + hypo->getInitialMomentum("pi+_1");
  Eigen::VectorXd tmpv(4);
  tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
  hypo->setInitialParticleParams("ks", tmpv);
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  const double chi2_val = hypo->getChiSquare();
  const double q_val = 0.5 * hypo->getdxTHdx() - hypo->getChiSquare();
  hists->chi2Hist->Fill(chi2_val);
  hists->qHist->Fill(q_val);
  hists->chi2_vs_q->Fill(chi2_val, q_val);
}

void TrPh::refitKMinus_(kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass *hypo,
                       Histograms *hists) {
  hypo->setInitialParticleParams("pi-_1",
                                 gaussgen::gaussgen_nfirst(5, params_.parPiMinus1,
                                                           params_.invCovPiMinus1,
                                                           &rnd_));
  hypo->setInitialParticleParams("k-",
                                 gaussgen::gaussgen_nfirst(5, params_.parChMinus0,
                                                           params_.invCovChMinus0,
                                                           &rnd_));
  hypo->setInitialParticleParams("pi+_1",
                                 gaussgen::gaussgen_nfirst(5, params_.parPiPlus1,
                                                           params_.invCovPiPlus1,
                                                           &rnd_));
  hypo->setInitialParticleParams("pi+_0",
                                 gaussgen::gaussgen_nfirst(5, params_.parChPlus0,
                                                           params_.invCovChPlus0,
                                                           &rnd_));
  auto ks_im = hypo->getInitialMomentum("pi-_1") + hypo->getInitialMomentum("pi+_1");
  Eigen::VectorXd tmpv(4);
  tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
  hypo->setInitialParticleParams("ks", tmpv);
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  const double chi2_val = hypo->getChiSquare();
  const double q_val = hypo->getdxTHdx() - 2 * hypo->getChiSquare();
  hists->chi2Hist->Fill(chi2_val);
  hists->qHist->Fill(q_val);
  hists->chi2_vs_q->Fill(chi2_val, q_val);
}

void TrPh::setEntry(int entry) {
  entry_ = entry;
}

void TrPh::setNEvents(int nevens) {
  nevents_ = nevens;
}

void TrPh::Loop(const std::string &outpath, double magneticField) {
  if (fChain == 0)
    return;
  fChain->GetEntry(0);
  Histograms hists;
  kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass hypo_plus(2.e-3 * emeas,
                                                      magneticField);
  kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass hypo_minus(2.e-3 * emeas,
                                                       magneticField);
  fChain->GetEntry(entry_);
  if (Cut(entry_) < 0) {
    std::cout << "[!] Entry " << entry_ << " didn't pass selection criteria"
              << std::endl;
    return;
  }
  if (!detectSimHypo_()) {
    std::cout << "[!] Simulation hypothesis cannot be detected"
              << std::endl;
    return;
  } else {
    switch (sim_hypo_) {
    case 0:
      std::cout << "sim hypo: KsK+pi-" << std::endl;
      break;
    case 1:
      std::cout << "sim hypo: KsK-pi+" << std::endl;
      break;
    default:
      break;
    }
  }
  hypo_plus.setBeamXY(xbeam, ybeam);
  hypo_minus.setBeamXY(xbeam, ybeam);
  hypo_plus.fixVertexParameter("vtx0", 0, xbeam);
  hypo_plus.fixVertexParameter("vtx0", 1, ybeam);
  hypo_minus.fixVertexParameter("vtx0", 0, xbeam);
  hypo_minus.fixVertexParameter("vtx0", 1, ybeam);
  switch (sim_hypo_) {
  case 0:
    if (!fitOnceKPlus_(&hypo_plus)) {
      std::cout << "[!] Fit isn't converged in KsK+pi- hypothesis"
                << std::endl;
      return;
    }
    for (int event = 0; event < nevents_; ++event) {
      refitKPlus_(&hypo_plus, &hists);
      //
      // if (event == 0) {
	//Eigen::MatrixXd hessian = hypo_plus.getHessian();
//	saveMatrix("hessian.csv", hessian);
  //    }
      //
    }
    hists.write(outpath);
    break;
    case 1:
      if (!fitOnceKMinus_(&hypo_minus)) {
        std::cout << "[!] Fit isn't converged in KsK-pi+ hypothesis"
                  << std::endl;
        return;
      }
      for (int event = 0; event < nevents_; ++event) {
        refitKMinus_(&hypo_minus, &hists);
	//
     // if (event == 0) {
	// Eigen::MatrixXd hessian = hypo_minus.getHessian();
	// saveMatrix("hessian.csv", hessian);
     // }
      //
      }
      hists.write(outpath);
      break;
    default:
      break;
  }
}
