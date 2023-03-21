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

Histograms::Histograms(){
  chi2Hist = new TH1F("kf_chi2", "", 128, 0, 24);
  qHist = new TH1F("qhist", "", 1024, -3, 7);
  sigma_x = new TH1F("sigma_x", "", 1024, 0, 2);
  sigma_y = new TH1F("sigma_y", "", 1024, 0, 2);
  sigma_z = new TH1F("sigma_z", "", 1024, 0, 2);
  vtx_dx = new TH1F("vtx_dx", "", 1024, -20, 20);
  vtx_dy = new TH1F("vtx_dy", "", 1024, -20, 20);
  vtx_dz = new TH1F("vtx_dz", "", 1024, -20, 20);
  vtx_dx_pull = new TH1F("vtx_dx_pull", "", 1024, -20, 20);
  vtx_dy_pull = new TH1F("vtx_dy_pull", "", 1024, -20, 20);
  vtx_dz_pull = new TH1F("vtx_dz_pull", "", 1024, -20, 20);
  chi2_vs_q = new TH2F("chi2_vs_q", "", 1024, 0, 24, 1024, -10, 10);
}

Histograms::~Histograms() {
  delete chi2Hist;
  delete qHist;
  delete sigma_x;
  delete sigma_y;
  delete sigma_z;
  delete vtx_dx;
  delete vtx_dy;
  delete vtx_dz;
  delete vtx_dx_pull;
  delete vtx_dy_pull;
  delete vtx_dz_pull;
  delete chi2_vs_q;
}

void Histograms::write(const std::string &path) {
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  chi2Hist->Write();
  qHist->Write();
  sigma_x->Write();
  sigma_y->Write();
  sigma_z->Write();
  vtx_dx->Write();
  vtx_dy->Write();
  vtx_dz->Write();
  vtx_dx_pull->Write();
  vtx_dy_pull->Write();
  vtx_dz_pull->Write();
  chi2_vs_q->Write();
  fl->Close();
  delete fl;
}

TrPh::TrPh(TTree *tree)
    : kfcmd::core::TrPh(tree), entry_(0), nevents_(100000) {}

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
  auto vtx0 = hypo.getFinalVertex("vtx0");
  params_.parVtx0 = Eigen::VectorXd(3);
  params_.parVtx0 << vtx0.X(), vtx0.Y(), vtx0.Z();
  params_.invCovPiMi = hypo.getParticleInverseCovarianceMatrix("pi-");
  params_.invCovPiPl = hypo.getParticleInverseCovarianceMatrix("pi+");
  params_.invCovGamma0 = hypo.getParticleInverseCovarianceMatrix("g0");
  params_.invCovGamma1 = hypo.getParticleInverseCovarianceMatrix("g1");
  params_.invCovVtx0 = hypo.getVertexInverseCovarianceMatrix("vtx0");
}

void TrPh::refit_(kfhypos::Hypo2ChPions2Photons *hypo, Histograms* hists) {
  hypo->setInitialParticleParams("pi-", gaussgen::gaussgen_nfirst(5, params_.parPiMi, params_.invCovPiMi, &rnd_));
  hypo->setInitialParticleParams("pi+", gaussgen::gaussgen_nfirst(5, params_.parPiPl, params_.invCovPiPl, &rnd_));
  hypo->setInitialParticleParams("g0", gaussgen::gaussgen(params_.parGamma0, params_.invCovGamma0, &rnd_));
  hypo->setInitialParticleParams("g1", gaussgen::gaussgen(params_.parGamma1, params_.invCovGamma1, &rnd_));
  hypo->setInitialVertexParams("vtx0", gaussgen::gaussgen_nfirst(1, params_.parVtx0, params_.invCovVtx0, &rnd_));
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  const double chi2_val = hypo->getChiSquare();
  const double q_val = 0.5 * hypo->getdxTHdx() - hypo->getChiSquare();
  hists->chi2Hist->Fill(chi2_val);
  hists->qHist->Fill(q_val);
  Eigen::MatrixXd cov = 2. * hypo->getInvHessian();
  long bi_vtx0 = hypo->getVertex("vtx0")->getBeginIndex();
  sigma_x_vtx_ = std::sqrt(std::fabs(cov(bi_vtx0, bi_vtx0)));
  sigma_y_vtx_ = std::sqrt(std::fabs(cov(bi_vtx0 + 1, bi_vtx0 + 1)));
  sigma_z_vtx_ = std::sqrt(std::fabs(cov(bi_vtx0 + 2, bi_vtx0 + 2)));
  hists->sigma_x->Fill(sigma_x_vtx_);
  hists->sigma_y->Fill(sigma_y_vtx_);
  hists->sigma_z->Fill(sigma_z_vtx_);
  auto vtx = hypo->getFinalVertex("vtx0");
  vtx.GetXYZ(kf_vtx_);
  hists->vtx_dx->Fill(kf_vtx_[0] - sim_ee_vtx_[0]);
  hists->vtx_dy->Fill(kf_vtx_[1] - sim_ee_vtx_[1]);
  hists->vtx_dz->Fill(kf_vtx_[2] - sim_ee_vtx_[2]);
  hists->vtx_dx_pull->Fill((kf_vtx_[0] - sim_ee_vtx_[0]) / sigma_x_vtx_);
  hists->vtx_dy_pull->Fill((kf_vtx_[1] - sim_ee_vtx_[1]) / sigma_y_vtx_);
  hists->vtx_dz_pull->Fill((kf_vtx_[2] - sim_ee_vtx_[2]) / sigma_z_vtx_);
  hists->chi2_vs_q->Fill(chi2_val, q_val);
}

void TrPh::setEntry(int entry) { entry_ = entry; }

void TrPh::setNEvents(int nevents) { nevents_ = nevents; }

// void TrPh::fillSimInfo_() {
//  for (int i = 0; i < nsim; ++i) {
//    switch (simorig[i]) {
//    case 0:
//      sim_ee_vtx_[0] = simvtx[i];
//      sim_ee_vtx_[1] = simvty[i];
//      sim_ee_vtx_[2] = simvtz[i];
//     default:
//      break;
//    }
//  }
// }

void TrPh::Loop(const std::string &outpath, double magneticField) {
  if (fChain == 0) return;
  fChain->GetEntry(0);
  Histograms hists;
  kfcmd::hypos::Hypo2ChPions2Photons hypo(2.e-3 * emeas, magneticField, 100, 1.e-8);
  Long64_t nb = fChain->GetEntry(entry_);
  if (Cut(entry_) < 0) {
    std::cout << "[!] Entry " << entry_ << " didn't pass selection criteria" << std::endl;
    return;
  }
  // fillSimInfo_();
  hypo.setBeamXY(xbeam, ybeam);
  //  hypo.fixVertexParameter("vtx0", 0, xbeam);
  // hypo.fixVertexParameter("vtx0", 1, ybeam);
  Eigen::MatrixXd icov = Eigen::MatrixXd::Zero(3, 3);
  icov(0, 0) = 1. / (1.10402760e-01 * 1.10402760e-01);
  // icov(1, 1) = 1. / (1.10402760e-01 * 1.10402760e-01);
  hypo.setVertexInverseCovarianceMatrix("vtx0", icov);

  if (!fitOnce_(&hypo)) {
    std::cout << "[!] Fit isn't converged" << std::endl;
    return;
  }
  auto vtx = hypo.getFinalVertex("vtx0");
  vtx.GetXYZ(sim_ee_vtx_);
  for (int event = 0; event < nevents_; ++event) {
    refit_(&hypo, &hists);
  }
  hists.write(outpath);
}
