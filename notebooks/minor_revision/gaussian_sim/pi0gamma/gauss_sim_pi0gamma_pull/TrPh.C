#include <cmath>
#include <iostream>
#include <gaussgen.hpp>
#include "Hypo3PhotonsCustom.hpp"
#include "TrPh.h"

const double TrPh::min_energy_ = 5.;
const std::set<std::string> TrPh::s_phpair0_ = {"g0", "g1"};
const std::set<std::string> TrPh::s_phpair1_ = {"g1", "g2"};
const std::set<std::string> TrPh::s_phpair2_ = {"g2", "g0"};
const std::set<std::string> TrPh::s_all_photons_ = {"g0", "g1", "g2"};

bool TrPh::cutPhotons_() {
  photonIndices_.clear();
  for (int i = 0; i < nph; ++i) {
    if (phen[i] > min_energy_) {
      photonIndices_.push_back(i);
    }
  }
  if (photonIndices_.size() != 3) {
    return false;
  }
  return true;
}

Histograms::Histograms(){
  chi2Hist = new TH1F("kf_chi2", "", 128, 0, 24);
  qHist = new TH1F("qhist", "", 1024, -3, 7);
  sigma_z = new TH1F("sigma_z", "", 1024, 0, 2);
  vtx_dz = new TH1F("vtx_dz", "", 1024, -20, 20);
  vtx_dz_pull = new TH1F("vtx_dz_pull", "", 1024, -20, 20);
  chi2_vs_q = new TH2F("chi2_vs_q", "", 1024, 0, 24, 1024, -10, 10);
}

Histograms::~Histograms() {
  delete chi2Hist;
  delete qHist;
  delete sigma_z;
  delete vtx_dz;
  delete vtx_dz_pull;
  delete chi2_vs_q;
}

void Histograms::write(const std::string &path) {
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  chi2Hist->Write();
  qHist->Write();
  sigma_z->Write();
  vtx_dz->Write();
  vtx_dz_pull->Write();
  chi2_vs_q->Write();
  fl->Close();
  delete fl;
}

TrPh::TrPh(TTree *tree)
    : kfcmd::core::TrPh(tree), entry_(0), nevents_(100000) {}

TrPh::~TrPh() {}

Int_t TrPh::Cut(Long64_t) {
  if (nph < 3 || nt > 0) return -1;
  if (!cutPhotons_()) return -1;
  return 1;
}

bool TrPh::fitOnce_(Hypo3PhotonsCustom *hypo) {
  bool result = false;
  double tchi2 = std::numeric_limits<double>::infinity();
  for (std::size_t iph = 0; iph + 2 < photonIndices_.size(); ++iph)
    for (std::size_t jph = iph + 1; jph + 1 < photonIndices_.size(); ++jph)
      for (std::size_t kph = jph + 1; kph < photonIndices_.size(); ++kph) {
	if (!hypo->fillPhoton("g0", photonIndices_[iph], *this))
	  continue;
	if (!hypo->fillPhoton("g1", photonIndices_[jph], *this))
	  continue;
	if (!hypo->fillPhoton("g2", photonIndices_[kph], *this))
	  continue;
	hypo->optimize();
	if (hypo->getErrorCode() != 0) continue;
	if (hypo->getChiSquare() > tchi2) continue;
	tchi2 = hypo->getChiSquare();
	fillParams_(*hypo);
	result = true;
      }
  return result;
}

void TrPh::fillParams_(const Hypo3PhotonsCustom &hypo) {
  params_.parGamma0 = hypo.getParticleFinalParams("g0");
  params_.parGamma1 = hypo.getParticleFinalParams("g1");
  params_.parGamma2 = hypo.getParticleFinalParams("g2");
  params_.invCovGamma0 = hypo.getParticleInverseCovarianceMatrix("g0");
  params_.invCovGamma1 = hypo.getParticleInverseCovarianceMatrix("g1");
  params_.invCovGamma2 = hypo.getParticleInverseCovarianceMatrix("g2");
}

void TrPh::refit_(Hypo3PhotonsCustom *hypo, Histograms* hists) {
  hypo->setInitialParticleParams("g0", gaussgen::gaussgen(params_.parGamma0, params_.invCovGamma0, &rnd_));
  hypo->setInitialParticleParams("g1", gaussgen::gaussgen(params_.parGamma1, params_.invCovGamma1, &rnd_));
  hypo->setInitialParticleParams("g2", gaussgen::gaussgen(params_.parGamma2, params_.invCovGamma2, &rnd_));
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  const double chi2_val = hypo->getChiSquare();
  const double q_val = 0.5 * hypo->getdxTHdx() - hypo->getChiSquare();
  hists->chi2Hist->Fill(chi2_val);
  hists->qHist->Fill(q_val);
  Eigen::MatrixXd cov = 2. * hypo->getInvHessian();
  long bi_vtx0 = hypo->getVertex("vtx0")->getBeginIndex();
  sigma_z_vtx_ = std::sqrt(std::fabs(cov(bi_vtx0 + 2, bi_vtx0 + 2)));
  hists->sigma_z->Fill(sigma_z_vtx_);
  auto vtx = hypo->getFinalVertex("vtx0");
  vtx.GetXYZ(kf_vtx_);
  hists->vtx_dz->Fill(kf_vtx_[2] - sim_ee_vtx_[2]);
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
  Hypo3PhotonsCustom hypo(2.e-3 * emeas);
  Long64_t nb = fChain->GetEntry(entry_);
  if (Cut(entry_) < 0) {
    std::cout << "[!] Entry " << entry_ << " didn't pass selection criteria" << std::endl;
    return;
  }
  // fillSimInfo_();
  hypo.setBeamXY(xbeam, ybeam);
  hypo.fixVertexParameter("vtx0", 0, xbeam);
  hypo.fixVertexParameter("vtx0", 1, ybeam);
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
