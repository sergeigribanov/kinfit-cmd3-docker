#include <iostream>
#include <algorithm>
#include <limits>
#include "TrPh.h"

const double TrPh::min_energy_ = 5.;
const std::set<std::string> TrPh::s_phpair0_ = {"g0", "g1"};
const std::set<std::string> TrPh::s_phpair1_ = {"g1", "g2"};
const std::set<std::string> TrPh::s_phpair2_ = {"g2", "g0"};
const std::set<std::string> TrPh::s_all_photons_ = {"g0", "g1", "g2"};

TrPh::TrPh(TTree *tree)
  : kfcmd::core::TrPh(tree) {}

TrPh::~TrPh() {}

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

Int_t TrPh::Cut(Long64_t) {
  if (nph < 3 || nt > 0) return -1;
  if (!cutPhotons_()) return -1;
  return 1;
}

void TrPh::setupOutputBranches_(TTree *tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("sigma_vtx0_z", &sigma_vtx0_z_, "sigma_vtx0_z/D");
  tree->Branch("in_mgg", in_mgg_, "in_mgg[3]/D");
  tree->Branch("kf_mgg", kf_mgg_, "kf_mgg[3]/D");
  tree->Branch("kf_vtx", kf_vtx_, "kf_vtx[3]/D");
  tree->Branch("in_total_p", in_total_p_, "in_total_p[4]/D");
  tree->Branch("kf_total_p", kf_total_p_, "kf_total_p[4]/D");
  tree->Branch("sim_ee_vtx_z", &sim_ee_vtx_z_, "sim_ee_vtx_z/D");
}

void TrPh::fillSimInfo_() {
  bool flag = false;
  for (int i = 0; i < nsim; ++i) {
    switch (simorig[i]) {
    case 0:
      sim_ee_vtx_z_ = simvtz[i];
      flag = true;
      break;
    default:
      break;
  }
    if (flag) break;
  }
}

void TrPh::fit_(Hypo3PhotonsCustom *hypo) {
  double tchi2;
  kf_chi2_ = std::numeric_limits<double>::infinity();
  fillSimInfo_();
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
    kf_err_ = 1;
    if (hypo->getErrorCode() != 0) {
        if (hypo->getErrorCode() != 1) kf_err_ = hypo->getErrorCode();
        continue;
    }
    kf_err_ = 0;
    tchi2 = hypo->getChiSquare();
    if (tchi2 >= kf_chi2_)
      continue;
    kf_chi2_ = tchi2;
    Eigen::MatrixXd cov = 2. * hypo->getInvHessian();
    long bi_vtx0 = hypo->getVertex("vtx0")->getBeginIndex();
    sigma_vtx0_z_ = std::sqrt(std::fabs(cov(bi_vtx0 + 2, bi_vtx0 + 2)));
    in_mgg_[0] = hypo->getInitialMomentum(s_phpair0_).M();
    in_mgg_[1] = hypo->getInitialMomentum(s_phpair1_).M();
    in_mgg_[2] = hypo->getInitialMomentum(s_phpair2_).M();
    kf_mgg_[0] = hypo->getFinalMomentum(s_phpair0_).M();
    kf_mgg_[1] = hypo->getFinalMomentum(s_phpair1_).M();
    kf_mgg_[2] = hypo->getFinalMomentum(s_phpair2_).M();
    auto vtx = hypo->getFinalVertex("vtx0");
    vtx.GetXYZ(kf_vtx_);
    auto in_p_total = hypo->getInitialMomentum(s_all_photons_) -
      hypo->getInitialMomentum("origin");
    in_p_total.GetXYZT(in_total_p_);
    auto kf_p_total =
      hypo->getFinalMomentum(s_all_photons_) -
      hypo->getFinalMomentum("origin");
    kf_p_total.GetXYZT(kf_total_p_);
  }
}

void TrPh::Loop(const std::string& outpath, double mfield) {
  if (fChain == 0) return;
  fChain->GetEntry(0);
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree *out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  Hypo3PhotonsCustom hypo(2.e-3 * emeas);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    hypo.setBeamXY(xbeam, ybeam);
    hypo.fixVertexParameter("vtx0", 0, xbeam);
    hypo.fixVertexParameter("vtx0", 1, ybeam);
    fit_(&hypo);
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  outfl->Close();
  delete outfl;
}
