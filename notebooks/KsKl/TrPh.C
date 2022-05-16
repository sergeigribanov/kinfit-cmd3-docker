#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>

#include "TrPh.h"

const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const std::set<std::string> TrPh::decaypds_ = {"pi-_1", "pi+_1"};

TrPh::TrPh(TTree *tree) :
  kfcmd::core::TrPh(tree) {}

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

Int_t TrPh::Cut(Long64_t) {
  if (nt < 2) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

double TrPh::getAngle_(const kfcmd::hypos::HypoKsKl &hypo) const {
  auto dv = hypo.getFinalVertex("vtx1") - hypo.getFinalVertex("vtx0");
  auto p = hypo.getFinalMomentum(decaypds_).Vect();
  return dv.Angle(p);
}

void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_mks", &in_mks_, "in_mks/D");
  tree->Branch("kf_mks", &kf_mks_, "kf_mks/D");
  tree->Branch("kf_vtx0", kf_vtx0_, "kf_vtx0[3]/D");
  tree->Branch("kf_vtx1", kf_vtx1_, "kf_vtx1[3]/D");
  tree->Branch("kf_vtx_dr", &kf_vtx_dr_, "kf_vtx_dr/D");
  tree->Branch("kf_vtx_drho", &kf_vtx_drho_, "kf_vtx_drho/D");
  tree->Branch("kf_dvtx_2pi_angle",
               &kf_dvtx_2pi_angle_,
               "kf_dvtx_2pi_angle/D");
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree *out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::HypoKsKl hypo(2.e-3 * emeas, magneticField);
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
    if (!hypo.fillTrack("pi-_1", trackIndices_[0], *this)) continue;
    if (!hypo.fillTrack("pi+_1", trackIndices_[1], *this)) continue;
    hypo.updateInitialParams();
    auto ks_im = hypo.getInitialMomentum(decaypds_);
    Eigen::VectorXd tmp_ks_pars(4);
    tmp_ks_pars << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
    hypo.setInitialParticleParams("ks", tmp_ks_pars);
    Eigen::VectorXd tmp_kl_pars(3);
    tmp_kl_pars << -ks_im.Px(), -ks_im.Py(), -ks_im.Pz();
    hypo.setInitialParticleParams("kl", tmp_kl_pars);
    hypo.optimize();
    kf_err_ = hypo.getErrorCode();
    kf_chi2_ = hypo.getChiSquare();
    in_mks_ = hypo.getInitialMomentum(decaypds_).M();
    kf_mks_ = hypo.getFinalMomentum(decaypds_).M();
    auto vtx1 = hypo.getFinalVertex("vtx1");
    auto vtx0 = hypo.getFinalVertex("vtx0");
    vtx0.GetXYZ(kf_vtx0_);
    vtx1.GetXYZ(kf_vtx1_);
    kf_vtx_dr_ = (vtx1 - vtx0).Mag();
    kf_vtx_drho_ = (vtx1 - vtx0).Perp();
    kf_dvtx_2pi_angle_ = getAngle_(hypo);

    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  outfl->Close();
  delete outfl;
}
