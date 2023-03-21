#include <cmath>
#include <iostream>
#include <set>
#include <kfcmd/hypos/Hypo2ChPions2Photons.hpp>
#include "TrPh.h"
#include "Hypo2PiVertex.hpp"

const double TrPh::dZ_ = 20;
const double TrPh::dRho_ = 1;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const double TrPh::minPhEn_ = 20;
const double TrPh::maxPhEn_ = 1000;

TrPh::TrPh(TTree* tree) : kfcmd::core::TrPh(tree) {}

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

Int_t TrPh::Cut(Long64_t entry) {
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

void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("evnum", &evnum, "evnum/I");
  tree->Branch("nph", &nph, "nph/I");
  tree->Branch("numiters", &numiters_, "numiters_/I");
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("sigma_x_vtx0", &sigma_x_vtx0_, "sigma_x_vtx0/D");
  tree->Branch("sigma_y_vtx0", &sigma_y_vtx0_, "sigma_y_vtx0/D");
  tree->Branch("sigma_z_vtx0", &sigma_z_vtx0_, "sigma_z_vtx0/D");
  tree->Branch("sim_ee_vtx", &sim_ee_vtx_, "sim_ee_vtx[3]/D");
  tree->Branch("kf_q", &kf_q_, "kf_q/D");
  tree->Branch("in_mgg", &in_mgg_, "in_mgg/D");
  tree->Branch("kf_mgg", &kf_mgg_, "kf_mgg/D");
  tree->Branch("kf_vtx_mgg", &kf_vtx_mgg_, "kf_vtx_mgg/D");
  tree->Branch("in_mpipi", &in_mpipi_, "in_mpipi/D");
  tree->Branch("kf_mpipi", &kf_mpipi_, "kf_mpipi/D");
  tree->Branch("kf_vtx", &kf_vtx_, "kf_vtx[3]/D");
  tree->Branch("in_total_p", in_total_p_, "in_total_p[4]/D");
  tree->Branch("kf_total_p", kf_total_p_, "kf_total_p[4]/D");

  tree->Branch("kf_ct_out_pipl", &kf_ct_out_pipl_, "kf_ct_out_pipl/D");
  tree->Branch("kf_ct_out_pimi", &kf_ct_out_pimi_, "kf_ct_out_pimi/D");
}

void TrPh::fillSimInfo_() {
  for (int i = 0; i < nsim; ++i) {
    switch (simorig[i]) {
      case 0:
        sim_ee_vtx_[0] = simvtx[i];
        sim_ee_vtx_[1] = simvty[i];
        sim_ee_vtx_[2] = simvtz[i];
      default:
        break;
    }
  }
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  const std::set<std::string> sAll = {"pi+", "pi-", "g0", "g1"};
  const std::set<std::string> sGG = {"g0", "g1"};
  const std::set<std::string> sPiPi = {"pi+", "pi-"};
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo2ChPions2Photons hypo(2.e-3 * emeas, magneticField, 100);
  kfcmd::hypos::Hypo2PiVertex hypo_vertex(2.e-3 * emeas, magneticField, 100);
  kfcmd::hypos::Hypo2ChPions2Photons hypo_aux(2.e-3 * emeas, magneticField, 100);

  Eigen::MatrixXd icov = Eigen::MatrixXd::Zero(3, 3);
  icov(0, 0) = 1. / (1.10402760e-01 * 1.10402760e-01);
  icov(1, 1) = 1. / (1.10402760e-01 * 1.10402760e-01);
  hypo.setVertexInverseCovarianceMatrix("vtx0", icov);

  double tchi2;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    hypo.setBeamXY(xbeam, ybeam);
    // hypo.fixVertexParameter("vtx0", 0, xbeam);
    // hypo.fixVertexParameter("vtx0", 1, ybeam);
    hypo_vertex.setBeamXY(xbeam, ybeam);
    // hypo_vertex.fixVertexParameter("vtx1", 0, xbeam);
    // hypo_vertex.fixVertexParameter("vtx1", 1, ybeam);
    kf_err_ = 1;
    numiters_ = -1;
    if (!hypo.fillTrack("pi-", trackIndices_[0], *this)) {
      out_tree->Fill();
      continue;
    }
    if (!hypo.fillTrack("pi+", trackIndices_[1], *this)) {
      out_tree->Fill();
      continue;
    }
    hypo_vertex.fillTrack("pi-_1", trackIndices_[0], *this);
    hypo_vertex.fillTrack("pi+_1", trackIndices_[1], *this);
    hypo_vertex.optimize();
      
    hypo_aux.setBeamXY(xbeam, ybeam);
    hypo_aux.fixVertexParameter("vtx0", 0, hypo_vertex.getFinalVertex("vtx1").X());
    hypo_aux.fixVertexParameter("vtx0", 1, hypo_vertex.getFinalVertex("vtx1").Y());
    hypo_aux.fixVertexParameter("vtx0", 2, hypo_vertex.getFinalVertex("vtx1").Z());  
    hypo_aux.fillTrack("pi-", trackIndices_[0], *this);
    hypo_aux.fillTrack("pi+", trackIndices_[1], *this);
      
    kf_chi2_ = std::numeric_limits<double>::infinity();
    for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
      if (!hypo.fillPhoton("g0", photonIndices_[iph], *this)) continue;
      for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
        if (!hypo.fillPhoton("g1", photonIndices_[jph], *this)) continue;
        hypo.optimize();
        tchi2 = hypo.getChiSquare();
        hypo_aux.fillPhoton("g0", photonIndices_[iph], *this);
        hypo_aux.fillPhoton("g1", photonIndices_[jph], *this);
        hypo_aux.optimize();
        if (kf_err_ != 0 && kf_err_ != 2) {
          numiters_ = hypo.getNumOfRequiredIters();
          kf_err_ = hypo.getErrorCode();
          kf_chi2_ = tchi2;
          kf_q_ = hypo.getdxTHdx() - 2. * tchi2;
          auto inP = hypo.getInitialMomentum(sAll);
          auto kfP = hypo.getFinalMomentum(sAll);
          inP.GetXYZT(in_total_p_);
          kfP.GetXYZT(kf_total_p_);
          in_mgg_ = hypo.getInitialMomentum(sGG).M();
          kf_mgg_ = hypo.getFinalMomentum(sGG).M();
          kf_vtx_mgg_ = hypo_aux.getInitialMomentum(sGG).M();
          in_mpipi_ = hypo.getInitialMomentum(sPiPi).M();
          kf_mpipi_ = hypo.getFinalMomentum(sPiPi).M();
          auto vtx = hypo.getFinalVertex("vtx0");
          vtx.GetXYZ(kf_vtx_);
          kf_ct_out_pipl_ = hypo.getParticleFinalParams("pi+")(5);
          kf_ct_out_pimi_ = hypo.getParticleFinalParams("pi-")(5);
        }
        if (hypo.getErrorCode() != 0) continue;
        if (tchi2 >= kf_chi2_) continue;
        Eigen::MatrixXd cov = 2. * hypo.getInvHessian();
        long bi_vtx0 = hypo.getVertex("vtx0")->getBeginIndex();
        sigma_x_vtx0_ = std::sqrt(std::fabs(cov(bi_vtx0, bi_vtx0)));
        sigma_y_vtx0_ = std::sqrt(std::fabs(cov(bi_vtx0 + 1, bi_vtx0 + 1)));
        sigma_z_vtx0_ = std::sqrt(std::fabs(cov(bi_vtx0 + 2, bi_vtx0 + 2)));
        numiters_ = hypo.getNumOfRequiredIters();
        kf_err_ = hypo.getErrorCode();
        kf_chi2_ = tchi2;
        kf_q_ = hypo.getdxTHdx() - 2. * tchi2;
        auto inP = hypo.getInitialMomentum(sAll);
        auto kfP = hypo.getFinalMomentum(sAll);
        inP.GetXYZT(in_total_p_);
        kfP.GetXYZT(kf_total_p_);
        in_mgg_ = hypo.getInitialMomentum(sGG).M();
        kf_mgg_ = hypo.getFinalMomentum(sGG).M();
        kf_vtx_mgg_ = hypo_aux.getInitialMomentum(sGG).M();
        in_mpipi_ = hypo.getInitialMomentum(sPiPi).M();
        kf_mpipi_ = hypo.getFinalMomentum(sPiPi).M();
        auto vtx = hypo.getFinalVertex("vtx0");
        vtx.GetXYZ(kf_vtx_);
        kf_ct_out_pipl_ = hypo.getParticleFinalParams("pi+")(5);
        kf_ct_out_pimi_ = hypo.getParticleFinalParams("pi-")(5);
      }
    }
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
