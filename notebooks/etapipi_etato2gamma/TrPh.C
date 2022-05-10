#include <cmath>
#include <iostream>
#include <set>
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
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_mgg", &in_mgg_, "in_mgg/D");
  tree->Branch("kf_mgg", &kf_mgg_, "kf_mgg/D");
  tree->Branch("in_mpipi", &in_mpipi_, "in_mpipi/D");
  tree->Branch("kf_mpipi", &kf_mpipi_, "kf_mpipi/D");
  tree->Branch("kf_vtx", kf_vtx_, "kf_vtx[3]/D");

  tree->Branch("in_total_px", &in_total_px_, "in_total_px/D");
  tree->Branch("in_total_py", &in_total_py_, "in_total_py/D");
  tree->Branch("in_total_pz", &in_total_pz_, "in_total_pz/D");
  tree->Branch("in_total_pe", &in_total_pe_, "in_total_pe/D");

  tree->Branch("kf_total_px", &kf_total_px_, "kf_total_px/D");
  tree->Branch("kf_total_py", &kf_total_py_, "kf_total_py/D");
  tree->Branch("kf_total_pz", &kf_total_pz_, "kf_total_pz/D");
  tree->Branch("kf_total_pe", &kf_total_pe_, "kf_total_pe/D");

  tree->Branch("in_g0_ct", &in_g0_ct_, "in_g0_ct/D");
  tree->Branch("in_g1_ct", &in_g1_ct_, "in_g1_ct/D");
  tree->Branch("kf_g0_ct", &kf_g0_ct_, "kf_g0_ct/D");
  tree->Branch("kf_g1_ct", &kf_g1_ct_, "kf_g1_ct/D");

  tree->Branch("in_g0_cpt", in_g0_cpt_, "in_g0_cpt[3]/D");
  tree->Branch("in_g1_cpt", in_g1_cpt_, "in_g1_cpt[3]/D");
  tree->Branch("kf_g0_cpt", kf_g0_cpt_, "kf_g0_cpt[3]/D");
  tree->Branch("kf_g1_cpt", kf_g1_cpt_, "kf_g1_cpt[3]/D");

  tree->Branch("in_g0_dir", in_g0_dir_, "in_g0_dir[3]/D");
  tree->Branch("in_g1_dir", in_g1_dir_, "in_g1_dir[3]/D");
  tree->Branch("kf_g0_dir", kf_g0_dir_, "kf_g0_dir[3]/D");
  tree->Branch("kf_g1_dir", kf_g1_dir_, "kf_g1_dir[3]/D");

  tree->Branch("in_g0_tdir", in_g0_tdir_, "in_g0_tdir[3]/D");
  tree->Branch("in_g1_tdir", in_g1_tdir_, "in_g1_tdir[3]/D");
  tree->Branch("kf_g0_tdir", kf_g0_tdir_, "kf_g0_tdir[3]/D");
  tree->Branch("kf_g1_tdir", kf_g1_tdir_, "kf_g1_tdir[3]/D");

  tree->Branch("g0_time", &g0_time_, "g0_time/D");
  tree->Branch("g1_time", &g1_time_, "g1_time/D");
  tree->Branch("g0_dr", &g0_dr_, "g0_dr/D");
  tree->Branch("g1_dr", &g1_dr_, "g1_dr/D");
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
  kfcmd::hypos::Hypo2ChPions2Photons hypo(2.e-3 * emeas, magneticField, 20, 1.e-4);
  // hypo.disableConstraint("#momentum-constraint-em-vtx0-pe");
  // hypo.disableConstraint("#momentum-constraint-em-vtx0-pz");
  // hypo.disableConstraint("#momentum-constraint-em-vtx0-px");
  // hypo.disableConstraint("#momentum-constraint-em-vtx0-py");
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
    hypo.fixVertexParameter("vtx0", 0, xbeam);
    hypo.fixVertexParameter("vtx0", 1, ybeam);
    kf_err_ = 1;
    if (!hypo.fillTrack("pi-", trackIndices_[0], *this)) {
      out_tree->Fill();
      continue;
    }
    if (!hypo.fillTrack("pi+", trackIndices_[1], *this)) {
      out_tree->Fill();
      continue;
    }
    kf_chi2_ = std::numeric_limits<double>::infinity();
    for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
      if (!hypo.fillPhoton("g0", photonIndices_[iph], *this)) continue;
      for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
        if (!hypo.fillPhoton("g1", photonIndices_[jph], *this)) continue;
        hypo.optimize();
        // if (hypo.getErrorCode() != 0) continue;
        tchi2 = hypo.getChiSquare();
        if (tchi2 >= kf_chi2_) continue;
        kf_err_ = hypo.getErrorCode(); // 0
        kf_chi2_ = tchi2;
        auto inP = hypo.getInitialMomentum(sAll);
        auto kfP = hypo.getFinalMomentum(sAll);
        in_total_px_ = inP.Px();
        in_total_py_ = inP.Py();
        in_total_pz_ = inP.Pz();
        in_total_pe_ = inP.E();
        kf_total_px_ = kfP.Px();
        kf_total_py_ = kfP.Py();
        kf_total_pz_ = kfP.Pz();
        kf_total_pe_ = kfP.E();
        in_mgg_ = hypo.getInitialMomentum(sGG).M();
        kf_mgg_ = hypo.getFinalMomentum(sGG).M();
        in_mpipi_ = hypo.getInitialMomentum(sPiPi).M();
        kf_mpipi_ = hypo.getFinalMomentum(sPiPi).M();
        auto vtx = hypo.getFinalVertex("vtx0");
        vtx.GetXYZ(kf_vtx_);
        auto g0 = dynamic_cast<kfcmd::core::Photon *>(hypo.getParticle("g0"));
        auto g1 = dynamic_cast<kfcmd::core::Photon *>(hypo.getParticle("g1"));
        auto in_g0_cpt = g0->getInitialConvPoint();
        in_g0_cpt.GetXYZ(in_g0_cpt_);
        auto in_g1_cpt = g1->getInitialConvPoint();
        in_g1_cpt.GetXYZ(in_g1_cpt_);
        auto kf_g0_cpt = g0->getFinalConvPoint();
        kf_g0_cpt.GetXYZ(kf_g0_cpt_);
        auto kf_g1_cpt = g1->getFinalConvPoint();
        kf_g1_cpt.GetXYZ(kf_g1_cpt_);

        auto in_g0_dir = g0->getInitialDirection();
        in_g0_dir.GetXYZ(in_g0_dir_);
        auto in_g1_dir = g1->getInitialDirection();
        in_g1_dir.GetXYZ(in_g1_dir_);
        auto kf_g0_dir = g0->getFinalDirection();
        kf_g0_dir.GetXYZ(kf_g0_dir_);
        auto kf_g1_dir = g1->getFinalDirection();
        kf_g1_dir.GetXYZ(kf_g1_dir_);

        auto in_g0_tdir = in_g0_cpt * (1. / in_g0_cpt.Mag());
        auto in_g1_tdir = in_g1_cpt * (1. / in_g1_cpt.Mag());

        auto kf_g0_tdir = kf_g0_cpt - vtx;
        kf_g0_tdir *= 1. / kf_g0_tdir.Mag();
        auto kf_g1_tdir = kf_g1_cpt - vtx;
        kf_g1_tdir *= 1. / kf_g1_tdir.Mag();

        in_g0_tdir.GetXYZ(in_g0_tdir_);
        in_g1_tdir.GetXYZ(in_g1_tdir_);
        kf_g0_tdir.GetXYZ(kf_g0_tdir_);
        kf_g1_tdir.GetXYZ(kf_g1_tdir_);

        g0_time_ = g0->getFinalParameters()(6);
        g1_time_ = g1->getFinalParameters()(6);
        g0_dr_ = (kf_g0_cpt - vtx).Mag();
        g1_dr_ = (kf_g1_cpt - vtx).Mag();
      }
    }
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
