#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <set>

#include <TStopwatch.h>

#include "TrPh.h"
#include "Hypo2PiVertex.hpp"

TrPh::TrPh(TTree *tree) :
  kfcmd::core::TrPh(tree) {
}

TrPh::~TrPh() {
}

Int_t TrPh::Cut(Long64_t) {
  if (nt > 1) return 1;
  return 0;
}


void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("evnum", &evnum, "evnum/I");
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("niter", &niter_, "ntier/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("nt", &nt, "nt/I");
  tree->Branch("kf_mks", &kf_mks_, "kf_mks/D");
  tree->Branch("sim_is_ks", &sim_is_ks_, "sim_is_ks/I");
  tree->Branch("trph_is_ks", &trph_is_ks_, "trph_is_ks/I");
  tree->Branch("nks", &nks, "nks/I");
  tree->Branch("ksminv", ksminv, "ksminv[nks]/F");
}


void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  std::set<std::string> pions = {"pi-_1", "pi+_1"};
  kfcmd::hypos::Hypo2PiVertex hypo(2.e-3 * emeas, magneticField, 20, 1.e-5);
  Long64_t nentries = fChain->GetEntriesFast();
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (Cut(ientry) == 0) continue;
    trph_is_ks_ = (nks > 0);
    for (int k = 0; k < nsim; ++k) {
        sim_is_ks_ = 0;
        if (simorig[k] == 310) sim_is_ks_ = 1;
      }
    hypo.setBeamXY(xbeam, ybeam);
    int ch_i = 0;
    int ch_j = 0;
    kf_chi2_ = std::numeric_limits<double>::infinity();
    kf_err_ = 1;
    niter_ = -1;
     for (int i = 0; i < nt; ++i) {
        ch_i = tcharge[i];
         for (int j = 0; j < nt; ++j) {
           ch_j = tcharge[j];
             if (ch_i == ch_j) continue;
             if (ch_i < 0) {
                if(!hypo.fillTrack("pi-_1", i, *this)) continue;
                if(!hypo.fillTrack("pi+_1", j, *this)) continue;
             } else {
                hypo.fillTrack("pi-_1", j, *this);
                hypo.fillTrack("pi+_1", i, *this);
             }
             hypo.optimize();
             if (hypo.getErrorCode() == 0 && kf_chi2_ > hypo.getChiSquare()) {
                 kf_err_ = 0;
                 niter_ = hypo.getNumOfRequiredIters();
                 kf_chi2_ = hypo.getChiSquare();
                 kf_mks_ = hypo.getFinalMomentum(pions).M();
             }
        }
    }

    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  outfl->Close();
  delete outfl;
}
