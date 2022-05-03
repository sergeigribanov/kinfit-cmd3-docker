//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 24 21:50:46 2019 by ROOT version 6.18/04
// from TTree tr_ph/Tree with the non-collinear events
// found on file: tr_ph_run021171.root
//////////////////////////////////////////////////////////

#ifndef _KFCmd_5pi_TrPh_H_
#define _KFCmd_5pi_TrPh_H_

#include <set>
#include <vector>
#include <kfcmd/core/TrPh.hpp>

class TrPh : public kfcmd::core::TrPh {
 public:
  TrPh(TTree* tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string& outpath,
                    double magneticField = 1.3) override final;

 private:
  bool cutTracks_();
  bool cutPhotons_();
  void setupOutputBranches_(TTree*);
  void fit_(kfcmd::hypos::Hypo4ChPions2Photons *);
  void fit_mpi0_(kfcmd::hypos::Hypo4ChPions2Photons *);
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  static const double minPhEn_;
  static const double maxPhEn_;
  static const std::set<std::string> sM3Pi_[4];
  static const std::set<std::string> sgg_;
  int kf_err_;
  double kf_chi2_;
  double in_m3pi_[4];
  double kf_m3pi_[4];
  double in_mgg_;
  double kf_mgg_;
  double vtx_[3];
  int kf_err_mpi0_;
  double kf_chi2_mpi0_;
  double in_m3pi_mpi0_[4];
  double kf_m3pi_mpi0_[4];
  double vtx_mpi0_[3];
  std::vector<std::size_t> trackIndices_;
  std::vector<std::size_t> photonIndices_;
};
#endif
