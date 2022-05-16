#ifndef _KsKlPi0_TrPh_H_
#define _KsKlPi0_TrPh_H_

#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/HypoKsKlPi0.hpp>

class TrPh : public kfcmd::core::TrPh {
 public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string&, double magneticField = 1.3) override final;
 private:
  bool cutTracks_();
  bool cutPhotons_();
  double getAngle_(const kfcmd::hypos::HypoKsKlPi0 &hypo) const;
  void setupOutputBranches_(TTree *);
  static const std::set<std::string> decaypds_;
  static const std::set<std::string> gg_;
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  std::vector<std::size_t> trackIndices_;
  std::vector<std::size_t> photonIndices_;
  int kf_err_;
  double kf_chi2_;
  double in_mgg_;
  double in_mks_;
  double kf_mgg_;
  double kf_mks_;
  double kf_vtx0_[3];
  double kf_vtx1_[3];
  double kf_vtx_dr_;
  double kf_vtx_drho_;
  double kf_dvtx_2pi_angle_;
};

#endif
