#ifndef _KFCmd_3Gamma_TrPh_
#define _KFCmd_3Gamma_TrPh_

#include <set>
#include <string>
#include <vector>
#include <kfcmd/core/TrPh.hpp>
#include "Hypo3PhotonsCustom.hpp"

class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t energy) override final;
  virtual void Loop(const std::string&, double mfield = 1.3) override final;
private:
  bool cutPhotons_();
  void setupOutputBranches_(TTree *);
  void fit_(Hypo3PhotonsCustom*);
  std::vector<std::size_t> photonIndices_;
  static const double min_energy_;
  static const std::set<std::string> s_phpair0_;
  static const std::set<std::string> s_phpair1_;
  static const std::set<std::string> s_phpair2_;
  static const std::set<std::string> s_all_photons_;
  double kf_err_;
  double kf_chi2_;
  double in_mgg_[3];
  double kf_mgg_[3];
  double kf_vtx_x_;
  double kf_vtx_y_;
  double kf_vtx_z_;
  double in_total_px_;
  double in_total_py_;
  double in_total_pz_;
  double in_total_pe_;
  double kf_total_px_;
  double kf_total_py_;
  double kf_total_pz_;
  double kf_total_pe_;
};

#endif
