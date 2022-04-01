#ifndef __KsChKaonChPion_TrPh__
#define __KsChKaonChPion_TrPh__

#include <kfcmd/core/TrPh.hpp>

class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string&, double magneticField = 1.3) override final;
private:
  bool cutTracks_();
  void setupOutptuBranches_(TTree*);
  void printStatus_(int, int);
  void printSummary_(int*, int, int);
  int getStatus_() const;
  void setStatus_(int);
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  int status_;
  double kf_mks_;
  double in_mks_;
  double kf_chi2_;
  double kf_vtx0_x_;
  double kf_vtx0_y_;
  double kf_vtx0_z_;
  double kf_vtx1_x_;
  double kf_vtx1_y_;
  double kf_vtx1_z_;
  double kf_vtx_dr_;
  double kf_vtx_drho_;
  double kf_dedx_vtx0_K_;
  double kf_dedx_vtx0_pi_;
  double kf_dedx_vtx1_pi_[2];
  double kf_p_vtx0_K_;
  double kf_p_vtx0_pi_;
  double kf_p_vtx1_pi_[2];
  TStopwatch time_;
  TStopwatch timePerEntry_;
  std::vector<std::size_t> trackIndices_;
};

#endif
