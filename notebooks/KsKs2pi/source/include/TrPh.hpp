#ifndef _KsKs2ChPi_TrPh_HPP_
#define _KsKs2ChPi_TrPh_HPP_

#include <set>
#include <kfcmd/core/TrPh.hpp>

typedef struct {
  std::string ifname;
  std::string ofname;
  double mfield;
} CmdOptions;

class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t energy) override final;
  virtual void Loop(const std::string&, double mfield = 1.3) override final;
private:
  bool cutTracks_();
  void setupOutputBranches_(TTree*);
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  std::vector<std::size_t> trackIndices_;
  static const std::set<std::string> sKs1_;
  static const std::set<std::string> sKs2_;
  int kf_err_;
  double kf_chi2_;
  double in_mks1_;
  double in_mks2_;
  double kf_mks1_;
  double kf_mks2_;
  double vtx0_[3];
  double vtx1_[3];
  double vtx2_[3];
  double vtx1_dr_;
  double vtx2_dr_;
};

#endif
