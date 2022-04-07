#ifndef _TRPH_5PI_H_
#define _TRPH_5PI_H_
#include <string>
#include <vector>
#include <TRandom.h>
#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/Hypo4ChPions2Photons.hpp>

namespace kfcore = kfcmd::core;
namespace kfhypos = kfcmd::hypos;

typedef struct {
  Eigen::VectorXd parPiMi0;
  Eigen::VectorXd parPiMi1;
  Eigen::VectorXd parPiPl0;
  Eigen::VectorXd parPiPl1;
  Eigen::VectorXd parGamma0;
  Eigen::VectorXd parGamma1;
  Eigen::MatrixXd invCovPiMi0;
  Eigen::MatrixXd invCovPiMi1;
  Eigen::MatrixXd invCovPiPl0;
  Eigen::MatrixXd invCovPiPl1;
  Eigen::MatrixXd invCovGamma0;
  Eigen::MatrixXd invCovGamma1;
} InputParams;

class TrPh : public kfcore::TrPh {
public:
  explicit TrPh(TTree* = 0);
  virtual ~TrPh();
  virtual void Loop(const std::string&, double) override final;
  void setEntry(int);
  void setNEvents(int);
private:
  bool cutTracks_();
  bool cutPhotons_();
  void setupOutptuBranches_(TTree*);
  bool cut_();
  bool fitOnce_(kfhypos::Hypo4ChPions2Photons*);
  void fillParams_(const kfhypos::Hypo4ChPions2Photons &);
  void refit_(kfhypos::Hypo4ChPions2Photons*, TH1F*);
  int entry_;
  int nevents_;
  std::vector<std::size_t> trackIndices_;
  std::vector<std::size_t> photonIndices_;
  TRandom rnd_;
  InputParams params_;
};

#endif
