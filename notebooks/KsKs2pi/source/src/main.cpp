#include "TrPh.hpp"
#include <TFile.h>
#include <TTree.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
namespace po = boost::program_options;

void setOptions(po::options_description *desc, CmdOptions *opts) {
  desc->add_options()
    ("help,h", "CMD-3 kinematic fit for hypothises: e+e- --> pi+pi-X")
    ("ifname,i", po::value<std::string>(&(opts->ifname)), "Path to input .root file.")
    ("ofname,o", po::value<std::string>(&(opts->ofname))->default_value("test.root"),
     "Path to output .root file.")
    ("mfield,m", po::value<double>(&(opts->mfield))->default_value(1.3), "Magnetic field");
}

void help(const po::options_description &desc) {
  std::cout << desc << std::endl;
}

template <class T> T *find_object(TFile *fl, const std::string &object_name) {
  T *object = dynamic_cast<T *>(fl->Get(object_name.c_str()));
  if (!object) {
    std::cerr << "[!] Object with name \"" << object_name << "\" of class "
              << T::Class_Name() << " is not found in file " << fl->GetName()
              << std::endl;
    exit(1);
  }
  return object;
}

int main(int argc, char *argv[]) {
  po::options_description desc("Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (!vmap.count("ifname")) {
    help(desc);
    return 1;
  }
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  auto tree = find_object<TTree>(fl, "tr_ph");
  TrPh fitter(tree);
  fitter.Loop(opts.ofname, opts.mfield);
  fl->Close();
  return 0;
}
