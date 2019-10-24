#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>       

// create N = n_subsamples based on the event number, to be used to train/apply the BDT 

int n_subsamples = 11;

void create_subsample(int year = 2016, int type = 0, int mc = 0){
  
  // define input files based on the year and the dataset
  std::string version = "Feb5_skimSoftMu";
  if      (year == 2016)  version = "Oct2019"; //"Feb5_skimSoftMu"; 
  else if (year == 2017)  version = "Apr16_fixIso"; 
  else if (year == 2018)  version = "Apr30"; 
  else {
    std::cout << "please enter valid year" << std::endl;
    return;
  }
  std::string     channel = "LMNR";
  if (type == 1)  channel = "Charmonium"; 

  std::string     isdata = "data";
  if (mc == 1)    isdata = "MC"; 

  std::cout << "creating sub_samples for year: " << year << "  channel: " << channel << std::endl;

  TChain* tree = new TChain("ntuple");
  if (mc ==1)
      tree -> Add(Form("/gwteray/users/fiorendi/p5prime/data%d/flat_ntuples/%s/%dMC_%s_part*.root", year, version.c_str(), year, channel.c_str() ));
  else
      tree -> Add(Form("/gwteray/users/fiorendi/p5prime/data%d/flat_ntuples/%s/ntuple_%s_%d*.root", year, version.c_str(), channel.c_str(), year ));

  // define output files
  std::vector<TFile*> out_files;
  for (int i=0; i < n_subsamples; i++){
      TFile* ifile = new TFile(Form("sub_samples/sample_%d_%s_%s_%d.root", year, isdata.c_str(), channel.c_str(), i),"RECREATE");
      out_files.push_back(ifile);
      std::cout << "output file : " << out_files.at(i) -> GetName() << std::endl;
  }
  
  tree -> SetBranchStatus("*",1);
  long eventN;
  tree -> SetBranchAddress( "eventN", &eventN);
  int nentries = tree->GetEntries();
  std::cout << "original n entries: " << nentries << std::endl;

  std::vector<TTree*> out_trees;
  std::vector<bool> keepentries;
  for (int i = 0; i < n_subsamples; i++){
      out_files.at(i) -> cd();
      TTree* itree = (TTree*) tree->CloneTree(0);
      out_trees.push_back(itree);
      out_trees.at(i)->SetAutoSave(1000000);
      keepentries.push_back(false);
  }

  double q   = 1.;
  double mod;
  double fractpart, intpart;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent = tree-> GetEntry(eventNo);
    q = fabs(eventN) / n_subsamples;
    fractpart = modf (q , &intpart);
    mod = fabs(eventN) - n_subsamples*intpart;
    out_trees.at(mod) -> Fill();
  }

  for (int i = 0; i < n_subsamples; i++){
      out_files.at(i) -> cd();
      out_trees.at(i) -> Write();
  }  
  return;
}