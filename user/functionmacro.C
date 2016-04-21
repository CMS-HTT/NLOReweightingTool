/*
  functionmacro.C : provide Higgs pT reweighting factor for the given mass and tan(beta)
  20 April 2016
  Yuta Takahashi
 */

#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include "TROOT.h"

TF1* func[31][65];
std::vector<int> marray = {80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200};

void ReadFile(){

  TFile *file = new TFile("/afs/cern.ch/user/y/ytakahas/public/forHTT/Reweight.root");

  const int num_of_tb = 60;  

  int imass = 0;
  for(auto mass: marray){
    for(int tanb=0; tanb < num_of_tb; tanb++){

      TString wname = "weight_mA";
      wname += mass;
      wname += "_tanb_";
      wname += tanb + 1;

      func[imass][tanb] = (TF1*) gROOT->FindObject(wname);
    }    
    imass++;
  }
  delete file;
}

float returnNLOweight(int mass, int tanb, float pt){

  if(pt > 800){
    //    std::cout << "[INFO] pT = " << pt << " exceeds the range --> set it to 800." << std::endl;    
    pt = 800;
  }

  auto iter = std::find(marray.begin(), marray.end(), mass);
  size_t index = std::distance(marray.begin(), iter);

  if(index == marray.size()){
    std::cout << "[WARNING] Invalid mass point ... " << mass << " -> return weight 1" << std::endl;    
    return 1;
  }
  
  if(tanb <1 || tanb > 60){
    std::cout << "[WARNING] Invalid tan(beta) point ... " << tanb << " -> return weight 1" << std::endl;
    return 1;
  }

  return func[index][tanb-1]->Eval(pt);

}

void functionmacro(){
  std::cout << std::endl;
  std::cout << "Initialize functionmacro.C ... " << std::endl;
  std::cout << std::endl;
  ReadFile();
}
