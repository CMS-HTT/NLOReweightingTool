using namespace std;

#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include "TROOT.h"
#include <typeinfo>


TF1* func[31][65];
std::vector<int> marray;

void ReadFile(){

  TFile *file = new TFile("Reweight.root");

  //  Int_t mp[] = {80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200};
  Int_t mp[] = {500, 800};
  const int num_of_mp = sizeof mp /sizeof mp[0];

  for(int mass=0; mass < num_of_mp; mass++){
    marray.push_back(mp[mass]);
  }

  const int num_of_tb = 60;  


  for(int mass=0; mass < (int)marray.size(); mass++){
    for(int tanb=0; tanb < num_of_tb; tanb++){

      TString wname = "weight_mA";
      wname += mp[mass];
      wname += "_tanb_";
      wname += tanb + 1;

      func[mass][tanb] = (TF1*) gROOT->FindObject(wname);
    }    
  }
  
  delete file;
}

float returnNLOweight(int mass, int tanb, float pt){

  if(pt > 700) pt = 700;

  if(typeid(mass) != typeid(int)) std::cout << "mass variable is not <int>" << std::endl;
  if(typeid(tanb) != typeid(int)) std::cout << "tanb variable is not <int>" << std::endl;
  if(typeid(pt) != typeid(float)) std::cout << "pt variable is not <float>" << std::endl;

  std::vector<int>::iterator iter = std::find(marray.begin(), marray.end(), mass);
  size_t index = std::distance(marray.begin(), iter);

  std::cout << marray.size() << std::endl;
  
  for(int ii=0; ii < (int)marray.size(); ii++){
    std::cout << marray.at(ii) << " / " << marray.size() << " " << index << std::endl;
  }

  if(index == marray.size()){
    //    std::cout << "Invalid mass point ... " << mass << " -> return weight 1" << std::endl;    
    return 1;
  }
  
  if(tanb <1 || tanb > 60){
    //    std::cout << "Invalid tan(beta) point ... " << tanb << " -> return weight 1" << std::endl;
    return 1;
  }

  //  std::cout << "mass index = " << index << ", tanb index = " << tanb-1 << ", pt = " << pt << std::endl;
  
  return func[index][tanb-1]->Eval(pt);

}

void functionmacro(){
  std::cout << "reading function macro ..." << std::endl;
  ReadFile();
}
