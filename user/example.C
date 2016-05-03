void example(){

  // compile functionmacro.C and create shared library
  // replace directory, as you wish (i.e. if you want to work on different directory)
  gROOT->Macro("./functionmacro.C+g");

  // read your ROOT file ... replace as you wish
  TFile *tfile = new TFile("/afs/cern.ch/user/y/ytakahas/public/forHTT/tree_HiggsSUSYGG500.root");
  TTree *tree = (TTree*) tfile->Get("tree");

  // define histogram
  TH1F* h_PY8 = new TH1F("h_PY8", "h_PY8", 30, 0, 500);
  TH1F* h_PY8_weighted = new TH1F("h_PY8_weighted", "h_PY8_weighted", 30, 0, 500);
   
  h_PY8->Sumw2();
  h_PY8_weighted->Sumw2();

  /*
    create histogram with Higgs pT, with additional weights from NLO:
    the function, returnNLOweight(mass, tanb, pthiggs) will provide additional weight, as shown below example.
    
    mass : Choose exactly the same mass point as the produced PY8 sample (in this example, 500 GeV)
           The following mass points are supported
	   80 90 100 110 120 130 140 160 180 200 250 300 350 400 450 500 600 700 800 900 1000 1200 1400 1500 1600 1800 2000

    tanb : Choose tanb from 1 ... 60

    pthiggs : reweighting currently scoping 0 - 800 GeV.
              pT > 800 GeV will be forced to be 800 GeV.   

    The example below provides additional weights for 500 GeV, tanb = 10.
    Change variable name as you wish (in this example, pthiggs is "the" name)
    You can also draw whatever variables with the new weight
  */ 

  tree->Draw("pthiggs >> h_PY8_weighted", "weight*returnNLOweight(500, 10, pthiggs)");

  // This is without Higgs pT reweighting
  tree->Draw("pthiggs >> h_PY8", "weight");


  // Draw
  h_PY8_weighted->SetLineColor(kRed);
  h_PY8->SetLineColor(kBlue);

  h_PY8_weighted->Draw();
  h_PY8->Draw("same");

  TLegend* leg = new TLegend(0.5,0.6,0.7,0.8);
  leg->AddEntry(h_PY8, "PY8", "lep");
  leg->AddEntry(h_PY8_weighted, "PY8 (NLO weighted)", "lep");

  leg->Draw();
}
