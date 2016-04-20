from ROOT import gSystem, TFile, TH1F, gROOT

if __name__ == '__main__':

    gROOT.Macro('./functionmacro.C+g')
    from ROOT import returnNLOweight

    tfile = TFile('../../gen_analysis/GEN_testrun-lhc-A-mA500_tb30_t_PY8/all.root')
    tree = tfile.Get('tree')
    
    hist_py8 = TH1F('h_PY8', 'h_PY8', 30, 0, 500)
    hist_py8_weighted = TH1F('h_PY8_weighted', 'h_PY8_weighted', 30, 0, 500)
    
    hist_py8.Sumw2()
    hist_py8_weighted.Sumw2()


    tree.Draw('gen_vpt >> h_PY8_weighted', '(1)*weight*returnNLOweight(500, 10, gen_vpt)')
    tree.Draw('gen_vpt >> h_PY8', '(1)*weight', 'same')
