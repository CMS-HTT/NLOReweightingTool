from ROOT import gSystem, TFile, TH1F


if __name__ == '__main__':

    gSystem.Load('./functionmacro_C.so')
    from ROOT import returnNLOweight

    tfile = TFile('../../gen_analysis/GEN_testrun-lhc-A-mA500_tb30_t_PY8/all.root')
    tree = tfile.Get('tree')
    
    hist_py8 = TH1F('h_PY8_500',
                    'h_PY8_500',
                    30, 0, 500)
    
    hist_py8.Sumw2()
 
    import pdb; pdb.set_trace()
    tree.Draw('gen_vpt >> h_PY8_500', '(1)*weight')



#    tree.Draw('gen_vpt >> h_PY8_500', '(1)*weight*returnNLOweight(500, 5, gen_vpt)')
