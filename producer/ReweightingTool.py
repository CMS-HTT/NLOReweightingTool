import ConfigParser
import shelve
import copy
from ROOT import gROOT, gStyle, TFile, TH1F

gROOT.SetBatch(True)


def histCreator_PY8(dir, tanb_PY8, var, mass):
    '''Simple function to derive PY8 spectrum'''

    tfile = TFile(dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + tanb_PY8 + '_t_PY8/all.root')
    tree = tfile.Get(var['tree'])
    
    hist_py8 = TH1F('h_PY8_' + str(mass),
                    'h_PY8_' + str(mass),
                    var['nbin'], var['xmin'], var['xmax'])
    
    hist_py8.Sumw2()
    
    tree.Draw(var['var'] + ' >> ' + hist_py8.GetName(), '(1)*weight')
    
    return copy.deepcopy(hist_py8)



class ReweightingManager(object):
    def __init__(self):
        
        init = ConfigParser.SafeConfigParser()
        init.read("./config.ini")

        self.mp = init.get('settings', 'mp').split()
        self.particle = init.get('settings', 'particle')
        self.dir = init.get('settings', 'dir')
        self.tanb_2HDM = init.get('settings', 'tanb_2HDM')
        self.tanb_PY8 = init.get('settings', 'tanb_PY8')
        self.Yukawa = shelve.open('Yukawa_' + self.particle + '.db')['Yukawa']
        
        # variables for the reweighting
        self.var = {'tree':'tree', 'var':'gen_vpt', 'nbin':50, 'xmin':0, 'xmax':700}


        print 
        print 'mass point :', self.mp
        print 'particle :', self.particle
        print 'tanb used for 2HDM :', self.tanb_2HDM
        print 'tanb used for PY8 :', self.tanb_PY8
        print 'ROOT file :', self.dir
        print 



    def derive(self):

        results = []

        for mass in self.mp:

            hist_PY8 = histCreator_PY8(self.dir, self.tanb_PY8, self.var, mass)

            # 2nd step : derive NLO 2HDM spectrum

            hdict = {
                'A':{'file':self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_2HDM + '_tb_interference/all.root'},
                'B':{'file':self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_2HDM + '_t_interference/all.root'},
                'C':{'file':self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_2HDM + '_b_interference/all.root'},
                'D':{'file':self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_2HDM + '_t_pure/all.root'},
                'E':{'file':self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_2HDM + '_b_pure/all.root'}
                }

            hists_2hdm = []

            for key, val in sorted(hdict.iteritems()):

                tfile = TFile(val['file'])
                tree = tfile.Get(self.var['tree'])
                
                hist = TH1F('h_' + key + '_' + str(mass), 
                            'h_' + key + '_' + str(mass), 
                            self.var['nbin'], self.var['xmin'], self.var['xmax'])
                hist.Sumw2()

                tree.Draw(self.var['var'] + ' >> ' + hist.GetName(), '(1)*weight')
                hists_2hdm.append(copy.deepcopy(hist))


            interference = copy.deepcopy(hists_2hdm[0])
            interference.Add(hists_2hdm[1], -1)
            interference.Add(hists_2hdm[2], -1)


            Yt_2HDM = 1./float(self.tanb_2HDM)
            Yb_2HDM = float(self.tanb_2HDM)


            # 3rd step : rescaling 2HDM spectrum to MSSM

            for tanb, val in sorted(self.Yukawa.iteritems()):
                
                print 'processing ... ', mass, tanb

                Yt_MSSM = val['Yt']
                Yb_MSSM = val['Yb']
                        
                hist_mssm = copy.deepcopy(interference)
                hist_mssm.Scale(Yt_MSSM*Yb_MSSM/(Yt_2HDM*Yb_2HDM))
                
                rescaled_top = copy.deepcopy(hists_2hdm[3])
                rescaled_top.Scale(Yt_MSSM*Yt_MSSM/(Yt_2HDM*Yt_2HDM))
                hist_mssm.Add(rescaled_top)

                rescaled_bottom = copy.deepcopy(hists_2hdm[4])
                rescaled_bottom.Scale(Yb_MSSM*Yb_MSSM/(Yb_2HDM*Yb_2HDM))
                hist_mssm.Add(rescaled_bottom)

                hist_mssm.GetYaxis().SetRangeUser(0., hist_mssm.GetBinContent(hist_mssm.GetMaximumBin())*1.1)

                hists =[hist_PY8, hist_mssm]

                for ii, ihist in enumerate(hists):
                    ihist.Scale(1./ihist.GetSumOfWeights())

                ratio = copy.deepcopy(hists[1])
                ratio.Divide(hists[0])

                ratio.GetYaxis().SetRangeUser(ratio.GetBinContent(ratio.GetMinimumBin())*0.8, ratio.GetBinContent(ratio.GetMaximumBin())*1.2)
                ratio.Draw()
                ratio.Fit('pol6')
                func = ratio.GetFunction('pol6')
                func.SetName('weight_mA' + str(mass) + '_' + tanb)

                results.append(copy.deepcopy(func))



        print 'Writing to files ...'
        ofile = TFile('../user/Reweight.root', 'recreate')
        for ii in results:
            print 'Writing ...', ii.GetName()
            ii.GetXaxis().SetTitle('Generated Higgs pT (GeV)')
            ii.GetYaxis().SetTitle('(NLO/PY8) weight')
            ii.Write()
        ofile.Write()
        ofile.Close()

if __name__ == '__main__':
    tool = ReweightingManager()
    tool.derive()

    print 
    print 'Done !'
