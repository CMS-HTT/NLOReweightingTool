import ConfigParser
import shelve
import copy
from array import array

from ROOT import gROOT, gStyle, TFile, TH1F, TF1, TCanvas, TLegend

gROOT.SetBatch(True)


def AvoidNegativeBin(hist):
    '''Set negative entries to zero'''
    for i in range (1,hist.GetNbinsX()+1) :
        if hist.GetBinContent(i) < 0:
            hist.SetBinContent(i, 0)

        if hist.GetBinContent(i) < 0 and abs(hist.GetBinContent(i)/hist.GetSumOfWeights()) > 0.001:
            print 'WARNING :', hist.GetName(), 'bin', i, 'is negative :', hist.GetBinContent(i), '/', hist.GetSumOfWeights(), ' -> set to zero'


def histCreator(file, var, mass, name):
    '''Simple function to derive PY8 spectrum'''

    tfile = TFile(file)
    tree = tfile.Get(var['tree'])
    
    hist = TH1F('h_' + name + '_' + str(mass),
                'h_' + name + '_' + str(mass),
                var['nbin'], var['xmin'], var['xmax'])
    
    hist.Sumw2()
    
    tree.Draw(var['var'] + ' >> ' + hist.GetName(), '(1)*weight')
    return copy.deepcopy(hist)




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
        self.var = {'tree':'tree', 'var':'gen_vpt', 'nbin':40, 'xmin':0, 'xmax':800}
        self.canvas = TCanvas('validation', 'validation', 1400,600)
        self.canvas.Print('validation.pdf[')

        print 
        print 'mass point :', self.mp
        print 'particle :', self.particle
        print 'tanb used for 2HDM :', self.tanb_2HDM
        print 'tanb used for PY8 :', self.tanb_PY8
        print 'ROOT file :', self.dir
        print         

    def overlay(self, hists, id):
        self.canvas.cd(id)
        ymax = max([hist.GetMaximum() for hist in hists])
        
        for ii, hist in enumerate(hists):
            hist.SetMaximum(ymax*1.4)
            hist.SetLineWidth(2)
            hist.SetLineColor(ii+1)

            if ii==0: hist.Draw()
            else: hist.Draw('same')


    def derive(self):

        results = []

        for mass in self.mp:

            # 1st step : derive PY8 spectrum

            hist_PY8 = histCreator(self.dir + '/GEN_testrun-lhc-A-mA' + str(mass) + '_tb' + self.tanb_PY8 + '_t_PY8/all.root', self.var, mass, 'PY8')

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

                hist = histCreator(val['file'], self.var, mass, key)
                hists_2hdm.append(copy.deepcopy(hist))


            interference = copy.deepcopy(hists_2hdm[0])
            interference.Add(hists_2hdm[1], -1)
            interference.Add(hists_2hdm[2], -1)


            Yt_2HDM = 1./float(self.tanb_2HDM)
            Yb_2HDM = float(self.tanb_2HDM)


            # 3rd step : rescaling 2HDM spectrum to MSSM

            for tanb, val in sorted(self.Yukawa.iteritems()):

                if tanb.find('0.')!=-1: continue
                
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
                hist_mssm.SetName('MSSM_' + str(mass) + '_' + tanb)
                AvoidNegativeBin(hist_mssm)

                hists =[hist_PY8, hist_mssm]

                for ii, ihist in enumerate(hists):
                    ihist.Scale(1./ihist.GetSumOfWeights())

                ratio = copy.deepcopy(hists[1])
                ratio.Divide(hists[1], hists[0],1,1,'b')
                
                ratio.GetYaxis().SetRangeUser(ratio.GetBinContent(ratio.GetMinimumBin())*0.8, ratio.GetBinContent(ratio.GetMaximumBin())*1.2)

                myfunc = TF1('myfunc', 'expo(0) + pol3(2)')

                ratio.Fit('myfunc','Q')
                func = ratio.GetFunction('myfunc')
                func.SetName('weight_mA' + str(mass) + '_' + tanb)
                
                if ratio.GetMinimum() < 0:
                    print 'Warning !'

                results.append(copy.deepcopy(func))

                # validation plots:
                
                leg = '(mA, tanb) = (' + str(mass) + ', ' + tanb.replace('tanb_','') + ')' 
                
                self.canvas.Clear()
                self.canvas.Divide(2)
                self.overlay(hists, 1)
                self.overlay([ratio], 2)
                self.canvas.Print('validation.pdf')


        ofile = TFile('../user/Reweight.root', 'recreate')
        for ii in results:
            ii.GetXaxis().SetTitle('Generated Higgs pT (GeV)')
            ii.GetYaxis().SetTitle('(NLO/PY8) weight')
            ii.Write()
        ofile.Write()
        ofile.Close()

        self.canvas.Print('validation.pdf]')

if __name__ == '__main__':
    tool = ReweightingManager()
    tool.derive()

    print 
    print 'Done !'
