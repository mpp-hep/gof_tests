# Imports
import os, sys
import numpy as np
import ROOT as rt
import uproot
import h5py
import matplotlib.pyplot as plt
import root_numpy as rtnp
import pandas

from math import log10, floor
from subprocess import Popen
import shlex 
import os, sys
import time


#Color style
scol = ['#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000']
scol.reverse()
bcol = ['#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#005824'] 
bcol.reverse()

sig_xsec       = 10. #In units of fb (== 0.01 pb) 
sig_inj_vals = [0,0.1, 1, 2,5,10] #These correspond to xsecs of 0fb(used for GOF w.o signal), 1fb, 10fb, 20fb, 50fb, 100fb
quants  = [0.01,0.1,0.3,0.5,0.7,0.9]
sig_inj_vals = [0]

def makeWS(data_rej,data_acc,signal_rej,signal_acc, outname, quantile,injectedSignal):
  
  efficiency = data_acc.Integral()/data_rej.Integral() #How much rej histogram needs to be scaled to match acc
  print('Using efficiency = {} '.format(efficiency))
  
  xaxis = data_acc.GetXaxis().GetXbins()
  min_bin = xaxis[0]
  max_bin = xaxis[len(xaxis)-1]
  n_bins  = len(xaxis)-1
  # set up workspace
  datacard_ws = rt.TFile.Open(outname,'recreate')
  w = rt.RooWorkspace('w','w')
  x = rt.RooRealVar('x','x',min_bin,max_bin)
  w.factory('x[%.1f,%.1f]'%(min_bin, max_bin))

  acc_bin_functions = rt.RooArgList()
  rej_bin_functions = rt.RooArgList()
  w.factory('eff[%f,-100.,100.]'%(efficiency))
  w.var('eff').setConstant(False)

  empty_hist = rt.TH1D('empty_hist','empty_hist', data_rej.GetNbinsX(), data_rej.GetXaxis().GetXbins().GetArray())
  # empty_hist = rt.TH1D('empty_hist','empty_hist', n_bins, min_bin, max_bin)
  for iBinX in range(1,data_rej.GetNbinsX()+1):
      empty_hist.SetBinContent(iBinX,1)
      w.factory('crBin%i_In[%.1f]'%(iBinX,data_rej.GetBinContent(iBinX)))
      w.factory('crBin%i[0,-100,100]'%(iBinX))
      w.var('crBin%i_In'%iBinX).setConstant(True)
      w.var('crBin%i'%iBinX).setConstant(False)
      if data_rej.GetBinContent(iBinX) !=  0.:
        power = 1/rt.TMath.Sqrt(data_rej.GetBinContent(iBinX))
      else:
        power = -1.0
        w.var('crBin%i'%iBinX).setConstant(True)
      #what is fit is actually (mjj+mjj/sqrt(mjj))^x_bin, meaning x_bin is representative of how many sigma acc is from rej
      w.factory("expr::crBin%iFunc('max(0,@0*pow(1.0+%f,@1))',crBin%i_In,crBin%i)"%(iBinX,power,iBinX,iBinX))
      # w.factory("expr::crBin%iFunc('max(0,@0*pow(1.0+%f,@1))',crBin%i_In,crBin%i)"%(iBinX,1/rt.TMath.Sqrt(data_rej.GetBinContent(iBinX)),iBinX,iBinX))
      w.factory("expr::bin%iFunc('max(0,@0*@1)',eff,crBin%iFunc)"%(iBinX,iBinX))
      rej_bin_functions.add(w.function('crBin%iFunc'%iBinX))
      acc_bin_functions.add(w.function('bin%iFunc'%iBinX))
    
  qcd_rph_rej = rt.RooParametricHist('background_rej','background_rej',w.var('x'),rej_bin_functions,empty_hist)
  qcd_rph_rej_norm = rt.RooAddition('background_rej_norm','background_rej_norm',rej_bin_functions)
  qcd_rph_acc = rt.RooParametricHist('background_acc','background_acc',w.var('x'),acc_bin_functions,empty_hist)
  qcd_rph_acc_norm = rt.RooAddition('background_acc_norm','background_acc_norm',acc_bin_functions)
  getattr(w,'import')(qcd_rph_rej, rt.RooCmdArg())
  getattr(w,'import')(qcd_rph_rej_norm, rt.RooFit.RecycleConflictNodes())
  getattr(w,'import')(qcd_rph_acc, rt.RooCmdArg())
  getattr(w,'import')(qcd_rph_acc_norm, rt.RooFit.RecycleConflictNodes())

  ds_signal_acc = rt.RooDataHist('signal_acc','signal_acc',rt.RooArgList(w.var('x')),signal_acc)
  ds_signal_rej = rt.RooDataHist('signal_rej','signal_rej',rt.RooArgList(w.var('x')),signal_rej)
  getattr(w,'import')(ds_signal_acc, rt.RooCmdArg())
  getattr(w,'import')(ds_signal_rej, rt.RooCmdArg())

  ds_data_acc = rt.RooDataHist('data_obs_acc','data_obs_acc',rt.RooArgList(w.var('x')),data_acc)
  ds_data_rej = rt.RooDataHist('data_obs_rej','data_obs_rej',rt.RooArgList(w.var('x')),data_rej)
  getattr(w,'import')(ds_data_acc, rt.RooCmdArg())
  getattr(w,'import')(ds_data_rej, rt.RooCmdArg())

  datacard_ws.cd()
  w.Write()
  datacard_ws.Close()

  # w.Print('v')
  
  datacard_ratio = \
  '''
  imax 1
  jmax 1
  kmax *
  ---------------
  shapes * * {WS} w:$PROCESS_$CHANNEL w:$PROCESS_$CHANNEL_$SYSTEMATIC
  ---------------
  bin {BIN}
  observation {OBS}
  ------------------------------
  bin             {BIN}      {BIN}
  process         signal     background
  process         0          1
  rate            {SIGRATE}    {BKGRATE}
  --------------------------------
  lumi lnN 1.01 -
  eff   flatParam
  '''

  for i in range(1,n_bins+1):
      datacard_ratio += 'crBin%i   flatParam\n'%i
  
  # write datacard
  datacard_ratio_acc = datacard_ratio.format(BIN='acc',
                            OBS=-1,
                            BKGRATE=1,
                            SIGRATE=signal_acc.Integral(),
                            WS=outname)
  print(datacard_ratio_acc)
  with open(outname.replace('.root','_acc.txt'),'w') as f:
      f.write(datacard_ratio_acc)
    
    
  datacard_ratio_rej = datacard_ratio.format(BIN='rej',
                            OBS=-1,
                            BKGRATE=1,
                            SIGRATE=signal_rej.Integral(),
                            WS=outname)
  print(datacard_ratio_rej)
  with open(outname.replace('.root','_rej.txt'),'w') as f:
      f.write(datacard_ratio_rej)
  os.system('combineCards.py rej={REJ} acc={ACC} > {RATIO}'.format(REJ=outname.replace('.root','_rej.txt'),ACC=outname.replace('.root','_acc.txt'),RATIO=outname.replace('.root','_ratio.txt')))
  
if __name__ == "__main__":
  
  # histfile = 'histograms_100GeV_5TeV.root'
  histfile = 'histograms_200GeV_5600.root'
  lumi = 60
  MX   = 3500. #Signal mass in units of GeV
  # sig_inj_vals = [0] #These correspond to xsecs of 0fb(used for GOF w.o signal), 1fb, 10fb, 20fb, 50fb, 100fb
  # quants       = [0.01]
  
  doCombine = True
  doPlots = False
  
  quants  = [0.1]
  
  
  for i,qkey in enumerate(quants):
  
    if i > 0:
      lower_bound = quants[i-1]
    else:
      lower_bound = 0.0
  
    for skey in sig_inj_vals:
      injectedSignal = int(skey*10.) #Stick to no injected signal for now
      quantile       = int(qkey*100)
      prefix = 'inj{}fb_q{}percent'.format(injectedSignal,quantile) #different name for all quantiles/injected signal

      # Get histograms
      f = rt.TFile.Open(histfile,"r")
      data_rej = f.Get('data_rej_{}'.format(prefix)); data_rej.SetDirectory(0);
      data_acc = f.Get('data_acc_{}'.format(prefix)); data_acc.SetDirectory(0);
      efficiency = data_acc.Integral()/data_rej.Integral()
    
      signal_acc     = f.Get('signal_q{}percent'.format(quantile)) ; signal_acc.SetDirectory(0);
      signal_rej     = f.Get('signal_template') ; signal_rej.SetDirectory(0);
      f.Close() 
    
      xaxis = data_acc.GetXaxis().GetXbins()
      min_bin = xaxis[0]
      max_bin = xaxis[len(xaxis)-1]
      n_bins  = len(xaxis)-1

      makeWS(data_rej,data_acc,signal_rej,signal_acc, outname='datacard_ws_{PREFIX}.root'.format(PREFIX=prefix), quantile=quantile,injectedSignal=injectedSignal)
      if doCombine:
        os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d datacard_ws_{PREFIX}_ratio.txt -n Ratio_gof_{PREFIX} --dataset data_obs -v 0'.format(PREFIX=prefix))
        for i in range(1):
          os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d datacard_ws_{PREFIX}_ratio.txt -t 200 --toysFreq -n Ratio_gof_toys_{PREFIX}  --dataset data_obs -s {S} -v 0'.format(PREFIX=prefix,S=40+i))
        os.system('hadd -f higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.ALLTOYS.root higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.4*.root'.format(PREFIX=prefix))
    
        # open file
        obs_gof_file = uproot.open('higgsCombineRatio_gof_{PREFIX}.GoodnessOfFit.mH120.root'.format(PREFIX=prefix))
        obs_gof = obs_gof_file['limit'].arrays('limit')['limit'][0]
        print("Obs gof = {}".format(obs_gof))
        exp_gof_file = uproot.open('higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.ALLTOYS.root'.format(PREFIX=prefix))
        exp_gof = exp_gof_file['limit'].arrays('limit')['limit']
        print("Exp gof (mean) = {}".format(np.mean(exp_gof)))
    
    
      if doPlots:

        # open file
        obs_gof_file = uproot.open('higgsCombineRatio_gof_{PREFIX}.GoodnessOfFit.mH120.root'.format(PREFIX=prefix))
        obs_gof = obs_gof_file['limit'].arrays('limit')['limit'][0]

        exp_gof_file = uproot.open('higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.ALLTOYS.root'.format(PREFIX=prefix)) #limit contains the value of the test-statistic in each toy
        exp_gof = exp_gof_file['limit'].arrays('limit')['limit']

        # get p-value from toys
        n_extreme = len(exp_gof[exp_gof > obs_gof])
        n_total = len(exp_gof)
        pval_toys = 1.*n_extreme/n_total

        # get p-value assuming chi2 dist (may not be valid)
        pval = rt.TMath.Prob(obs_gof,n_bins)
        # print('sig inj = %.1f fb, gof = %.2f, p-value (from chi2) = %.6f, p-value (from toys) = %.6f'%(sig_inj*sig_xsec,obs_gof,pval,pval_toys))

        bin_width = (max(exp_gof+[obs_gof])+np.std(exp_gof)-(min(exp_gof)-np.std(exp_gof)))/30.
        exp_gof_hist = rt.TH1D('gof','gof',30,min(exp_gof)-np.std(exp_gof), max(exp_gof+[obs_gof])+np.std(exp_gof))
        exp_gof_hist_gt = rt.TH1D('gof_gt','gof_gt',30,min(exp_gof)-np.std(exp_gof), max(exp_gof+[obs_gof])+np.std(exp_gof))
        for g in exp_gof:
            exp_gof_hist.Fill(g)
            if g > obs_gof:
                exp_gof_hist_gt.Fill(g)


        d = rt.TCanvas("ratio", "", 1000, 800)
        d.SetLeftMargin(0.13)
        # signal_hist_template    .Draw('same HIST')

        rt.gStyle.SetOptTitle(0)
        rt.gStyle.SetOptStat(0)
        f = rt.TF1("chi2","%f*ROOT::Math::chisquared_pdf(x,%i,0)"%(exp_gof_hist.Integral()*bin_width,n_bins),min(exp_gof)-np.std(exp_gof),max(exp_gof+[obs_gof])+np.std(exp_gof))

        tleg = rt.TLegend(0.48, 0.6, 0.89, 0.85)
        tleg.SetTextSize(0.05)
        tleg.SetBorderSize(0)
        tleg.SetFillStyle(0)
        tleg.SetTextSize (0.03)
        tleg.SetTextFont( 62 )
        tleg.SetTextSize (0.03)
        tleg.SetTextFont( 42 )
        exp_gof_hist.Draw('hist')
        exp_gof_hist.SetXTitle('Test statistic -2ln#lambda')
        exp_gof_hist.SetYTitle('N toys')
        exp_gof_hist.SetTitle("")
        exp_gof_hist.GetYaxis().SetLabelSize(0.05)
        exp_gof_hist.GetYaxis().SetTitleSize(0.05)
        f.SetLineColor((rt.TColor.GetColor(scol[1])))
        exp_gof_hist.SetLineWidth(2)
        exp_gof_hist.SetLineColor((rt.TColor.GetColor(bcol[0])))
        exp_gof_hist_gt.SetLineColor((rt.TColor.GetColor(bcol[0])))
        exp_gof_hist_gt.SetFillColorAlpha((rt.TColor.GetColor(bcol[0])), 0.30)
        exp_gof_hist_gt.Draw('fhistsame')
        f.Draw('same')
        line = rt.TLine(obs_gof,0,obs_gof,exp_gof_hist.GetMaximum())
        line.SetLineWidth(2)
        line.Draw()
        tleg.AddEntry(exp_gof_hist_gt,'p-value (from toys) = %.2f'%pval_toys)
        tleg.AddEntry(f,'p-value (from #chi^{2}) = %.2f'%pval,'l')
        tleg.AddEntry(line,'Observed (Best fit = {:.1f})'.format(obs_gof),'l')
        tleg.Draw()
        latex = rt.TLatex()
        latex.SetNDC ()
        latex.SetTextSize (0.03)
        latex.SetTextFont( 62 )
        latex.DrawLatex (0.67 ,0.27 , "%.1f fb^{-1} dijet events"%lumi)
        latex.SetTextSize (0.03)
        latex.SetTextFont( 42 )
        latex.DrawLatex(0.67 ,0.23 , 'q = {} - {}'.format(lower_bound,qkey))
        # latex.DrawLatex(0.67 ,0.19 , "Acc./rej. = {:.2f}".format(efficiency))
        # latex.DrawLatex(0.67 ,0.15 , "50 GeV binning")
        # d.Draw()
        d.SaveAs('gof_{}.pdf'.format(prefix))
      
        # # by hand GOF calculation
        os.system('combine -M FitDiagnostics -d datacard_ws_{PREFIX}_ratio.txt -n _fit_result_{PREFIX} --saveShapes --saveWithUncertainties --dataset data_obs'.format(PREFIX=prefix))
        fitDiag = rt.TFile.Open('fitDiagnostics_fit_result_{PREFIX}.root'.format(PREFIX=prefix),'r')

        byhand_gof = 0
        bw = 100
        f = rt.TCanvas('f','f',1000,800)
        f.cd()
        for cat in ['rej', 'acc']:
            bkgd = fitDiag.Get('shapes_fit_b/{cat}/background'.format(cat=cat))
            if not bkgd:
              continue
            bkgd.Scale(bw) # need to multiply by bin width for some reason?
            data = fitDiag.Get('shapes_fit_b/{cat}/data'.format(cat=cat))
            if cat=='rej':
                bkgd.Draw('hist')
                bkgd.SetMinimum(5000)
                bkgd.SetMaximum(15000)
                bkgd.GetXaxis().SetRangeUser(100,1700)
            else:
                bkgd.Draw("histsame")
            data.SetMarkerStyle(20)
            data.SetMarkerColor(rt.kBlack)
            for i in range(0,bkgd.GetNbinsX()):
                bw = bkgd.GetBinWidth(i)
                data.SetPointEXlow(i,0)
                data.SetPointEXhigh(i,0)
                data.SetPoint(i,data.GetX()[i], bw*data.GetY()[i]) # need to multiply by bin width for some reason?
                data.SetPointEYlow(i,bw*data.GetErrorYlow(i)) # need to multiply by bin width for some reason?
                data.SetPointEYhigh(i,bw*data.GetErrorYhigh(i)) # need to multiply by bin width for some reason?
            data.Draw('samepez')


            for i in range(0,bkgd.GetNbinsX()):
                x = bkgd.GetBinCenter(i+1)
                fi = bkgd.GetBinContent(i+1)
                di = data.GetY()[i]
                gofi = 2*(fi - di + di*rt.TMath.Log(di/fi)) # see eq. 14 of http://cousins.web.cern.ch/cousins/ongoodness6march2016.pdf

                # expect each bin to give GOF contribution ~ O(1)
                if gofi>5:
                    print('Bin center {}'.format(x))
                    print('{cat} bin {i}: fi={fi}, di={di}, gofi={gofi}'.format(cat=cat,i=i,fi=fi,di=di,gofi=gofi))
                    print(" -> BIG GOF CONTRIBUTION in this bin")
                byhand_gof += gofi
            f.SetLogy()
            f.Draw()
            fitDiag.Close()
        print("by hand obs GOF = {}".format(byhand_gof))
        print("combine obs GOF = {}".format(obs_gof))
        print("combine exp GOF = {}".format(np.mean(exp_gof)))

        print("Done with quantile ",qkey )
        print("p-value = {}".format(pval_toys))
        print("Observed GOF = {}".format(obs_gof))
