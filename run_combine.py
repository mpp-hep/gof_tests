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
import math
from subprocess import Popen
import shlex 
import cmsstyle

def makeWS(data_rej,data_acc,signal_rej,signal_acc, outname,quantile,qr_unc_file=None):
  
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
  
  # Name efficiency per quantile to allow for combination
  w.factory('eff_%s[%f]'%(quantile,efficiency))
  w.var('eff_%s'%quantile).setConstant(True)

  
  empty_hist = rt.TH1D('empty_hist','empty_hist', data_rej.GetNbinsX(), data_rej.GetXaxis().GetXbins().GetArray())

  for iBinX in range(1,data_rej.GetNbinsX()+1):
      empty_hist.SetBinContent(iBinX,1)
      rej_bin = data_rej.GetBinContent(iBinX)
      
     # If adding a per-bin uncertianty on the CR or the efficiency estimate
      # w.factory('eff_%s_%i[%f,%f,%f]'%(quantile,iBinX,efficiency,efficiency*0.98,efficiency*1.02,))
      # w.var('eff_%s_%i'%(quantile,iBinX)).setConstant(True)
      
      # WIP that was abandoned: Adding Uncertainty on yield in reject based on QR envelope: qr_rms/qr_cutValue. Use function option qr_unc_file=None which needs a json file as input
#       rms_qr = qr_unc_file['q90'][iBinX-1][2]/qr_unc_file['q90'][iBinX-1][1]
#       print('crBin%i_In[%.1f,%.1f,%.1f]'%(iBinX,rej_bin,rej_bin*(1-rms_qr),rej_bin*(1+rms_qr)))
#       w.factory('crBin%i_In[%.1f,%.1f,%.1f]'%(iBinX,rej_bin,rej_bin*(1-rms_qr),rej_bin*(1+rms_qr)))

      
      w.factory('crBin%i_In[%.1f]'%(iBinX,rej_bin))
      # Fix control region bin count
      w.var('crBin%i_In'%iBinX).setConstant(True)
      
      w.factory('crBin_q%s_%i[0,-10,10]'%(quantile,iBinX))
      w.var('crBin_q%s_%i'%(quantile,iBinX)).setConstant(False)
      if rej_bin !=  0.:
        power = 1/rt.TMath.Sqrt(rej_bin)
      else:
        power = -1.0

        w.var('crBin_q%s_%i'%(quantile,iBinX)).setConstant(True)
        
      #what is fit is actually (mjj+mjj/sqrt(mjj))^x_bin, meaning x_bin is representative of how many sigma acc is from rej
      w.factory("expr::crBin%iFunc_q%s('max(0,@0*pow(1.0+%f,@1))',crBin%i_In,crBin_q%s_%i)"%(iBinX,quantile,power,iBinX,quantile,iBinX))
      w.factory("expr::bin%iFunc_q%s('max(0,@0*@1)',eff_%s,crBin%iFunc_q%s)"%(iBinX,quantile,quantile,iBinX,quantile))
      rej_bin_functions.add(w.function('crBin%iFunc_q%s'%(iBinX,quantile)))
      acc_bin_functions.add(w.function('bin%iFunc_q%s'%(iBinX,quantile)))
    
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
  '''
  datacard_ratio += 'eff_%s   flatParam\n'%quantile
  # for i in range(1,n_bins+1): #If you have bin-ny-bin uncertianties on the QR efficiency
  #     datacard_ratio += 'eff_%s_%i   flatParam\n'%(quantile,i) 
  for i in range(1,n_bins+1):
      datacard_ratio += 'crBin_q%s_%i   flatParam\n'%(quantile,i)
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
  
  # Which crossection you are running the GOF for
  xsecs   = [0,10,20,40,60,80,100]
  xsecs   = [0]
  
  # Which quantiles to run GOF for
  #columns = {0.5: 21, 0.1: 18, 0.3: 22, 0.9: 17, 0.7: 20, 0.01: 19}
  # columns = {0.01: 19, 0.1: 18, 0.3: 22,0.5: 21, 0.7: 20}
  columns = {0.01: 19, 0.1: 18, 0.3: 22,0.5: 21}
  
  # RUn combine on all quantiles
  doAllQuantiles = True
  
  # You have already run the above, and now want to combine the quantiles
  combineAll     = False
  
  # Informative name for the Combine outputs
  tag = "sig_GtoWW35naReco"  
  
  NTOYS = 5 # Multiplied by 10, running in a i=10 loop with different seed each time  
  if doAllQuantiles:
  
    for xs in xsecs:
      print("\nFor xsec={}:".format(xs) )
      for key in sorted(columns.keys()):
          print("For q = {}".format(key) )
          quantile = int(key*100)
          # Get histograms
          f = rt.TFile.Open('histograms.root',"r")
          data_acc = f.Get('data_acc_{}fb_q{}p'.format(xs,quantile)); data_acc.SetDirectory(0);
          data_rej = f.Get('data_rej_{}fb_q100p'.format(xs)); data_rej.SetDirectory(0);

          signal_acc     = f.Get('signal_acc_{}fb_q{}p'.format(xs,quantile)); signal_acc.SetDirectory(0);
          signal_rej     = f.Get('signal_rej_{}fb_q100p'.format(xs)); signal_rej.SetDirectory(0);
          f.Close()

          print('XS={}: SIG PASS = {} SIG FAIL = {} BKG PASS = {} BKG FAIL = {}'.format(xs,signal_acc.Integral(),signal_rej.Integral(), data_acc.Integral(),data_rej.Integral()))
          xaxis = data_acc.GetXaxis().GetXbins()
          min_bin = xaxis[0]
          max_bin = xaxis[len(xaxis)-1]
          n_bins  = len(xaxis)-1

          prefix = "q{}_xs{}_{}".format(quantile,xs,tag)
          makeWS(data_rej,data_acc,signal_rej,signal_acc,'datacard_ws_{PREFIX}.root'.format(PREFIX=prefix),quantile)

    for xs in xsecs:
        print("\nFor xsec={}:".format(xs) )
        for key in sorted(columns.keys()):
            print("For q = {}".format(key) )
            quantile = int(key*100)
            prefix = "q{}_xs{}_{}".format(quantile,xs,tag)
            os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d datacard_ws_{PREFIX}_ratio.txt -n Ratio_gof_{PREFIX} --dataset data_obs -v 2'.format(PREFIX=prefix))

            for i in range(10):
              os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d datacard_ws_{PREFIX}_ratio.txt -t {NTOYS} --toysFreq -n Ratio_gof_toys_{PREFIX}  --dataset data_obs -s {S} -v 0'.format(PREFIX=prefix,NTOYS=NTOYS,S=40+i))
            os.system('hadd -f higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.ALLTOYS.root higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.4*.root'.format(PREFIX=prefix))

    for xs in xsecs:
        print("\nxsec={}:\n".format(xs) )
        for key in sorted(columns.keys()):
            print("q = {}".format(key) )
            quantile = int(key*100)
            prefix = "q{}_xs{}_{}".format(quantile,xs,tag)
            # open file
            obs_gof_file = uproot.open('higgsCombineRatio_gof_{PREFIX}.GoodnessOfFit.mH120.root'.format(PREFIX=prefix))
            obs_gof = obs_gof_file['limit'].arrays('limit')['limit'][0]

            exp_gof_file = uproot.open('higgsCombineRatio_gof_toys_{PREFIX}.GoodnessOfFit.mH120.ALLTOYS.root'.format(PREFIX=prefix))
            exp_gof = exp_gof_file['limit'].arrays('limit')['limit']
            
            print("Obs.   {:.1f}".format(obs_gof))
            print("Exp.   {:.1f}\n".format(np.mean(exp_gof)))    
          
  if combineAll:
      for xs in xsecs:
        print("\nFor xsec={}:".format(xs) )
        os.system('combineCards.py q1=datacard_ws_q1_xs{XS}_ratio.txt q10=datacard_ws_q10_xs{XS}_ratio.txt q30=datacard_ws_q30_xs{XS}_ratio.txt q50=datacard_ws_q50_xs{XS}_ratio.txt q70=datacard_ws_q70_xs{XS}_ratio.txt q90=datacard_ws_q90_xs{XS}_ratio.txt &> combined_{XS}.txt'.format(XS=0))
        os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d combined_{XS}.txt -n Ratio_combined_{XS} --dataset data_obs -v 2'.format(XS=xs))  
        os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d combined_{XS}.txt  -t 5 --toysFreq -n Ratio_combined_{XS} --dataset data_obs -v 2 -s 42'.format(XS=0))    