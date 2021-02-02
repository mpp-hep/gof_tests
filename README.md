# Goodness-of-fit tests for model-independent searches

This code is based on J. Duartes script https://github.com/mpp-hep/JetAnomaly_CHEP19/Simple_combine.ipynb. It uses Higgs combine tool via CMSSW:

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.13
scramv1 b clean; scramv1 b
```

Create the Conda environment

```
conda env create -f environment.yml
conda activate jet-anomaly
cmsenv
```

Every time you open a new shell
```
conda activate jet-anomaly
cd CMSSW_10_2_13/src/
cmsenv
cd ../../
```

The script contains two almost identical scripts: One Jupyter notebook for illustrative purposes and quick checks, and one python script that can be imported and used in another script

To run the notebook from lxplusABC
```
[your-user-name@lxplusABC gof_tests]$ jupyter notebook --no-browser
[I 20:37:09.457 NotebookApp] The Jupyter Notebook is running at:
[I 20:37:09.457 NotebookApp] http://localhost:XXXX/?token=some-token
```

Note the port number XXXX and on your laptop do
```
 ssh -N -f -L YYYY:localhost:XXXX your-user-name@lxplusABC.cern.ch
 where YYYY is some local port of your choosing, eg 8888
 ```
 Then open the link
  ```
http://localhost:XXXX/?token=some-token
  ```
and you should be in the notebook
  
Open gof.ipynb and work through the cells