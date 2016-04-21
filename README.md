# How to use

Go to user directory

`cd user`

## If you already have your final trees

Load the shared library and use the returnNLOweight function. 

The example code is available below.

C++:

`root -l example.C`

Python:

`python -i example.py`


The example code creates the shared library out of [functionmacro.C](https://github.com/CMS-HTT/NLOReweightingTool/blob/master/user/functionmacro.C) and load it at the beginning of the ROOT session. 

This way, ROOT can recognize the reweighting function.


**DO NOT MODIFY functionmacro.C unless you need your own defined function**


## If you want to produce Ntuples with the NLO weights (not recommended)

You can read, whenever you want: 

`file = TFile('/afs/cern.ch/user/y/ytakahas/public/forHTT/Reweight.root')`

`func = file.Get('your desired weight function for mA, tanb point')`

`func->Eval(higgs pT)`

This should give you the new weights and you can store the values.


# How to produce weight file (for administrative persons)

Go to producer directory

`cd producer`

`python shelf.py`

This will create Python shelve object, which stores tan(beta) vs Yukawa couplings from Stefan's file, locating in the producer/data directory. The resulting shelve object is already uploaded to the git (Yukawa_A.db) and not necessary to redo this step, unless we want to do the heavy Higgs (H).

`python ReweightingTool.py`

This macro will Produce reweighting file with the settings in config.ini (which should be easy to understand).

This macro compare the PY8 and MSSM distribution and derive the weight (MSSM/PY8 ratio) as a function of Higgs pT. The fitting function is exponential + pol3. (empirically, this was the best).

The macro also produces validation.pdf, where you can check if your fitting went OK or not.
You can change output ROOT file location as you wish.
