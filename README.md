# NLOReweightingTool

Production of the Reweighting file (only administrative person should do):

> cd producer

> python shelf.py (will create Python shelve object, which stores tan(beta) vs Yukawa couplings from Stefan's file, locating in the producer/data directory) 
-> The resulting shelve object is already uploaded to the git (Yukawa_A.db) and not necessary to redo this step, unless we want to do the heavy Higgs (H) case

> python ReweightingTool.py
-> Produce reweighting file with the settings in config.ini (which should be easy to understand)
-> Compare the PY8 and MSSM distribution and derive the weight (MSSM/PY8 ratio) as a function of Higgs pT
-> Fit the ratio with 6th order polinominal
-> This will create Reweight.root in the user directory



Plotting macros for the user:

> cd user

The idea is to load the functionmacro.C macro, which provide the weight file, as a function of mass, tan(beta) and Higgs pT (in each event)

This will show you the example in python:
> python example.py

This will show you the example in C++:
> root -l example.C

