# COPS analysis

Currently, this package requires pluq to function. Make sure pluq is installed by following the instructions here. 
https://github.com/kfritzsc/pluq

To use: 
a) as a python package: download files, cd to this directory, then run: 

pip install .
pip install -r requirements.txt 

b) as a poky plugin:
download files, make sure all requirements from requirements.txt are installed. 
One way to do this is to cd to this directory, then run: 
pip install -r requirements.txt

Then, 
    1) Start POKY
    2) Select spectrum
    3) use the 2-letter command "py", then click Load Module... and open cops_init.py
    4) initialize with spectrum names and experiment settings
    5) click on a peak, then click on "Calculate"


Spectrum requirements: 
NMRpipe- or UCSF-format. All must be aligned. 

pluq references: 

1. K. J. Fritzsching, Mei Hong, K. Schmidt-Rohr. "Conformationally Selective Multidimensional Chemical Shift Ranges in Proteins from a PACSY Database Purged Using Intrinsic Quality Criteria " J. Biomol. NMR 2016, 64 (2), 115–130 doi:10.1007/s10858-016-0013-5
2. Lee, W.; Yu, W.; Kim, S.; Chang, I.; Lee, W. PACSY, a Relational Database Management System for Protein Structure and Chemical Shift Analysis. J Biomol NMR 2012, 54 (2),169–179. doi:10.1109/BIBMW.2012.6470267
