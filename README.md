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
