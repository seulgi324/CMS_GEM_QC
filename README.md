# Setup

This code works on lxplus9 (default lxplus) with installations of the packages below.
```
pip3 install pandas
pip3 install numpy==1.22.4
pip3 install mplhep
pip3 install fpdf
pip3 install openpyxl
```
It's based on the raw data format in 904, .xlsm for QC3 and .txt for QC4. 
Date is written in yyyy/mm/dd.


# How to run
1. Change /DATA/PATH/ to your path.
2. Change /FONT/PATH/ to your path (fonts are available in ./font/ dir in this repository).
3. Run the following command:
```
python3 QC34_all.py -mt [] -mn [] -d3 [] -d4 []
```
For GE2/1, 'mt' should be the type of module (ex.) M6, M7, ...)
For ME0, 'mt' should be 'ME0'
d3 and d4 are date of the QC3 and QC4 test respectively.
