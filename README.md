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

If there is problem in installing mplhep, please change matplotlib version to 3.7.5.


# How to run
1. Change /DATA/PATH/ to your path.
2. Change /FONT/PATH/ to your path (fonts are available in ./font/ dir in this repository).
3. Run the following command:
```
python3 QC34_all.py -mt [Module_type] -mn [Module_number] -d3 [Date_of_QC3_data] -d4 [Data_of_QC4_data]
```
For GE2/1, 'mt' should be the type of module (ex.) M6, M7, ...)
For ME0, 'mt' should be 'ME0'
d3 and d4 are date of the QC3 and QC4 test respectively.

For example, if you want to make a report for GE21-MODULE-M6-0080 with QC3 data produced on 2024 Oct 08 and QC4 data generated on 2024 Oct 10, run this command:
```
python3 QC34_all.py -mt M6 -mn 0080 -d3 20241008 -d4 20241010
```
If you want to make a report for ME0-MODULE-0053 with QC3 data produced on 2024 Mar 23 and QC4 data generated on 2024 Mar 24, run this command:
```
python3 QC34_all.py -mt ME0 -mn 0053 -d3 20240323 -d4 20240324
```
