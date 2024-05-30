# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 17:21:24 2021

@author: luife
"""
import os
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
from fpdf import FPDF
from datetime import datetime
from fpdf.enums import XPos, YPos

"Edit this part of the acode according to the location of your QC2 files"
"---------------------------------------------------------------------------"
data_folder='QC2_results/data_GE21_foils_20230910'
subfolder='/'
part1_file=''
part1_file='QC2LONG_PART1_GE21-FOIL-M3-G12-KR-B02-0018_20230910_13-06'
megger_file='QC2FAST_GE21-FOIL-M3-G12-KR-B02-0018_20230910'
all_channels_file=''
all_channels_file='QC2_all_channels_monitor_20230911_10-25'
include_IV_plot=True   #Set to False if you do not want to include the I-V plot in the pdf report and the /plots folder.
"----------------------------------------------------------------------------"


#Create folder for plots and pdf reports if it doesn't exit:
if not os.path.exists(data_folder+'/plots'):
    os.makedirs(data_folder+'/plots')
#Create folder for reports if it doesn't exit
if not os.path.exists(data_folder+'/pdf_reports'):
    os.makedirs(data_folder+'/pdf_reports')

############################Open all the associated files
with open(data_folder+'/megger'+'/'+megger_file+'.txt') as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    megger_list = list(csv_content)


with open(data_folder+subfolder+part1_file+'.txt') as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    data_list_part1 = list(csv_content)
description_list=data_list_part1[0:5]
CH_number=int(description_list[0][1][2]) #The channel number is only one digit in position [0][1] after the substring "CH"
del data_list_part1[0:6] #Delete the headers
voltage_list_part1=[]
current_list_part1=[]
time_list_part1=[]

with open(data_folder+subfolder+part1_file+'_IVplot.txt') as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    data_list_IV = list(csv_content)
del data_list_IV[0]
IV_voltage=[]
IV_current=[]
IV_current_error=[]

with open(data_folder+subfolder+all_channels_file+'.txt') as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    data_list_allCH = list(csv_content)
del data_list_allCH[0:2]
voltage_list_part2=[]
current_list_part2=[]
time_list_part2=[]

lenpart1=len(part1_file)
with open(data_folder+subfolder+'QC2NOTES_'+part1_file[14:lenpart1-15]+'.txt') as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    notes_list = list(csv_content)
del notes_list[0]

###################################Generate the plots
for i in range(len(data_list_part1)):
    voltage_list_part1.append(float(data_list_part1[i][0]))
    current_list_part1.append(float(data_list_part1[i][1]))
    time_list_part1.append(float(data_list_part1[i][2]))
for i in range(len(data_list_IV)):
    IV_voltage.append(float(data_list_IV[i][0]))
    IV_current.append(float(data_list_IV[i][1]))
    IV_current_error.append(float(data_list_IV[i][2]))
for i in range(len(data_list_allCH)):
    voltage_list_part2.append(float(data_list_allCH[i][CH_number+1]))
    current_list_part2.append(float(data_list_allCH[i][CH_number+9]))
    time_list_part2.append(float(data_list_allCH[i][0])/3600.0)

fig, axc1 = plt.subplots()
fig.set_figheight(9)
fig.set_figwidth(10)
hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=24)
axc1.set_ylabel('Current [nA]', fontsize=24, color='blue', loc='top')
axc1.set_xlabel('Time [s]', fontsize=24, loc='right')
axc1.tick_params(axis="x", direction='in', labelsize=20, length=8)
axc1.tick_params(axis="y", direction='in', labelsize=20, colors='blue', length=8)
axc1.plot(time_list_part1,current_list_part1,'s-', color='blue')

axv1= axc1.twinx()
axv1.set_ylabel('Voltage [V]', fontsize=24, color='red', loc='top')
axv1.tick_params(axis="y", direction='in', labelsize=20, length=8, colors='red')
axv1.plot(time_list_part1,voltage_list_part1,'s-', color='r')
plt.savefig(data_folder+'/plots/'+part1_file+'-VI-t.png',bbox_inches='tight')
plt.close()

fig2, ax = plt.subplots()
fig2.set_figheight(9)
fig2.set_figwidth(10)
hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=24)
ax.errorbar(x=IV_voltage,y=IV_current,yerr=IV_current_error,fmt='s', color='k')
ax.set_yscale("log")
ax.set_yticks([1, 7, 10])
ax.tick_params(axis="x", direction='in', labelsize=20, length=8)
ax.tick_params(axis="y", direction='in', labelsize=20, length=8)
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel('Voltage [V]', fontsize=24, loc='right')
ax.set_ylabel('Current [nA]', fontsize=24, loc='top')

ax.axhline(y=7, color='r', linestyle='-')
if(include_IV_plot):plt.savefig(data_folder+'/plots/'+part1_file+'-I-V.png',bbox_inches='tight')
plt.close()

fig3, axc2 = plt.subplots()
fig3.set_figheight(9)
fig3.set_figwidth(10)
hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=24)
axc2.set_ylabel('Current [nA]', fontsize=24, color='blue', loc='top')
axc2.set_xlabel('Time [hr]', fontsize=24, loc='right')
axc2.tick_params(axis="x", direction='in', labelsize=20, length=8)
axc2.tick_params(axis="y", direction='in', labelsize=20, colors='blue', length=8)
axc2.plot(time_list_part2,current_list_part2,'s-', color='blue')

axv2= axc2.twinx()
axv2.set_ylabel('Voltage [V]', fontsize=24, color='r', loc='top')
axv2.tick_params(axis="y", direction='in', labelsize=20, length=8, colors='red')
axv2.plot(time_list_part2,voltage_list_part2,'s-', color='r')
plt.savefig(data_folder+'/plots/'+part1_file+'-VI-t-long.png',bbox_inches='tight')


###################################Generate PDF report
pdf = FPDF()
pdf.add_page()
pdf.set_auto_page_break(True, margin = 1.0)
#### Generate header
pdf.set_font('helvetica', 'B', 22)
pdf.cell(40, 10, 'GE21 QC2 Report',new_x=XPos.LMARGIN, new_y=YPos.NEXT)
pdf.set_font('helvetica', 'B', 12)
for k in range (len(description_list)):
    pdf.cell(40,5, description_list[k][0]+' '+description_list[k][1],new_x=XPos.LMARGIN, new_y=YPos.NEXT)
pdf.set_font('helvetica', 'B', 16)
pdf.cell(300, 20, 'Acceptance test I (Megger test)')

#### Generate Megger table
pdf.set_font('helvetica', '', 10)
pdf.cell(300, 15, '',new_x=XPos.LMARGIN, new_y=YPos.NEXT)
line_height = pdf.font_size * 1.5
col_width = pdf.epw / 3  # distribute content evenly

for row in megger_list:
    for entry in row:
        pdf.multi_cell(col_width, line_height, entry, border=1, new_x=XPos.RIGHT, new_y=YPos.TOP, max_line_height=pdf.font_size)
    pdf.ln(line_height)
    
#### Place V-time, I-time and I-V plots of first part
pdf.set_font('helvetica', 'B', 16)
pdf.cell(300, 20, 'Acceptance test II -Part 1-')
pdf.set_font('helvetica', '', 10)
pdf.cell(300, 15, '',new_x=XPos.LMARGIN, new_y=YPos.NEXT)
pdf.image(data_folder+'/plots/'+part1_file+'-VI-t.png',pdf.get_x(),pdf.get_y(),pdf.epw*0.4)
#plt.xticks(fontsize=16); plt.yticks(fontsize=16);
if(include_IV_plot):pdf.image(data_folder+'/plots/'+part1_file+'-I-V.png',pdf.get_x()+pdf.epw*0.45,pdf.get_y(),pdf.epw*0.365)
pdf.ln(62)
pdf.set_font('helvetica', 'B', 16)
pdf.cell(300, 20, 'Acceptance test II -Part 2-')
pdf.set_font('helvetica', '', 10)
pdf.cell(300, 15, '', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
pdf.image(data_folder+'/plots/'+part1_file+'-VI-t-long.png',pdf.get_x(),pdf.get_y(),pdf.epw*0.4)
pdf.ln(62)
pdf.set_font('helvetica', 'B', 16)
pdf.cell(300, 20, 'Notes')
pdf.set_font('helvetica', '', 10)
pdf.cell(300, 15, '', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
for row in notes_list:
     pdf.cell(300,pdf.font_size,row[0],new_x=XPos.LMARGIN, new_y=YPos.NEXT)
pdf.output(data_folder+'/pdf_reports/'+'QC2REPORT_'+part1_file[14:44]+'_'+datetime.now().strftime("%Y%m%d_%H-%M")+'.pdf')

################################### Generate individual .txt files for the QC2 long Part 2
QC2_part2_list_to_save=[]
for i in range (len(description_list)-1):
    temp_header=[]
    temp_header.append(description_list[i][0])
    temp_header.append(description_list[i][1])
    temp_header.append('\t')
    QC2_part2_list_to_save.append(temp_header)
time_stamp_part2=all_channels_file[25:39]
QC2_part2_list_to_save.append(['Time_stamp:',time_stamp_part2,'\t'])
QC2_part2_list_to_save.append(['Voltage (V)','Current (uA)','Time (s)'])
for i in range(len(data_list_allCH)):
    temp_list=[]
    temp_list.append(data_list_allCH[i][CH_number+1])
    temp_list.append(data_list_allCH[i][CH_number+9])
    temp_list.append(data_list_allCH[i][0])
    QC2_part2_list_to_save.append(temp_list)
np.savetxt(data_folder+'/QC2LONG_PART2'+part1_file[13:44]+all_channels_file[24:39]+'.txt',QC2_part2_list_to_save, delimiter='\t', fmt='%s')
