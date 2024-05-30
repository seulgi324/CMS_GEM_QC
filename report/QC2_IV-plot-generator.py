# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:10:57 2021

@author: Felipe Ramirez
Works
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit
import numpy as np
import csv
import math
import statistics
import os
import fnmatch
import sys
"------------------USER INPUT"
data_folder='QC2_results/data_GE21_foils_20231001/' #Choose data folder
index_to_open=int(sys.argv[1])
 #Choose the index of the file to be opened
threshold=7 #Threshold (nA) to add a point to the I-V plot. This is used to avoid including current spikes in the I-V plot
savetofile=True #True if you want to save the plot contents in the txt file.
"----------------------------"


print('The following files were found: ')
file_names=[]
files_list=os.listdir(data_folder)
for file in files_list:
    if(fnmatch.fnmatch(file, 'QC2LONG_PART1_*.txt')):
        isIV_plot=file.find('IVplot')
        if(isIV_plot==-1): #This is done to exclude from the list the already generated I-V plots
            print(file)
            file_names.append(file)
    
with open(data_folder+'/'+file_names[index_to_open]) as csv_file:
    csv_content = csv.reader(csv_file, delimiter='\t')
    data_list = list(csv_content)
del data_list[0:6] #Delete the headers

voltage_list=[]
current_list=[]
time_list=[]
start=False
ramp_up=False
store=False
update_list=False

voltage_list_to_plot=[]
current_list_to_plot=[]
err_current_list_to_plot=[]

for i in range(len(data_list)):
    voltage_list.append(float(data_list[i][0]))
    current_list.append(float(data_list[i][1]))
    time_list.append(float(data_list[i][2]))
##This function extracts the I-V plot, considering that the raw I-V data folows the sequence:
##Ramp-up=True(I!=0) -> Ramp-up=False(I!=0) -> Scale change (I=0) -> Data for the I-Vplot (I!=0)-> Scale change (I=0) -> Ramp-up =True (I!=0) ->.... 

for j in range(800):
    if(voltage_list[j+2]-voltage_list[j]>2):
        #print("Ramping up")
        ramp_up=True
        start=False
        store=False
    else:
        ramp_up=False
    #if(start==False and ramp_up==False and current_list[j]!=0): #---->A
        #print("A")
        #print(voltage_list[j],current_list[j])
    if(start==False and ramp_up==False and current_list[j]==0 and voltage_list[j]!=0):#--->B
        #print("B")
        #print(voltage_list[j],current_list[j])
        start=True
        voltage_list_to_average=[]
        current_list_to_average=[]
    elif(start==True and current_list[j]!=0):
        #print("C")
        #print(voltage_list[j],current_list[j])
        store=True
        voltage_list_to_average.append(voltage_list[j])
        current_list_to_average.append(current_list[j]*1000.0) #Convert to nA
        update_list=True
    elif(start==True and store==True and current_list[j]==0):
        #print("D")
        if(len(voltage_list_to_average)>0 and update_list==True):
            update_list=False
            try:
                if(statistics.mean(current_list_to_average)<threshold):
                    voltage_list_to_plot.append(statistics.mean(voltage_list_to_average))
                    current_list_to_plot.append(statistics.mean(current_list_to_average))
                    err_current_list_to_plot.append(statistics.stdev(current_list_to_average)/math.sqrt(len(current_list_to_average)))
            except ValueError:
                voltage_list_to_plot.pop(-1)
                current_list_to_plot.pop(-1)
                #list_to_save=["Voltage (V)"+"\t"+"Current (nA)"+"\t""Error_current (nA)"]
                #np.savetxt(data_folder+'/'+file_names[index_to_open].replace('.txt','_IVplot.txt'),list_to_save,  delimiter='\t',fmt='%s')
            else:
                print('Imon= ',statistics.mean(current_list_to_average),' nA. The current was higher than the threshold. The point was not added to the I-V plot')
            voltage_list_to_average=[]
            current_list_to_average=[]

#Remove the points with very high current (sometimes there is a spike druring the I-V plot generation)
# voltage_list_to_plot.pop(12)
# current_list_to_plot.pop(12)
# err_current_list_to_plot.pop(12)
# voltage_list_to_plot.pop(8)
# current_list_to_plot.pop(8)
# err_current_list_to_plot.pop(8)
print('Imon values added to the I-V plot:')
print(len(voltage_list_to_plot))
print(voltage_list_to_plot)
print(current_list_to_plot)
print(err_current_list_to_plot)

list_idx = []
for i in range(0,len(voltage_list_to_plot)-1):
    if abs(voltage_list_to_plot[i] - voltage_list_to_plot[i+1]) < 5:
       list_idx.append(i+1)

for idx in sorted(list_idx, reverse = True):
    del voltage_list_to_plot[idx]
    del current_list_to_plot[idx]
    del err_current_list_to_plot[idx]

print(len(voltage_list_to_plot))
print(voltage_list_to_plot)
print(current_list_to_plot)
print(err_current_list_to_plot)

#Plot
plt.plot(time_list,voltage_list,'*')
plt.xlabel('Voltage (V)')
plt.ylabel('Current (uA)')
plt.close()
plt.plot(time_list,current_list,'*')
plt.xlabel('Voltage (V)')
plt.ylabel('Current (uA)')
plt.close()
plt.errorbar(x=voltage_list_to_plot,y=current_list_to_plot,yerr=err_current_list_to_plot,fmt='*')
plt.xlabel('Voltage (V)')
plt.ylabel('Current (nA)')
#plt.savefig(data_folder+file_names[index_to_open].replace('.txt','_IVplot.png'))
plt.close()

#Save data to file

list_to_save=[]
list_to_save.append(['Voltage (V)','Current (nA)','Error_current (nA)'])
for i in range(len(voltage_list_to_plot)):
    temp_list=[]
    temp_list.append(voltage_list_to_plot[i])
    temp_list.append(current_list_to_plot[i])
    temp_list.append(err_current_list_to_plot[i])       
    list_to_save.append(temp_list)
if(savetofile==False):  list_to_save=[['Voltage (V)','Current (nA)','Error_current (nA)']] #If the user choose not to save the list, only the headers will be saved.
np.savetxt(data_folder+'/'+file_names[index_to_open].replace('.txt','_IVplot.txt'),list_to_save, delimiter='\t', fmt='%s')

if (len(current_list_to_plot)==0):
    list_to_save=["Voltage (V)"+"\t"+"Current (nA)"+"\t""Error_current (nA)"]
    np.savetxt(data_folder+'/'+file_names[index_to_open].replace('.txt','_IVplot.txt'),list_to_save,  delimiter='\t',fmt='%s')


