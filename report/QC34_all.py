import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.font_manager as font_manager
import mplhep as hep
from fpdf import FPDF
import argparse
import os
'''
python3 QC34_all.py -mt M2 -mn 0003 -d3 20230908 -d4 20230908  # GE21
python3 QC34_all.py -mt ME0 -mn 0003 -d3 20230908 -d4 20230908  # ME0
'''

plt.style.use(hep.style.CMS)

# Fit function for QC3
def func(x, m, t):
    return m*np.exp(-t*x)

def label(mt, mn):
    if 'ME0' in mt: return 'ME0-MODULE-{}'.format(mn)
    else: return 'GE21-MODULE-{}-{}'.format(mt, mn)

def qc3_data(d3, label):    
    data_qc3 = pd.read_excel(r'/DATA/PATH/QC3_{}_{}.xlsm'.format(label, d3), engine='openpyxl')
    return data_qc3

def qc4_data(d4, label):    
    data_qc4 = pd.read_csv('/DATA/PATH/QC4_{}_{}.txt'.format(label, d4), sep="\t", skiprows=[0, 1, 2, 3, 4, 5, 7])
    return data_qc4

def qc3_plot(mn, d3, label):
    fig, ax = plt.subplots()
    fig.set_figheight(9)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab")
    data = qc3_data(d3, label)  
    time = np.array(data['Seconds'].tolist())
    time_hr = time/3600
    pressure = np.around(np.array(data['Pressure (mBar)'].tolist()), 2)
    temperature = data['Temperature (C)'].tolist()[0]
    atm = data['Atm Pressure (mBar)'].tolist()[0]
    t0 = int(np.where(time == 1.00)[0])
    t1 = int(np.where(time == 3600.00)[0])
    p0 = (26, 0.001)
    popt, pcov = curve_fit(func, time_hr[t0+3:t1-3], pressure[t0+3:t1-3])
    a, b = popt
    ymin = pressure[t1]-2.0
    ymax = pressure[t0]+2.0
    plt.ylim([ymin, ymax])
    plt.plot(time_hr[1:t1+1], pressure[1:t1+1], 'ok', markersize=2, label=label)
    plt.plot(time_hr[1:t1+1], func(time_hr[1:t1+1], *popt), 'b-', linewidth=1.5, label=r'$P(t) = [p0]e^{-t/\tau}$')
    t = 1/b
    plt.text(0.6, ymax-(ymax-ymin)/4-(ymax-ymin)/20, 'p0              %.2f mbar' % a, fontsize = 20, color='b')
    plt.text(0.6, ymax-(ymax-ymin)/4-(ymax-ymin)*30/(20*14), r'$\tau$                %.2f h' % t, fontsize = 20, color='b', fontweight='bold')
    if 'ME0' in label:
        plt.text(0, ymin+(ymax-ymin)/60+(ymax-ymin)*30/(20*13), 'ME0 Module Production', fontsize=20) 
    else:
        plt.text(0, ymin+(ymax-ymin)/60+(ymax-ymin)*30/(20*13), 'GE2/1 Module Production', fontsize=20) 
    plt.text(0, ymin+(ymax-ymin)/60+(ymax-ymin)/20, 'Gas = $CO_{2}$', fontsize=20) 
    lg = font_manager.FontProperties(#weight='bold',
                                     style='normal', size=20)
    plt.xlabel('Time [h]')
    plt.ylabel('Pressure [mbar]')
    legend = plt.legend(loc='upper right', prop=lg)
    plt.savefig('./plot/QC3_{}_{}.png'.format(label, d3), dpi=50)
    plt.clf()
    return b

def qc4_plot(mn, d4, label):
    fig, ax = plt.subplots()
    fig.set_figheight(9)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab")
    dt = qc4_data(d4, label) 
    voltage = (np.array(dt['Vmon'])/1000).tolist()
    current = np.array(dt['Imon']).tolist()
    coeff = np.polyfit(current[:-2], voltage[:-2], 1)
    r_m = 1000*coeff[0]
    poly1d_fn = np.poly1d(coeff)
    plt.plot(current, voltage, 'ok', label=label)
    plt.plot(current, poly1d_fn(current), '-r', label='Fit')
    if 'ME0' in label:
        plt.text(50, 4.7, 'ME0 Module Production', fontsize=20) #, fontproperties=font)
    else:
        plt.text(50, 4.7, 'GE2/1 Module Production', fontsize=20) #, fontproperties=font)
    plt.text(50, 4.33, 'Gas = $CO_{2}$', fontsize=20) #, fontproperties=font)
    plt.text(50, 3.96, '$R_{n}$ = 5.0 $M\Omega$' % r_m, fontsize=20) #, fontproperties=font)
    plt.text(50, 3.59, '$R_{m}$ = %.3f $M\Omega$' % r_m, fontsize=20) #, fontproperties=font)
    lg = font_manager.FontProperties(#weight='bold',
                                     style='normal', size=20)
    plt.xticks([200, 400, 600, 800, 1000])
    plt.xlabel('Divider Current $I_{divider} \ [\mu$A]')
    plt.ylabel('Applied Voltage V [kV]')
    plt.legend(loc='lower right', prop=lg)
    plt.savefig('./plot/QC4_{}_{}.png'.format(label, d4), dpi=50)
    return r_m

def qc34_report(mn, d3, d4, b, r_m, label):
    data = qc3_data(d3, label)  
    pressure = np.around(np.array(data['Pressure (mBar)'].tolist()), 2)
    temperature = data['Temperature (C)'].tolist()[0]
    atm = data['Atm Pressure (mBar)'].tolist()[0]
    time = np.array(data['Seconds'].tolist())
    #print(time)
    time_hr = time/3600
    t0 = int(np.where(time == 1.00)[0])
    t1 = int(np.where(time == 3600.00)[0])
    t = 1/b
    dt = qc4_data(d4, label) 
    # Generate PDF file
    pdf = FPDF()
    pdf.add_page()
    pdf.set_auto_page_break(True, margin = 1.0)
    pdf.add_font('FreeSans', '', '/FONT/PATH/FreeSans.ttf', uni=True)
    pdf.add_font('FreeSansB', '', '/FONT/PATH/FreeSansBold.ttf', uni=True)
    
    omega = str('\u03A9')
    print(omega)
    m = '\u2098'
    
    # Header
    pdf.set_font('FreeSansB', '', 22)
    pdf.cell(0, 10, 'QC3 & QC4 Report on {}'.format(label),ln=1, align='C')
    pdf.set_font('FreeSansB', '', 15)
    pdf.image('./plot/QC3_{}_{}.png'.format(label, d3), x=0, y=43,  w=100)
    pdf.ln(10)
    pdf.cell(100, 20, 'QC3 Result', ln=1)
    pdf.set_font('FreeSansB', '', 10)
    pdf.ln(12)
    pdf.cell(85) 
    pdf.cell(70, 12, 'Test Date')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '{}-{}-{}'.format(d3[:4], d3[4:6], d3[6:8]), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(85)
    pdf.cell(70, 12, 'Temperature  ({}C) / Pressure (mBar)'.format(chr(176)))
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '{} / {}'.format(temperature, atm), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(85)
    pdf.cell(70, 12, 'Pressure drop (mBar/hr)')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '%.2f' % (pressure[t0]-pressure[t1]), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(85)
    pdf.cell(70, 12, 'Time constant (hr)')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '%.2f' % t, ln=1, align='R')
    
    
    pdf.ln(7)
    pdf.set_font('FreeSansB', '', 15)
    pdf.cell(100, 20, '', ln=1)
    pdf.image('./plot/QC4_{}_{}.png'.format(label, d4), x=0, y=150, w=100)
    pdf.cell(100, 20, 'QC4 Result', ln=1)
    pdf.set_font('FreeSansB', '', 10)
    pdf.ln(17)
    pdf.cell(85)
    pdf.cell(70, 12, 'Test Date')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '{}-{}-{}'.format(d4[:4], d4[4:6], d4[6:8]), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(85)
    
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(45.8, 12, 'The measured resistance R')
    pdf.set_font('FreeSansB', '', 5)
    pdf.cell(3, 13.7, 'm')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(4.5, 12, '(M')
    pdf.cell(16.7, 12, '{})'.format(omega))
    
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '%.3f' % r_m, ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(85)
    pdf.cell(70, 12, 'Resistance Deviation (%)')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 12, '%.2f' % (100*(5.0-r_m)/5.0), ln=1, align='R')
    
    pdf.output('./pdf/QC34_report_{}.pdf'.format(label))
    


if __name__=="__main__":
    plt.style.use(hep.style.CMS)
    parser = argparse.ArgumentParser()
    parser.add_argument("-mt", "--module_type", dest="module_type", help="module type")
    parser.add_argument("-mn", "--module_number", dest="module_number", help="module number")
    parser.add_argument("-d3", "--qc3_date", dest="qc3_date", help="qc3 test date (YYYYMMDD)")
    parser.add_argument("-d4", "--qc4_date", dest="qc4_date", help="qc4 test date (YYYYMMDD)")
    args = parser.parse_args()
    print(args.module_number, args.qc3_date, args.qc4_date)

    os.makedirs('./plot', exist_ok=True)
    os.makedirs('./pdf', exist_ok=True)
    label = label(args.module_type, args.module_number)
    b = qc3_plot(args.module_number, args.qc3_date, label)
    r_m = qc4_plot(args.module_number, args.qc4_date, label)
    qc34_report(args.module_number, args.qc3_date, args.qc4_date, b, r_m, label)

