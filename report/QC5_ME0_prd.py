import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.font_manager as font_manager
import mplhep as hep
from fpdf import FPDF
import argparse, os, ROOT
import warnings
import math
from copy import copy
from decimal import Decimal
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

'''
python3 QC5_ME0_prd.py -mn 0003 -d5 20230908 -pr 1700 -lr 1703 -dc 600
'''

plt.style.use(hep.style.CMS)

def func(x, a, b):
    return a*np.exp(b*np.array(x))

def gauss(x, a, b, c):
    return a*np.exp(-0.5*pow(((x-b)/c), 2))

def correction_factor(P, T):
    p0 = 964.4
    t0 = 297.1
    factor0 = (p0/P) * (T/t0)
    return pow(factor0, 0.43)

def gain(rate, current):
    primary_electron = 346
    e = 1.6e-19
    return current / (primary_electron * e * rate)

def qc5_eff_vals(mn, d51):
    data_dir = '/yourpath/'
    dr = pd.read_csv(data_dir+'QC5_ME0-MODULE-{}_{}.txt'.format(mn, d51), sep="\t", skiprows=[1, 2])
    dc = pd.read_csv(data_dir+'QC5_ME0-MODULE-{}_{}_currents_OFF_ON.txt'.format(mn, d51), sep="\t", header=None)
    pressure = dr.iloc[0, 3]
    temp = dr.iloc[0, 4] + 273.15
    imon = dr.iloc[:, 1].tolist()
    count_off = dr.iloc[:, 5].tolist()
    count_on = dr.iloc[:, 7].tolist()
    rates, ro_currents = [], []
    for i in range(len(count_on)):   # Rate
        rate = (count_on[i] - count_off[i]) / 10
        rates.append(rate)
    for im_ in range(len(imon)):   # Readout current
        ro_current = np.mean(np.array(dc.iloc[:, im_])) - np.mean(np.array(dc.iloc[:, im_+len(imon)]))
        ro_currents.append(abs(ro_current))
    im = np.where(np.array(imon) == 720)[0][0]
    r = rates[im]
    gains = [gain(r, ro_currents[i]) for i in range(len(imon))]
    return pressure, temp, imon, rates, ro_currents, gains

def qc5_eff_plot(mn, qc5_vals):
    pressure, temp, imon, rates, ro_currents, gains = qc5_vals

    # plot1: Rate & RO current
    fig, ax1 = plt.subplots()
    fig.set_figheight(8)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=24)
    warnings.filterwarnings('ignore')
    eff_ro = ax1.plot(imon, ro_currents, 'ok', markersize=7, color='green', label='Readout Current')
    fig.canvas.draw()
    offset_text = ax1.yaxis.get_offset_text().get_text()
    exp = int(offset_text.split('-')[1]) 
    ax1.set_xlabel('$I_{divider} \ [\mu$A]', fontsize=24, labelpad=7)
    ax1.set_ylabel(r'Readout Current [A] $\times$ $10^{{{}}}$'.format(exp), fontsize=24, labelpad=9)
    ax1.tick_params(axis='y', which='major', labelsize=20)
    ax1.tick_params(axis='x', which='major', labelsize=20)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Rate [Hz]', fontsize = 24, labelpad=10)
    ax2.tick_params(axis='y', which='major', labelsize=20)
    rate = ax2.plot(imon, rates, '^k', markersize=7, color='darkviolet', label='Rate')
    leg = rate + eff_ro
    labs = [l.get_label() for l in leg]
    ax1.legend(leg, labs, frameon=False, loc='lower right', fontsize=18)
    ax1.set_ylim(ymax=max(ro_currents)*1.14)
    ax2.set_ylim(ymax=max(rates)*1.1)
    plt.text(min(imon), max(rates)*1.02, 'ME0-MODULE-{}'.format(mn), fontsize=18)
    plt.text(min(imon), max(rates)*0.96, 'Gas = $Ar/CO_{2}$ (70/30)', fontsize=18)
    plt.text(min(imon), max(rates)*0.9, 'X-ray target: Ag', fontsize=18)
    ax1.grid()
    ax1.yaxis.get_offset_text().set_visible(False)
    plt.tight_layout()
    plt.savefig('./plot/QC5Part1_ME0-MODULE-{}_rate_rocurrent.png'.format(mn), dpi=150)
    plt.clf()  

    # plot2: Eff Gain before/after T/P correction
    fig, ax1 = plt.subplots()
    fig.set_figheight(8.5)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab")
    print("Correction Factor : ", correction_factor(pressure, temp))
    imon_corr = [correction_factor(pressure, temp) * i for i in imon]
    warnings.filterwarnings('ignore')
    popt1, pcov1 = curve_fit(func, imon[3:], gains[3:], p0=(1e-08, 3e-02))
    popt2, pcov2 = curve_fit(func, imon_corr[3:], gains[3:], p0=(1e-08, 3e-02))
    c20000_1 = math.log(20000.0/popt1[0]) / popt1[1]
    g700_1 = func(700, popt1[0], popt1[1])
    cg700_1 = math.log(g700_1/popt1[0]) / popt1[1]
    c20000_2 = math.log(20000.0/popt2[0]) / popt2[1]
    g700_2 = func(700, popt2[0], popt2[1])
    cg700_2 = math.log(g700_2/popt2[0]) / popt2[1]
    print("before TP: gain: {}, {}, {}".format(c20000_1, g700_1, cg700_1))
    print("before TP: g0 : {}, g1 : {}".format(popt1[0], popt1[1]))
    print("after TP: gain: {}, {}, {}".format(c20000_2, g700_2, cg700_2))
    print("after TP: g0 : {}, g1 : {}".format(popt2[0], popt2[1]))
    plt.yscale("log")
    plt.xlabel('$I_{divider} \ [\mu$A]')
    plt.ylabel('Effective Gas Gain')
    eff_gain = plt.plot(imon, gains, 'ok', color='blue', label='Effective Gas Gain before T/P correction')
    eff_gain_fit = plt.plot(imon, func(imon, *popt1), 'b-', linewidth=1.5, label=r'$G(I_{div}) = [g0]e^{[g1]I_{div}}$')
    eff_gain_corr = plt.plot(imon_corr, gains, '^k', color='red', label='Effective Gas Gain after T/P correction')
    eff_gain_corr_fit = plt.plot(imon_corr, func(imon_corr, *popt2), 'r-', linewidth=1.5, label=r'$G(I_{div}^{corr}) = [g0]e^{[g1]I_{div}^{corr}}$')
    leg = eff_gain + eff_gain_fit + eff_gain_corr + eff_gain_corr_fit
    labs = [l.get_label() for l in leg]
    plt.legend(leg, labs, frameon=False, loc='upper left', fontsize=18)
    plt.ylim(ymin=300, ymax=5*max(gains))
    ymin, ymax = plt.ylim()
    xx = 0.9333*max(imon)
    plt.text(xx, 1.7*ymin*pow(1.54, 4), 'ME0-MODULE-{}'.format(mn), fontsize=18)
    plt.text(xx, 1.7*ymin*pow(1.54, 3), 'Gas = $Ar/CO_{2}$ (70/30)', fontsize=18)
    plt.text(xx, 1.7*ymin*pow(1.54, 2), 'X-ray target: Ag', fontsize=18)
    plt.text(xx, 1.7*ymin*1.54, 'Gain at $I_{{divider}} 700 \mu$A = {}'.format(round(g700_1)), fontsize=18)
    plt.text(xx, 1.7*ymin, 'Gain at $I_{{divider}}^{{corr}} 700 \mu$A = {}'.format(round(g700_2)), fontsize=18)
    plt.grid()
    plt.grid(which='minor', alpha=0.3)
    plt.tight_layout()
    plt.savefig('./plot/QC5Part1_ME0-MODULE-{}_effective_gain.png'.format(mn), dpi=150)
    plt.clf()

    return popt1, round(g700_1), popt2, round(g700_2)

def qc5_uniformity_plot(mn, pr, lr, dc):
    fig, ax = plt.subplots()
    fig.set_figheight(9)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab")
    inFile = ROOT.TFile('/yourpath/ME0_MODULE_{}_{}uA_{}to{}_MMDAQ.root'.format(mn, dc, pr, lr))
 
    # plot1
    tree1 = inFile.Get("Summary/mgraph_ME0-MODULE-{}_ResponseFitPkPos_AllEta".format(mn))
    graphs = tree1.GetListOfGraphs()
    colors = ['red', 'green', 'blue', 'magenta', 'olive', 'orange', 'dimgrey', 'brown']
    for i in range(len(graphs)):
        a, b, c, d = np.array(graphs[i].GetX()), np.array(graphs[i].GetY()), np.array(graphs[i].GetEX()), np.array(graphs[i].GetEY())
        index_ = np.where(b > 2750)
        a, b, c, d = np.delete(a, index_), np.delete(b, index_), np.delete(c, index_), np.delete(d, index_)
        index_ = np.where(b < 1000)
        a, b, c, d = np.delete(a, index_), np.delete(b, index_), np.delete(c, index_), np.delete(d, index_)
        plt.errorbar(a, b, xerr=c, yerr=d, capsize=2, color=colors[i], ls='none', label='i$\eta$={}'.format(i+1))
    plt.text(0, 900, 'ME0-MODULE-{}'.format(mn), fontsize=16, ha='center') #, ha='right', va='top',)
    plt.text(0, 710, 'Gas = Ar/CO$_{2}$ (70/30)', fontsize=16, ha='center') #, ha='right', va='top',)
    plt.xlabel('Cluster Position [mm]', loc='right', fontsize=18)
    plt.ylabel('MPV [ADC count]', loc='top', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks([500, 1000, 1500, 2000, 2500, 3000], fontsize=18)
    plt.ylim(ymin=0, ymax=3000)
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0, 4, 1, 5, 2, 6, 3, 7]
    plt.legend([handles[i] for i in order], [labels[i] for i in order], loc='lower center', ncol=4, fontsize=20)
    plt.savefig('./plot/QC5Part2_ME0-MODULE-{}_1.png'.format(mn), dpi=150)
    plt.clf()

    # plot2
    fig, ax = plt.subplots()
    fig.set_figheight(7)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=22)
    y_bar = []
    for i in range(len(graphs)):
        y_bar.append(graphs[i].GetY())
    y_bar = np.array(y_bar)
    y_bar[y_bar < 1000] = -1
    y_bar = y_bar / np.mean(y_bar)
    my_cmap = copy(plt.cm.get_cmap('viridis'))
    my_cmap.set_under('white')
    plt.pcolormesh(y_bar, cmap=my_cmap, vmin=0.0000001)
    plt.text(46, 7.5, 'ME0 Module Production', fontsize=16, ha='right', va='top')
    plt.text(46, 6.6, 'ME0-MODULE-{}'.format(mn), fontsize=16, ha='right') #, ha='right', va='bottom')
    plt.xticks([8, 24, 40], [1, 2, 3], fontsize=17)
    plt.yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5], [1, 2, 3, 4, 5, 6, 7, 8], fontsize=17)
    plt.xlabel("$i\phi$", fontsize=17)
    plt.ylabel("$i\eta$", fontsize=17)
    plt.minorticks_off()
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    cb.set_label(label='Normalized peak pos. of cluster charge fits (ADC counts)', size=14.8, labelpad=15)
    plt.tight_layout()
    plt.savefig('./plot/QC5Part2_ME0-MODULE-{}_2.png'.format(mn), dpi=150)
    plt.clf()

    # plot3
    fig, ax = plt.subplots()
    fig.set_figheight(9)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab")
    canvas = inFile.Get("Summary/canv_ME0-MODULE-{}_ResponseFitPkPosDataset_AllEta;1".format(mn))
    tree3 = canvas.GetPrimitive("h_Summary_ResponseFitPkPosDataset")
    tree4 = canvas.GetPrimitive("fit_Summary_ResponseFitPkPosDataset")
    fit = tree3.Fit(tree4, "test")
    x, y, y_err = [], [], []
    for i in range(tree3.GetNbinsX()):
        y.append(tree3.GetBinContent(i))
        y_err.append(tree3.GetBinError(i))
        x.append(tree3.GetBinCenter(i))
    v = fit.Chi2() / fit.Ndf()
    v1 = fit.Parameter(0)
    v2 = fit.Parameter(1)
    v2_err = fit.ParError(1)
    v3 = fit.Parameter(2)
    v3_err = fit.ParError(2)
    v32 = v3/v2
    v32_err = np.sqrt( pow(v3_err / v2, 2) + pow((v3 * v2_err) / pow(v2, 2), 2))
    xerr = [(x[i+1] - x[i])/2 for i in range(len(x)-1)]
    xerr = xerr + [xerr[0]]
    x = np.array(x)
    y = np.array(y)
    x_ = np.linspace(min(x), max(x), 100)
    y_ = gauss(x_, v1, v2, v3)
    plt.errorbar(x, y, xerr=xerr, yerr=y_err, ls='none', color='blue')
    plt.plot(x_, y_, color='red')
    plt.text(min(x), 0.9*plt.ylim()[1], '$\chi^{2}$ / NDF = %.4f' % v, fontsize=16)
    plt.text(max(x), 0.9*plt.ylim()[1], 'ME0-MODULE-{}'.format(mn), fontsize=16, ha='right') #, ha='right', va='top',)
    plt.text(min(x), 0.83*plt.ylim()[1], '$\mu = %.2f \pm %.2f $' % (v2, v2_err), fontsize=16)
    plt.text(max(x), 0.83*plt.ylim()[1], 'Gas = Ar/CO$_2$ (70/30)', fontsize=16, ha='right') #, ha='right', va='top',)
    plt.text(min(x), 0.76*plt.ylim()[1], '$\sigma = %.2f \pm %.2f $' % (v3, v3_err), fontsize=16)
    plt.text(min(x), 0.69*plt.ylim()[1], '$\sigma/\mu = %.4f \pm %.5f $' % (v32, v32_err), fontsize=16)
    plt.ylim(ymin=0)
    plt.xlabel("Cluster MPV [ADC count]", fontsize=18)
    plt.ylabel("Number of Fits", fontsize=18)
    plt.savefig('./plot/QC5Part2_ME0-MODULE-{}_3.png'.format(mn), dpi=150)
    plt.clf()
     
    # plot4
    fig, ax = plt.subplots()
    fig.set_figheight(7)
    fig.set_figwidth(10)
    hep.cms.label(llabel="Preliminary", rlabel="CERN 904 Lab", fontsize=22)
    y_bar = []
    occ = []
    for i in range(1, 9):
        h = inFile.Get("SectorEta{0}/h_iEta{0}_hitPos;1".format(i))
        n_bins = h.GetNbinsX()
        a = np.array([h.GetBinContent(b) for b in range(1, n_bins + 1)])
        occ.append(np.array([h.GetBinContent(b) for b in range(1, n_bins + 1)]))
    y_bar = np.array(occ)
    nonzero = y_bar[y_bar > 0]
    vmean = np.mean(nonzero)
    my_cmap = copy(plt.cm.get_cmap('viridis'))
    my_cmap.set_under('white')
    plt.pcolormesh(y_bar, cmap=my_cmap, vmin=0.0000001, vmax=1.5*vmean)
    plt.text(370, 7.5, 'ME0 Module Production', fontsize=16, ha='right', va='top', color='white')
    plt.text(370, 6.6, 'ME0-MODULE-{}'.format(mn), fontsize=16, ha='right', color='white') #, ha='right', va='bottom')
    plt.yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5], [1, 2, 3, 4, 5, 6, 7, 8], fontsize=17)
    plt.xlabel("$strip$", fontsize=17)
    plt.ylabel("$i\eta$", fontsize=17)
    plt.minorticks_off()
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    cb.set_label(label='Hit Occupancy', size=14.8, labelpad=15)
    plt.tight_layout()
    print("check")
    plt.savefig('./plot/QC5Part2_ME0-MODULE-{}_4.png'.format(mn), dpi=150)
    plt.clf()

    return v32

def qc5_report(mn, dc, qc5_vals, g, ru):
    pressure, temp, imon, rates, ro_currents, gains = qc5_vals
    g_wo, g700_wo, g_w, g700_w = g

    # Generate PDF file
    pdf = FPDF()
    pdf.add_page()
    pdf.set_auto_page_break(True, margin = 1.0)
    pdf.add_font('FreeSans', '', '/yourpath/FreeSans.ttf', uni=True)
    pdf.add_font('FreeSansB', '', '/yourpath/FreeSansBold.ttf', uni=True)
    pdf.add_font('FreeSerif', '', '/yourpath/FreeSerif.ttf', uni=True)
    
    omega = str('\u03A9')
    print(omega)
    m = '\u2098'
    
    # Header
    pdf.set_font('FreeSansB', '', 22)
    pdf.cell(0, 10, 'QC5 Report on ME0-MODULE-{}'.format(mn),ln=1, align='C')
    pdf.set_font('FreeSansB', '', 15)
    pdf.image('./plot/QC5Part1_ME0-MODULE-{}_rate_rocurrent.png'.format(mn), x=20, y=47,  h=60)
    pdf.image('./plot/QC5Part1_ME0-MODULE-{}_effective_gain.png'.format(mn), x=110, y=47,  h=60)
    pdf.ln(8)
    pdf.cell(100, 7, 'Result of QC5 Part1 - Effective Gain Measurement', ln=1)
    pdf.set_font('FreeSans', '', 10)
    pdf.ln(2)
    pdf.cell(100, 12, 'Before T/P Correction (left)')
    pdf.cell(100, 12, 'After T/P Correction (right)', ln=1)
    pdf.ln(57)
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'Temperature  ({}C) / Pressure (mBar)'.format(chr(176)))
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{:.1f} / {}'.format(temp-273.15, pressure), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'g0 (before/after)')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{:.3e} / {:.3e}'.format(Decimal(g_wo[0]), Decimal(g_w[0])), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'g1 (before/after)')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{:.5f} / {:.5f}'.format(g_wo[1], g_w[1]), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'Gain at 700 {}A (before/after)'.format(chr(181)))
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{:.0f} / {:.0f}'.format(g700_wo, g700_w), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'Tested Sector')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, 'C20', ln=1, align='R')
    pdf.ln(8)
    pdf.set_font('FreeSansB', '', 15)
    pdf.cell(100, 7, 'Result of QC5 Part2 - Gain Uniformity Measurement', ln=1)
    pdf.image('./plot/QC5Part2_ME0-MODULE-{}_1.png'.format(mn), x=20, y=175,  h=60)
    pdf.image('./plot/QC5Part2_ME0-MODULE-{}_3.png'.format(mn), x=110, y=175,  h=61)
    pdf.ln(67)
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'Operating Divider Current'.format(chr(181)))
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{} {}A'.format(dc, chr(181)), ln=1, align='R')
    pdf.set_font('FreeSansB', '', 10)
    pdf.cell(150, 10, 'Response Uniformity')
    pdf.set_font('FreeSans', '', 10)
    pdf.cell(30, 10, '{:.2f} %'.format(100*ru), ln=1, align='R')
    pdf.set_font('FreeSans', '', 10)
    pdf.add_page()
    pdf.ln(80)
    pdf.image('./plot/QC5Part2_ME0-MODULE-{}_2.png'.format(mn), x=20, y=30,  h=60)
    pdf.image('./plot/QC5Part2_ME0-MODULE-{}_4.png'.format(mn), x=110, y=30,  h=61)
    note = "The left plot shows the normalized most probable value of cluster ADC aross the detector obtained from fitting. Empty bins mostly correspond to failed fits. A few slices show clear difference without consistent patterns on the same plot and these are mainly due to incorrect fit. The right plot shows a 2D hit occupancy plot used to check for missing or broken protection resistors on foil segments. ME0-MODULE-{} shows no protection resistor issues. Strips with abnormally high counts are caused by noise.".format(mn)
    pdf.multi_cell(0, 6, note, align='L')
    
    pdf.output('./pdf/QC5_report_ME0-MODULE-{}.pdf'.format(mn))


if __name__=="__main__":
    plt.style.use(hep.style.CMS)
    parser = argparse.ArgumentParser()
    parser.add_argument("-mn", "--module_number", dest="module_number", help="module number")
    parser.add_argument("-d5", "--qc5_date", dest="qc5_date", help="qc5 test date (YYYYMMDD)")
    parser.add_argument("-pr", "--pedestal_run", dest="pedestal_run", help="pedestal run number")
    parser.add_argument("-lr", "--last_run", dest="last_run", help="last run number")
    parser.add_argument("-dc", "--divider_current", dest="divider_current", help="operating divider current")
    args = parser.parse_args()

    os.makedirs('./plot', exist_ok=True)
    os.makedirs('./pdf', exist_ok=True)
    a = qc5_eff_vals(args.module_number, args.qc5_date)
    effgain = qc5_eff_plot(args.module_number, a)
    ru = qc5_uniformity_plot(args.module_number, args.pedestal_run, args.last_run, args.divider_current)
    qc5_report(args.module_number, args.divider_current, a, effgain, ru)
