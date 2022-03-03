#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:59:49 2022

@author: diesel
"""

import numpy as np
import matplotlib.pyplot as plt

from .ptracmod import PtracMod
from .utilities import get_keys, get_n, calc_mod_dist
from .plotter import r_plot_new


def main():
    """Script for processing ptrac moderator files."""
    # set data parameters
    folder = ('C:/Users/Avram//Dropbox (MIT)/MIT/research/collaborations/IVAC'
              '/github/IVAC/MCNP/data/')
    # 'C:/Users/Avram//Dropbox (MIT)/MIT-PNNL/MCNP/moderator-studies/'
    keys = get_keys(folder)
    keys = ['D2O', 'D2O', 'HDPE', 'HDPE', 'BHDPE', 'BHDPE', 'BHDPE',
            'H2O', 'H2O', 'frontmod']
    # ['cyl-10cm-10cmback', 'finite-10cm-8cmback', 'poly-b-10']
    thicks = ['100mm', '200mm', '100mm', '200mm', '100mm', '200mm', '250mm',
              '100mm', '200mm', '100mm']
    keys = ['BHDPE', 'BeBHDPE', 'BeBHDPE',
            'PbBHDPE', 'PbBHDPE', 'PbHDPE',
            'BeHDPE', 'PbHDPE', 'BeHDPE', 'PbBHDPE', 'BeBHDPE',
            'PbHDPE', 'PbHDPE']
    thicks = ['200mm', '100mm', '50mm', '100mm', '50mm', '25mm',
              '25mm', '50mm', '50mm', '25mm-circa', '25mm-circa',
              '10mmcore', '25mmcore']
    keys = ['HDPE', 'BHDPE', 'BHDPE', 'PbBHDPE', 'BeBHDPE',
            'PbHDPE', 'PbHDPE', 'PbHDPE', 'PbHDPE']
    thicks = ['100mm', '100mm', '200mm', '25mm-circa', '25mm-circa',
              '10mmcore', '25mmcore', '50mmcore']
    labels = [key + '-' + thick for key, thick in zip(keys, thicks)]#['D2O', 'LDPE', 'Borated LDPE', 'Lithiated LDPE']
    # keys = ['poly-b-10', 'poly-b-12', 'poly-b-12alt', 'poly-b-12alt2', 'poly-b-15', 'poly-b-15alt', 'poly-b-15alt2', 'poly-b-15alt3']
    # labels = ['10 cm', '12cm', '12cm alt', '15cm', '15cm alt', '15cm alt2', '15cm alt3']
    # # keys = ['B_HDPE_1E-1', 'B_HDPE_5E-1', 'boron-loaded', 'finite-10-10cmback', 'Li-LDPE', 'B-LDPE', 'B5-LDPE', 'D2O']
    # keys = ['cyl-10cm-10cmback', 'finite-10', 'finite-20', 'finite-10-center']
    mode = 'n'
    l_tof = 200 # in cm
    mod_dict = {}
    for key in keys:
        mod_dict[key] = {}
    for key, thick in zip(keys, thicks):
        ptrac = PtracMod(key=key, thick=thick, mode=mode, folder=folder)
        ptrac.get_MCNP_params()                     # get run paramters
        ptrac.read_ptrac(l_tof=l_tof, save='hdf')   # process run data
        mod_dict[key][thick] = ptrac                 # store in dictionary
    
    # ptrac = PtracMod(key='cyl-12_5cm', mode='n', folder=folder)
    # ptrac.get_MCNP_params()                     # get run paramters
    
    colors = ['black',  'blue', 'red', 'green', 'purple', 'orange', 'brown']
    
    # subset neutrons within axial radius at given TOF distance
    rads = np.linspace(1, 10, 10)
    areas = [np.pi*rad**2 for rad in rads]
    
    # define experimental parameters
    ID_mod = 5 # inner radius
    OD_mod = 15 # outer radius
    OD_shield = 35 # box diameter
    
    # read in data
    # flag_process = 0
    # if flag_process == 1:
    #     data = read_pickle_dict('ptracs')
    #     keys = list(data.keys())
    # else:
    #     data = {}
    #     keys = []
    # process new runs
    # # runs_unproc =  [run for run in runs if run not in keys]
    # data, keys = process_ptrac(run_params, folder, data, l_tof)
    # print(runs_unproc)
    # write_pickle(data)
    
    # keys = ['finite-10-10cmback', 'B_HDPE_1E-1', 'B_HDPE_5E-1', 'boron-loaded', 'Li-LDPE', 'B-LDPE', 'B5-LDPE', 'D2O']
    # labels = ['unborated HDPE', '0.1\% borated HDPE', '0.5\% borated HDPE', '1.0\% borated HDPE', '7.5\% lithiated HDPE', '1\% borated HDPE', '5\% borated HDPE', 'D2O']
    
    # keys = ['finite-10-10cmback', 'Li-LDPE', 'B-LDPE', 'B5-LDPE']
    # labels = ['pure LDPE', '7.5\% lithiated LDPE', '1\% borated LDPE', '5\% borated LDPE']
    
    # make plots of counts vs. radial distance
    # r_plot(mod_dict, data='data', keys=keys, labels=labels,
    #         modes=modes, w_bin=20, save=False, err=False, save_label='modified_HDPE')
    
    # r_plot(mod_dict, data='data_axial', keys=keys, labels=labels,
    #         modes=modes, w_bin=5, save=False, err=True, save_label='modified_HDPE')
    
    # get ratio of neutrons to gammas
    # ratios_ng = {key:get_ng_ratio(key, mod_dict) for key in mod_dict['n'] if get_ng_ratio(key, mod_dict) != (0, 0)}
    
    effs, errs = get_n(mod_dict, keys, thicks)
    
    # bvals, __ = np.histogram(mod_dict['BHDPE']['200mm'].data_axial['E'], range =[0, 100], bins=10,
    #                           weights=np.ones(len(mod_dict['BHDPE']['200mm'].data_axial['E']))/mod_dict['BHDPE']['200mm'].nps)
    # pvals, __ = np.histogram(mod_dict['HDPE']['200mm'].data_axial['E'], range =[0, 100], bins=10,
    #                           weights=np.ones(len(mod_dict['HDPE']['200mm'].data_axial['E']))/mod_dict['HDPE']['200mm'].nps)
    # bvals/pvals
    
    means, sigmas = calc_mod_dist(mod_dict['PbHDPE']['10mmcore'].data_axial,
                                  d_hi=200, d_res=0.5, e_res=100.0, plot=1)
    
    r_plot_new(mod_dict, data='data_axial', keys=keys, thicks=thicks, labels=labels,
               w_bin=5, save=False, err=True)
    
    plt.figure(figsize=(16, 9))
    for key, thick in zip(keys, thicks):
        means, sigmas = calc_mod_dist(mod_dict[key][thick].data_axial,
                                      d_hi=200, d_res=0.5, e_res=10.0, plot=0)
        plt.hist(mod_dict[key][thick].data_axial['E'], histtype='step', range =[0, 100], bins=100,
                 weights=np.ones([len(mod_dict[key][thick].data_axial)])/mod_dict[key][thick].nps,
                 lw=2, label=key + '-' + thick)
        print(key, thick, mod_dict[key][thick].nps)
    plt.xlim(1, 100)
    #plt.xscale('log')
    plt.xlabel('ENERGY [eV]')
    #plt.yscale('log')
    plt.ylabel('COUNTS')
    plt.legend()
    plt.tight_layout()
    
    # key_axial = 'finite-10-10cmback'
    # df_axial = mod_dict['n'][key_axial].data_axial
    
    """
    
    plt.figure(figsize=(16, 9))
    plt.errorbar([7.0, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 12.0],
                 [7.17e-6, 7.47e-6, 7.66e-6, 7.56e-6, 7.87e-6, 7.89e-6, 7.43e-6, 7.25e-6, 7.19e-6],
                 yerr=[1.89e-7, 2.73e-7, 1.24e-7, 1.94e-7, 1.98e-7, 1.26e-7, 1.93e-7, 2.69e-7, 2.68e-7],
                 marker='s', capsize=2, linestyle='None', elinewidth=2, color='black')
    plt.xlabel(r'AXIAL POSITION [cm]')
    plt.ylabel(r'$\phi_{epi}$')
    plt.tight_layout()
    
    # plt.figure(figsize=(16, 9))
    # plt.plot([2.5, 5.0, 7.5, 10.0, 12.5, 13.5, 15.0],
    #          [6.68e-6, 7.97e-6, 9.57e-6, 8.29e-6, 7.08e-6, 6.45e-6, 5.85e-6],
    #          marker='s', color='black', linestyle='None')
    # plt.tight_layout()
    
    # cumulative plot of radial position of axial neutrons at source-moderator assembly front face
    r_lin = np.linspace(0, 17.5, 1001)
    cumsums = {}
    plt.figure(figsize=(16, 9))
    for key in keys:
        df = df_axial[key].copy()
        r_mod = run_params[key]['r_mod']
        hist, bin_edges = np.histogram(df.r, bins=int(10*r_mod), weights=df.weight*1e6/run_params[key]['nps'], range=[0, r_mod]) # weights are arbitrary scaled
        cumsums[key] = np.cumsum(hist)
        f_int = interp1d(bin_edges[1:], cumsums[key], bounds_error=False)
        plt.plot(r_lin, f_int(r_lin),lw=2, label=run_params[key]['label'])
    plt.legend()
    plt.xlim(0, 10)
    plt.ylim(1, 20)
    plt.xlabel('RADIAL POSITION (CM)', labelpad=10)
    plt.ylabel('CUM. COUNTS (a.u.)', labelpad=10)
    plt.tight_layout()
    
    
    
    # plot moderation time vs. energy
    df = df_axial.copy()
    plt.figure(figsize=(16, 9))
    plt.scatter(df.E, df.t, lw=1)
    
    
    hist_plot(df_axial, mode='E', d_lo=100, d_hi=150, d_res=0.25) # for radial
    
    # plot time distribution of gamma
    def f_linear(x, m, b):
        return m*x + b
    
    def gamma_fit(run_params, df=df_axial, keys=[], plot=0, t_lo=0, t_hi=200, n_bins=200):
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
        if keys == []:
            keys = df.keys()
        if plot == 1:
            plt.figure(figsize=(16, 9))
        for key, color in zip(keys, colors):
            vals, bin_edges = np.histogram(df[key]['p'].t, bins=n_bins, range=[t_lo, t_hi],
                                           weights=df[key]['p'].weight/run_params[key]['data']['p']['nps']*1e9)
            bin_centers = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
            popt, pcov = curve_fit(f_linear, bin_centers, np.log(vals), p0 = [-1, np.max(vals)])
            perr = np.diag(pcov)
            t_lin = np.linspace(t_lo, t_hi, int((t_hi-t_lo)/0.1)+1)
            print(f'Fit: y =  {np.exp(popt[1]):1.2e}*exp({popt[0]:1.2e}*x)')
            print(f'Tau = ({-1/popt[0]:3.0f} +/- {np.sqrt(perr[0])/popt[0]**2:1.0f}) us')
            if plot == 1:
                plt.errorbar(x=bin_centers, y=vals, yerr=np.sqrt(vals),
                             drawstyle='steps-mid', capsize=2, elinewidth=1, lw=2, color=color,
                             label=run_params[key]['data']['p']['label'])
                plt.plot(t_lin, np.exp(f_linear(t_lin, *popt)), color=color, label=key + ' - fit')
    
        if plot == 1:
            plt.xlim(t_lo, t_hi)
            plt.ylim(bottom=1)
            plt.yscale('log')
            plt.xlabel('TIME (us)', labelpad=10)
            plt.ylabel('COUNTS')
            plt.legend()
            plt.tight_layout()
            plt.gca().set_yticks([20, 50, 100, 200, 300])
            plt.gca().get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    
    gamma_fit(run_params, plot=1, t_lo=5, t_hi=200, n_bins=390)
    
    
    df = mod_dict['n']['cyl-10cm-10cmback'].data
    
    df_close = [df.rad2 <= 10**2]
    
    
    
    
    
    # make plot of on-axis epithermal neutrons per radius at TOF distance
    plt.figure(figsize=(25, 12))
    for key, color, label, n in zip(keys, colors, labels, nps):
        df = data[key].copy()
        vals = []
        for rad, area in zip(rads, areas):
            df_close = df.loc[df.rad2 <= rad**2]
            vals.append(sum(df_close.weight)*(1e8/nps[n]/area))
        plt.plot(rads, vals, marker='s', markersize=6, color=color, label=label)
    plt.xlim([1-0.1, max(rads)+1])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(f'RADIAL DISTANCE (cm) FROM TOF AXIS @ {l_tof/100:.2f} m', labelpad=30)
    plt.ylabel('ON-AXIS EPITHERMAL NEUTRONS / CM^2 / 1E8 N', labelpad=10)
    plt.legend()
    plt.tight_layout()
    
    rad_close = 5 # in cm
    
    
    
    
    
    
    
    
    
    # set parameters
    d_lo = 100
    d_hi = 200
    d_res = 2. # in mm
    
    t_lo = 0.0
    t_hi = 10.0
    t_res = 0.05 # in us
    
    e_lo = 1.0
    e_hi = 100.0
    e_res = 1.0 # in eV
    e_lin = np.linspace(e_lo, e_hi, int((e_hi - e_lo)/e_res) + 1)
    
    # filter to neutrons exiting HDPE moderator
    df_mod = df.loc[(df.x**2 + df.y**2 > ID_mod**2) & (df.x**2 + df.y**2 < OD_mod**2)]
    df_mod_axial =df_mod.loc[df_mod.rad <= rad_close]
    
    polar_flag = 0
    if polar_flag == 1:
        # calculate polar histogram
        rbins = np.linspace(0, df['r'].max(), 7)
        abins = np.linspace(-np.pi, np.pi, 25)
        hist, _, _ = np.histogram2d(df_axial.phi, df_axial.r, bins=(abins, rbins))
        A, R = np.meshgrid(abins, rbins)
    
        # generate polar plot
        fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
        pc = ax.pcolormesh(A, R, hist.T, cmap="magma_r")
        fig.colorbar(pc)
    
    # plot x,y position of neutrons exiting front face of source assembly
    [mod_plot(mod_dict, key, OD_shield=OD_shield, ID_mod=ID_mod) for key in keys]
    # mod_plot(df_mod, OD_shield=OD_shield, ID_mod=ID_mod, OD_mod=OD_mod)
    mod_plot(df_axial, OD_shield=OD_shield, ID_mod=ID_mod, OD_mod=OD_mod, bins=[20, 20])
    # rad_plot(df_axial, weighted=1)
    
    # plot moderation time or TOF distance as a function of output neutron energy
    hist_plot(mod_dict['n']['D2O'].data_axial, mode=mode, d_lo=100, d_hi=200, d_res=0.25) # for radial
    # hist_plot(df_axial, mode=mode, d_hi=150, d_res=1.0) # for axial
    
    means, sigmas = calc_mod_dist(df_axial, d_lo=0, d_hi=200, d_res=0.5, e_res=1.0, plot=1)
    mod_dist_plot(means, mode, unit='cm')
    mod_dist_plot(sigmas, mode, unit='cm')
    
    # for each 1 eV neutron group, plot histogram of total collisions
    IQR_coll = []
    avg_coll = []
    counts_coll = []
    for E in e_lin:
        data = df.loc[(df.E < E + e_res) & (df.E > E - e_res)]
        result = np.percentile(data.coll.values, q=[0, 25, 50, 75, 100])
        counts_coll.append(len(data.coll.values))
        avg_coll.append(np.mean(data.coll.values)) # multiply by 100 to convert m to cm
        IQR_coll.append((result[3] - result[1]))
    
    
    sigmas_coll =[i/(2*np.sqrt(c)) for i, c in zip(IQR_coll, counts_coll)]
    weights_raw = [1/s**2 for s in sigmas_coll]
    
    plt.figure(figsize=(16, 9))
    plt.errorbar(x=e_lin+e_res/2, y=avg_coll,
                  yerr=sigmas_coll,
                  linestyle='none', elinewidth=1, capsize=2, marker='s', color='black')
    plt.xlabel('ENERGY (eV)')
    plt.ylabel('COLLISIONS')
    plt.xlim([0, 100])
    plt.tight_layout()
    plt.legend()
    
    # make a 3D plot of moderation time or equivalent TOF distance distribution
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x_vec = np.arange(t_lo, t_hi+t_res, t_res) # time (us)
    z_vec = np.arange(e_lo, e_hi, e_res) # energy (eV)
    verts = []
    IQR = []
    avg = []
    counts = []
    for z in z_vec:
        data = df[(df.E >= z-e_res/2) & (df.E < z+e_res/2)]
        # hist bin values for mod time within e cut
        if mode == 'd':
            vals, bin_edges = np.histogram(data.d.values, range=[d_lo, d_hi],
                                            bins=int((d_hi-d_lo)/d_res))
            bin_edges = [b-100 for b in bin_edges]
            # with open(str(int(z)).zfill(2)+'_eV.dat', 'w+') as f:
            #     np.savetxt(f, data.d.values, fmt='%f', delimiter='\n')
            y_vec, __ = np.histogram(data.d.values, range=[d_lo, d_hi], bins=int((d_hi-d_lo)/d_res)+1)
            result = np.percentile(data.d.values, q=[0, 25, 50, 75, 100])
            ax.set_xlabel(r'DISTANCE (mm)', labelpad=40)
            print(f'avg for E = {z:.1f} eV: {result[2]-d_lo:.3f} mm')
            print(f'IQR for E = {z:.1f} eV: {result[3] - result[1]} mm'+'\n')
        elif mode == 't':
            y_vec, __ = np.histogram(data.t.values, range=[t_lo, t_hi], bins=int((t_hi-t_lo)/t_res)+1)
            result = np.percentile(data.t.values, q=[0, 25, 50, 75, 100])
            ax.set_xlabel(r'TIME ($\mu$s)', labelpad=40)
            print(f'avg for E = {z:.1f} eV: {result[2]:.3f} us')
            print(f'IQR/2 for E = {z:.1f} eV: {(result[3] - result[1])/2} us'+'\n')
        verts.append(list(zip(x_vec, y_vec)))
        # calculate and print means and IQRs for selected distribution
        counts.append(len(data.t.values))
        avg.append(result[2])
        IQR.append((result[3] - result[1]))
    poly = PolyCollection(verts,facecolor='white')
    poly.set_edgecolor('black')
    ax.add_collection3d(poly, zs=z_vec, zdir='y')
    ax.set_xlim3d(0, 3)
    ax.set_ylabel('ENERGY (eV)', labelpad=40)
    ax.set_ylim3d(e_lo, e_hi)
    ax.set_zlabel('COUNTS', labelpad=40)
    ax.set_zlim3d(10, 10000)
    ax.set_zscale('log')
    
    sigmas =[i/(2*np.sqrt(c)) for i, c in zip(IQR, counts)]
    weights_raw = [1/s**2 for s in sigmas]
    weights = [w/sum(weights_raw) for w in weights_raw]
    
    # print distribution of average
    plt.figure(figsize=(16, 9))
    plt.errorbar(x=e_lin[:-1]+e_res/2, y=avg,
                  yerr=sigmas,
                  linestyle='none', elinewidth=1, capsize=2, marker='s', color='black')
    plt.xlabel('FINAL NEUTRON ENERGY (eV)')
    if mode == 'd':
        plt.ylabel('TOF DISTANCE (mm)')
    elif mode == 't':
        plt.ylabel('MODERATION TIME (us)')
    plt.xlim([e_lo, e_hi])
    plt.ylim([min(avg)/1.1, 1.1*max(avg)])
    plt.tight_layout()
    plt.xticks(np.arange(0, 110, step=10))
    
    plt.figure(figsize=(16, 9))
    plt.plot(e_lin[:-1]+e_res/2, [i/2 for i in IQR], marker='s', linestyle='None', color='black')
    plt.xlim([e_lo, e_hi])
    plt.ylim([0, 20])
    plt.xlabel('FINAL NEUTRON ENERGY (eV)')
    plt.ylabel('DELTA TOF DISTANCE (mm)')
    
    val = sum([i*w for i, w in zip(avg, weights)])/sum(weights)
    sigma = np.sqrt(sum([s**2 * w**2 for s, w in zip(sigmas, weights)]))
    
    plt.text(20, 15, r'$\bar{x}$ = ' + f'{val-d_lo:.2f} +/- {sigma:.2f} cm')
    
    # xlin = np.linspace(1, 100, 1000)
    # # Calculate the slope and y-intercept of the trendline
    # fit = np.polyfit(1/(E_lin**0.75), IQR, 1)
    # yfit = fit[0]/xlin**0.75 + fit[1]
    
    # plt.figure()
    # plt.scatter(e_lin[:-1]+e_res/2, IQR, s=60, color='black')
    # # plt.plot(xlin, yfit,'black')
    # plt.xlim([min(e_lin)-0.5, max(e_lin)+0.5])
    # plt.xlabel ('Energy (eV)')
    # if mode == 'd':
    #     plt.ylabel(r'IQR, $\Delta$d (cm)')
    # elif mode =='t':
    #     plt.ylabel(r'IQR, $\Delta$t ($\mu$s)')
    """
    
    # sigma mod plot for IVAC review
    e_lin = [1.27, 2.12, 3.53, 5.72, 9.27, 14.7, 23.1, 35.7, 54.7, 81.6]
    plt.figure(figsize=(16, 9))
    plt.errorbar(x=e_lin, y=sigmas[:, 0]/10, yerr=sigmas[:, 1]/10,
                 linestyle='none', elinewidth=1, capsize=2, marker='s', ms=10, color='black')
    plt.hlines(np.mean(sigmas[:,0])/10, 1, 95, linestyle='dashed', lw=2, color='red', label='fit')
    plt.xlabel('FINAL NEUTRON ENERGY (eV)', labelpad=10)
    plt.ylabel(r'$\sigma_{d}$ (cm)', labelpad=10)
    plt.xlim([0.9, 100])
    plt.xscale('log')
    plt.ylim([1, 1.8])
    plt.tight_layout()
    plt.legend()
    # plt.xticks(np.arange(0, 110, step=10))
    
    # radial mod vs. axial mod plot for IVAC review
    plt.figure(figsize=(16, 9))
    labels = ['radial', 'axial']
    colors = ['blue', 'red']
    for key, thick, label, color in zip(keys, thicks, labels, colors):
        x = mod_dict[key][thick].data_axial['E']
        nps = mod_dict[key][thick].nps
        hist, bin_edges = np.histogram(x, range=[1, 100], bins=99)
        plt.errorbar(bin_edges[:-1], hist/nps, np.sqrt(hist)/nps,
                     capsize=2,
                     drawstyle='steps-mid', lw=2, color=color, label=label)
    plt.xlim(1, 100)
    plt.xscale('log')
    plt.xlabel(r'Energy (eV)')
    plt.yscale('log')
    plt.ylabel(r'Counts (Norm.)')
    plt.legend()
    plt.tight_layout()


if __name__ == '__main__':
    main()
