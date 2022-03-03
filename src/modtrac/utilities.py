#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:32:33 2022
@author: diesel
"""

import os
import pickle
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# define constants
J_eV = 6.242E18
eV_J = 1/J_eV
m_n = 1.675E-27  # kg/mol


def get_keys(folder):
    """Find run keys in moderator studies folder."""
    dirlist = [item for item in os.listdir(folder) if
               os.path.isdir(os.path.join(folder, item))]
    if 'scripts' in dirlist:
        dirlist.remove('scripts')
    return dirlist


def read_pickle_dict(filename, folder='C:/Users/Avram/Dropbox (MIT)/MIT'
                     '/Research/Thesis/Sections/mod_studies/'):
    """Read ptrac data dictionary from pickle file."""
    if folder == '':
        folder = Path.cwd()    # set folder ='' to use current directory
    path = Path(folder, filename).with_suffix('.pkl')
    try:
        print(f'Reading {path}...')
        with open(path, 'rb') as f:
            data = pickle.load(f)
    except IOError:
        print('Pickle file not found.')
        data = {}
    return data


def write_pickle(data=[], filename='ptracs', folder=''):
    """Store data as pickle file."""
    if folder == '':
        folder = Path.cwd()
    path = Path(folder, filename).with_suffix('.pkl')
    print(f'Writing variables to \'{str(path)}\'', end="...")
    with open(path, 'wb') as f:
        pickle.dump(data, f)
        print('Done!')


def t_to_d(t, E):
    """Convert TOF to distance."""
    return np.sqrt(2*E*eV_J/m_n)*t


def f_chisq(d, l, d_0):
    """Compute chi-squared function."""
    return np.exp(-(d-d_0)/l)*(((d-d_0)/l)**2)/(2*l)


def calc_mod_dist(df, e_lo=1.0, e_hi=100.0, e_res=1.0, d_offset=0.,
                  d_hi=200., d_res=2.5, plot=0):
    """Calculate moderator distance distribution."""
    if e_res < 0:
        e_lin = range(0, int(len(df)/np.abs(e_res))+1)
        means = np.zeros([len(e_lin), 2])
        sigmas = np.zeros([len(e_lin), 2])
    else:
        e_lin = np.linspace(e_lo, e_hi, int((e_hi - e_lo)/e_res) + 1)
        # initialize arrays to store mean and sigma with errors
        means = np.zeros([len(e_lin), 2])
        sigmas = np.zeros([len(e_lin), 2])
    # for each energy group, fit function to the equivalent
    # moderation distance distribution
    for i, E in enumerate(e_lin):
        if e_res < 0:
            if i != e_lin[-1]:
                data = df.nlargest((i+1)*np.abs(int(e_res)), 'E')[
                    np.abs(int(e_res))*i:np.abs(int(e_res))*(i+1)]
            else:
                data = df.nlargest((i+1)*np.abs(int(e_res)), 'E')[
                    np.abs(int(e_res))*i:]
            print(np.median(data.E))
        else:
            # subset data within energy
            data = df.loc[(df.E < E + e_res) & (df.E > E - e_res)]
        d_lo = df['d'].describe()['min'] - d_offset
        vals, bin_edges = np.histogram(data.d - d_offset,
                                       bins=int((d_hi-d_lo)/d_res),
                                       range=[d_lo, d_hi], weights=data.weight,
                                       density=True)
        bin_centers = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
        popt, pcov = curve_fit(f_chisq, bin_centers, vals, p0=[2, d_lo],
                               bounds=[[0, d_lo], [100, d_lo + 10]])
        # calculate stdev on optimal fit parameters
        perr = np.sqrt(np.diag(pcov))
        lamb = popt[0]
        mean = 3*lamb
        mean_err = 3*perr[0]
        means[i, 0] = mean
        means[i, 1] = mean_err
        sigma = np.sqrt(3*lamb**2)
        sigma_err = np.sqrt(3)*perr[0]
        sigmas[i, 0] = sigma
        sigmas[i, 1] = sigma_err
        print(f'Counts: {len(data)}; '
              f'E: {E:5.1f}; '
              f'λ: {mean/10:.2f} +/- {perr[0]/10:.2f} cm; '
              f'σ: {sigma/10:.2f} +/- {perr[1]/10:.2f} cm')
        # for low energy
        if plot == 1:
            plt.figure()
            plt.plot(bin_centers, vals, linestyle='none', marker='s',
                     color='black', label='MCNP')
            elin = np.linspace(d_lo, d_hi, 1+int((d_hi-d_lo)*10))
            plt.plot(elin, f_chisq(elin, *popt), color='blue', label='fit')
            plt.xlabel('EQUIVALENT DISTANCE (mm)', labelpad=10)
            plt.xlim([d_lo, d_hi])
            plt.ylabel('NORMALIZED COUNTS', labelpad=10)
            plt.ylim(bottom=0)
            plt.legend(loc='upper right')
    return means, sigmas


def get_ng_ratio(key, mod_dict):
    """Get ratio of epithermal (1-100 eV) neutrons to 2.2 MeV gammas."""
    vals = []
    vals_err = []
    for mode in ['n', 'p']:
        ptrac = mod_dict[mode][key]
        nps_DT = 5e8
        norm = (ptrac.nps/nps_DT)*ptrac.geometry['area_mod']
        if not ptrac.data_axial.empty:
            vals.append(sum(ptrac.data_axial['weight'])/norm)
            vals_err.append(np.sqrt(sum(ptrac.data_axial['weight']))/norm)
            # print(f'{val:1.2e}' '\t' f'{ptrac.label}')
    if len(vals) == 2:
        ratio = vals[0]/vals[1]
        ratio_err = ratio*np.sqrt((vals_err[0]/vals[0])**2 + (
            vals_err[1]/vals[1])**2)
        print(f'ratio: {ratio:1.2e} +/- {ratio_err:1.0e}' '\t' f'{key}')
        return ratio, ratio_err, vals
    else:
        return 0, 0, vals


def get_n(mod_dict, keys, thicks):
    """Get neutron efficiency."""
    effs = []
    errs = []
    for key, thick in zip(keys, thicks):
        n_axial = sum(mod_dict[key][thick].data_axial['weight'])
        nps = mod_dict[key][thick].nps
        eff = n_axial/nps
        err = np.sqrt(n_axial)/nps
        print(f'Efficiency: {eff:1.2e} +/- {err:1.2e} for {key}-{thick}')
        effs.append(eff)
        errs.append(err)
    return effs, errs
