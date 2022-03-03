#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:32:42 2022

@author: diesel
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as patches
import numpy as np


def set_plot_params():
    """Set plotting parameters."""
    plt.set_cmap('plasma')
    plt.close('all')
    plt.rcParams.update({'font.size': 24,
                         'font.family': 'serif',
                         'font.serif': ['Computer Modern Roman'],
                         'font.weight': 'bold'})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'


def r_plot_new(mod_dict, data='data', keys=[], thicks=[], mode='n',
               labels=[], colors=[], w_bin=0, err=False,
               save=False, save_label=''):
    """Plot neutron radial position at exit of source-moderator assembly."""
    r_DT = 5
    plt.figure(figsize=(16, 9))
    if keys == []:
        keys = mod_dict.keys()
    if thicks == []:
        thicks = [mod_dict[key].keys() for key in keys]
    if labels == []:
        labels = ['-'.join(key, thick) for key, thick in zip(keys, thicks)]
    if colors == []:
        colors = ['black', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
    for key, thick, label, color in zip(keys, thicks, labels, colors):
        ptrac = mod_dict[key][thick]
        df = getattr(ptrac, data)
        # see if r_mod in geometry to define range
        if 'r_mod' in ptrac.geometry:
            r_mod = ptrac.geometry['r_mod']  # outer radius of moderator
        else:
            r_mod = 10.0
            print(f'Geometry not found for {key}. '
                  f'Using default value of {r_mod}')
        # set number of bins
        if w_bin == 0:
            n_bins = int(r_mod*20)
        else:
            n_bins = int(r_mod*w_bin)
        # if radial data present, plot
        if 'r' in df.columns:
            if err is True:
                vals_raw, bin_edges = np.histogram(
                    df['r'], bins=n_bins, weights=df['weight'],
                    range=[0, r_mod + r_DT]
                )
                bin_centers = bin_edges[:-1] + (bin_edges[1] - bin_edges[0])/2
                vals = [i/ptrac.nps/(
                    2*np.pi*bin_centers[idx]*(r_mod+r_DT)/n_bins)
                        for idx, i in enumerate(vals_raw)]
                errs = [np.sqrt(i)/ptrac.nps/(2*np.pi*bin_centers[idx]*(
                    r_mod+r_DT)/n_bins) for idx, i in enumerate(vals_raw)]
                plt.errorbar(bin_centers, vals, yerr=errs,
                             elinewidth=1, capsize=2, drawstyle='steps-mid',
                             lw=3, color=color, label=label)
            else:
                # arbitrarily normalized to 1E9 particles
                plt.hist(df['r'], bins=n_bins,
                         weights=df['weight']/ptrac.nps/(
                             2*np.pi*df['r']*(r_mod+r_DT)/n_bins),
                         histtype='step', range=[0, r_mod + r_DT],
                         lw=3, color=color, label=label)
        else:
            print(f'WARNING: r not found for {key}-{mode}')
    plt.xlim(0, r_mod + r_DT)
    plt.ylim(bottom=0)
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    if True:
        plt.xlabel(r'$r \: [cm]$', labelpad=10)
    else:
        plt.xlabel('RADIAL POSITION (CM)', labelpad=10)
    if data == 'data':
        plt.ylabel(r'$\phi$', labelpad=10)
    elif data == 'data_axial':
        plt.ylabel(r'$\phi_{axial}$', labelpad=10)
    else:
        plt.ylabel('COUNTS', labelpad=10)
    plt.tight_layout()
    plt.legend(loc='upper right')
    if save is True:
        plt.savefig('C:/Users/Avram//Dropbox (MIT)/MIT-PNNL'
                    '/figures/moderator-studies/' + save_label +
                    f'{data}-{mode}' '.pdf',
                    format='pdf')  # , dpi=1200)


def r_plot(mod_dict, data='data', keys=[], labels=[], colors=[], modes=['n'],
           w_bin=0, err=False, save=False, save_label=''):
    """Plot neutron radial position at exit of source-moderator assembly."""
    r_DT = 5
    for mode in modes:
        plt.figure(figsize=(16, 9))
        if keys == []:
            keys = mod_dict[mode].keys()
        if labels == []:
            labels = keys
        if colors == []:
            colors = ['black', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
        for key, label, color in zip(keys, labels, colors):
            ptrac = mod_dict[mode][key]
            df = getattr(ptrac, data)
            # see if r_mod in geometry to define range
            if 'r_mod' in ptrac.geometry:
                r_mod = ptrac.geometry['r_mod']  # outer radius of moderator
            else:
                r_mod = 10.0
                print(f'Geometry not found for {key}. '
                      f'Using default value of {r_mod}')
            # set number of bins
            if w_bin == 0:
                n_bins = int(r_mod*20)
            else:
                n_bins = int(r_mod*w_bin)
            # if radial data present, plot
            if 'r' in df.columns:
                if err is True:
                    vals_raw, bin_edges = np.histogram(df['r'], bins=n_bins,
                                                       weights=df['weight'],
                                                       range=[0, r_mod + r_DT])
                    bin_centers = bin_edges[:-1] + (
                        bin_edges[1] - bin_edges[0])/2
                    vals = [i/ptrac.nps/(
                        2*np.pi*bin_centers[idx]*(r_mod+r_DT)/n_bins)
                        for idx, i in enumerate(vals_raw)]
                    errs = [np.sqrt(i)/ptrac.nps/(
                        2*np.pi*bin_centers[idx]*(r_mod+r_DT)/n_bins)
                        for idx, i in enumerate(vals_raw)]
                    plt.errorbar(bin_centers, vals, yerr=errs,
                                 elinewidth=1, capsize=2,
                                 drawstyle='steps-mid',
                                 lw=3, color=color, label=label)
                else:
                    # arbitrarily normalized to 1E9 particles
                    plt.hist(df['r'], bins=n_bins,
                             weights=df['weight']/ptrac.nps/(
                                 2*np.pi*df['r']*(r_mod+r_DT)/n_bins),
                             histtype='step', range=[0, r_mod + r_DT],
                             lw=3, color=color, label=label)
            else:
                print(f'WARNING: r not found for {key}-{mode}')
        plt.xlim(0, r_mod + r_DT)
        plt.ylim(bottom=0)
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        if True:
            plt.xlabel(r'$r \: [cm]$', labelpad=10)
        else:
            plt.xlabel('RADIAL POSITION (CM)', labelpad=10)
        if data == 'data':
            plt.ylabel(r'$\phi$', labelpad=10)
        elif data == 'data_axial':
            plt.ylabel(r'$\phi_{axial}$', labelpad=10)
        else:
            plt.ylabel('COUNTS', labelpad=10)
        plt.tight_layout()
        plt.legend(loc='upper right')
        if save is True:
            plt.savefig('C:/Users/Avram//Dropbox (MIT)/MIT-PNNL/'
                        'figures/moderator-studies/' + save_label +
                        f'{data}-{mode}' '.pdf',
                        format='pdf')  # , dpi=1200)


def mod_plot(dict_mod, key, data='data_axial', mode='n',
             OD_shield=40, ID_mod=5, bins=[0, 0],
             save=0, save_prefix='modplot-axial'):
    """Plot moderator."""
    df = getattr(dict_mod[mode][key], data)
    nps = dict_mod[mode][key].nps
    r_mod = dict_mod[mode][key].geometry['r_mod']
    plt.figure(figsize=(11, 9))
    if bins != [0, 0]:
        plt.hist2d(df.x, df.y, bins=bins,
                   range=[[-OD_shield/2, OD_shield/2],
                          [-OD_shield/2, OD_shield/2]],
                   weights=df.weight)
    else:
        plt.hist2d(df.x, df.y, bins=[OD_shield, OD_shield],
                   range=[[-OD_shield/2, OD_shield/2],
                          [-OD_shield/2, OD_shield/2]],
                   weights=df.weight/nps*3E9, vmin=0, vmax=100)
    plt.xlabel('X (cm)')
    plt.ylabel('Y (cm)')
    circ1 = patches.Circle((0, 0), ID_mod,
                           facecolor='None', edgecolor='w', lw=2)
    circ2 = patches.Circle((0, 0), ID_mod+r_mod,
                           facecolor='None', edgecolor='w', lw=2)
    rect = patches.Rectangle((-OD_shield/2, -OD_shield/2),
                             OD_shield, OD_shield,
                             facecolor='none', edgecolor='black', lw=2)
    plt.gca().add_patch(circ1)
    plt.gca().add_patch(circ2)
    plt.gca().add_patch(rect)
    plt.xlim(-(OD_shield/2 + 1), OD_shield/2 + 1)
    plt.ylim(-(OD_shield/2 + 1), OD_shield/2 + 1)
    plt.title(f'Neutron flux plot for {key}')
    plt.colorbar()
    if save == 1:
        plt.savefig(save_prefix + f'-{key}.png')


def rad_plot(df, r_lo=0., r_hi=25., weighted=0):
    """Plot radial distribution."""
    plt.figure(figsize=(16, 9))
    if weighted == 1:
        plt.hist(df.r, bins=100, range=[r_lo, r_hi], weights=df.weight/df.r**2,
                 histtype='step', color='black', linewidth=2)
        plt.ylabel(r'COUNTS/CM$^2$')
    else:
        plt.hist(df.r, bins=100, range=[r_lo, r_hi], histtype='step',
                 color='black', linewidth=2)
        plt.ylabel('COUNTS')
    plt.xlim([r_lo, r_hi])
    plt.ylim(bottom=0)
    plt.xlabel('RADIUS (CM)')
    plt.tight_layout()


def hist_plot(df, mode='d', e_lo=1.0, e_hi=100.0, e_res=1.0,
              t_lo=0.0, t_hi=10.0, t_res=0.05,
              d_lo=85., d_hi=200., d_res=2.5):
    """Plot histogram."""
    plt.figure(figsize=(16, 9))
    if mode == 'd':
        plt.hist2d(df.E, df.d,
                   bins=[int((e_hi - e_lo)/e_res), int((d_hi-d_lo)/d_res)],
                   range=[[e_lo, e_hi], [d_lo, d_hi]],
                   weights=df.weight)
        plt.ylabel('TOF DISTANCE (mm)', labelpad=20)
    elif mode == 't':
        plt.hist2d(df.E, df.t,
                   bins=[int((e_hi - e_lo)/e_res), int((t_hi-t_lo)/t_res)],
                   range=[[e_lo, e_hi], [t_lo, t_hi]],
                   weights=df.weight)
        plt.ylabel('MODERATION TIME (us)', labelpad=20)
    plt.xlabel('ENERGY (eV)', labelpad=10)


def mod_dist_plot(vals, mode, e_lo=1.0, e_hi=100.0, e_res=1.0, unit='mm'):
    """Print distribution of average."""
    e_lin = np.linspace(e_lo, e_hi, int((e_hi - e_lo)/e_res) + 1)
    plt.figure(figsize=(16, 9))
    if unit == 'mm':
        pass
    elif unit == 'cm':
        vals /= 10
    plt.errorbar(x=e_lin + e_res/2, y=vals[:, 0], yerr=vals[:, 1],
                 linestyle='none', elinewidth=1, capsize=2,
                 marker='s', color='black')
    plt.xlabel('FINAL NEUTRON ENERGY (eV)')
    if mode == 'd':
        if unit == 'mm':
            plt.ylabel('TOF DISTANCE (mm)')
        elif unit == 'cm':
            plt.ylabel('TOF DISTANCE (cm)')
    elif mode == 't':
        plt.ylabel('MODERATION TIME (us)')
    plt.xlim([e_lo, e_hi])
    plt.ylim([min(vals[:, 0])/1.1, 1.1*max(vals[:, 0])])
    plt.tight_layout()
    plt.xticks(np.arange(0, 110, step=10))
