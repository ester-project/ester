#!env python
# -*- coding: utf-8 -*-

"""
Created on Thu Apr  4 16:47:05 2019

@author: Kévin Bouchaud
updated by Michel Rieutord
made compatible python2.7

 Usage:
 run SED
 a=Star('model_filename',inclination in degrees)
 a.Teff_mean() gives the mean Teff

 all Ester model data are in a.mod.Your_Param

"""

import os

import numpy as np
from numpy import pi, cos, sin, sqrt
from scipy.integrate import quad
import itertools
import matplotlib.pyplot as plt
from matplotlib import ticker
import warnings

from ester import star2d

class Star:
    """Class representing the visible grid of the model star"""

    # Important constants, in cgs units
    k    = 1.38064852e-16                  # Boltzmann's constant
    h    = 6.62607015e-27                  # Planck's constant
    c    = 29979245800.                    # Speed of light in vacuum
    pc   = 648000 / np.pi * 14959787070000 # One parsec in centimeters
    dist = 5.13 * pc                       # Distance to the star in centimeters
    Msun = 1.9891e33
    Rsun = 6.95508e10
    Lsun = 3.8396e33

    def __init__(self, mod, incl, nth):
        self.incl = incl*pi/180
        self.mod  = star2d(mod)
        self.nth  = nth
        self.init_grid()
        self.visgrid()

    def init_grid(self):
        """Initiate physical and grid parameters"""
        self.M                  = self.mod.M / Star.Msun
        self.L                  = self.mod.L / Star.Lsun
        self.nphi0              = 5
        self.dth                = pi / self.nth
        self.theta              = [self.dth/2. + j * self.dth for j in range(self.nth)]
        self.Rs                 = self.mod.Rp
        self.eval_func          = [self.mod.leg_eval_matrix(j) for j in self.theta]
        self.eval_func_antisym  = [self.mod.leg_eval_matrix_antisym(j) for j in self.theta]
        self.r                  = np.dot(self.mod.r[-1, 1:-1], self.eval_func).ravel()
        self.rt                 = np.dot(self.mod.rt[-1, 1:-1], self.eval_func_antisym).ravel()
        self.Teff               = np.dot(self.mod.Teff[:, 1:-1], self.eval_func).ravel()
        self.logg               = np.array([np.log10(i) for i in np.dot(self.mod.gsup[:, 1:-1],
                                            self.eval_func)[0]]).ravel()
        self.ds0                = self.ds_func(0, self.nphi0)
        self.dphi, self.nphi    = self.grid_phi(nphi0=self.nphi0, nth=self.nth, ds0=self.ds0)
        self.phi                = [[self.dphi[i] / 2. + j * self.dphi[i] for j in\
                                    range(int(self.nphi[i]))] for i in range(self.nth)]
        self.ds                 = [self.ds_func(i, self.nphi[i]) for i in range(self.nth)]
        self.ngrid              = sum(self.nphi)
        self.init_mu()
        self.flat()

    def ds_func(self, i, j):
        """
        Function to compute the surface element's area at given theta (theta[i]) and dphi (2*pi/j).
        """
        return self.r[i]**2 * sqrt(1 + (self.rt[i]**2/self.r[i]**2)) * np.sin(self.theta[i])\
               * self.dth * 2 * pi / j

    def grid_phi(self, nphi0, nth, ds0):
        """
        Compute the phi step at every theta so that the surface elements' areas are as
        homogeneous as possible across the star.
        """
        dphi = [2 * pi / nphi0]
        nphi = [nphi0]
        for i in range(1, nth):
            temp = float('inf')
            j = 1
            while abs(self.ds_func(i, j) - ds0) <= temp:
                temp = abs(self.ds_func(i, j) - ds0)
                index = j
                j += 1
            dphi.append(2 * pi / index)
            nphi.append(int(index))

        return np.array(dphi), np.array(nphi)

    def init_mu(self):
        """
        Compute mu, cosine of the angle between the normal to the surface and the line of sight
        (depends on inclination angle, and necessary to compute the visible grid).
        """
        self.mu       = [[(cos(self.phi[i][j])*sin(self.incl)*(self.r[i]*sin(self.theta[i])
                          - self.rt[i]*cos(self.theta[i])) + cos(self.incl)
                          * (self.r[i]*cos(self.theta[i]) + self.rt[i]*sin(self.theta[i])))
                          / (self.r[i]*np.sqrt(1 + (self.rt[i]**2 / self.r[i]**2)))
                          for j in range(int(self.nphi[i]))]
                         for i in range(self.nth)]
        """
        self.mu       = [[round((cos(self.phi[i][j])*sin(self.incl)*(self.r[i]*sin(self.theta[i])
                          - self.rt[i]*cos(self.theta[i])) + cos(self.incl)
                          * (self.r[i]*cos(self.theta[i]) + self.rt[i]*sin(self.theta[i])))
                          / (self.r[i]*np.sqrt(1 + (self.rt[i]**2 / self.r[i]**2))), 6)
                          for j in range(int(self.nphi[i]))]
                         for i in range(self.nth)]
        """
        self.mu_flat  = np.array(list(itertools.chain(*self.mu)))
        self.vis_mask = np.where(self.mu_flat >= 0.)
        self.mu_vis = self.mu_flat[self.vis_mask]

    def flat(self):
        """Flatten all arrays"""
        ds   = [list(itertools.repeat(self.ds[i], self.nphi[i])) for i in range(self.nth)]
        Teff = [list(itertools.repeat(self.Teff[i], self.nphi[i])) for i in range(self.nth)]
        logg = [list(itertools.repeat(self.logg[i], self.nphi[i])) for i in range(self.nth)]
        self.ds_flat   = np.array(list(itertools.chain(*ds)))
        self.Teff_flat = np.array(list(itertools.chain(*Teff)))
        self.logg_flat = np.array(list(itertools.chain(*logg)))
        return

    def visgrid(self):
        """Only keep values for the visible surface of the star"""
        self.ds_vis   = self.ds_flat[self.vis_mask]
        self.Teff_vis = self.Teff_flat[self.vis_mask]
        self.logg_vis = self.logg_flat[self.vis_mask]

    def Lapp(self):
        inc=self.incl*180/pi
        return self.mod.apparent_luminosity(inc)
    def Teff_mean(self):
        return sum(self.Teff_vis*self.mu_vis*self.ds_vis)/sum(self.mu_vis*self.ds_vis)

    def surf_proj(self):
        return sum(self.mu_vis*self.ds_vis)

    def logg_mean(self):
        return sum(self.logg_vis*self.mu_vis*self.ds_vis)/sum(self.mu_vis*self.ds_vis)

    def B(self, l, T):
        return 2*Star.h*Star.c**2 / (l**5*(np.exp(Star.h*Star.c/(l*Star.k*T))-1))

    def SED(self, wavelength=np.linspace(1e-5, 2e-2, 10000)):
        flux = []
        for T in self.Teff_vis:
            flux.append(self.B(wavelength, T))
        return wavelength, np.array([sum(f*self.mu_vis*self.ds_vis*self.Rs**2/Star.dist**2)
                                     for f in np.array(flux).T])

    def plot_SED(self):
        wave, sed = self.SED()
        fig, ax = plt.subplots(figsize=(19.20, 10.80))

        ax.loglog(wave*1e4, sed*1e-4)

        ax.set_title('SED')
        ax.tick_params(labelsize=16, direction='inout', which='both', length=5, width=1)
        ax.set_xlabel('$\lambda$ ($\mu$m)', fontsize=20)
        ax.set_ylabel('F$_\lambda$ (erg$\cdot$s$^{-1}\cdot$cm$^{-2}\cdot\mu$m$^{-1}$)', fontsize=20)
        ax.yaxis.set_minor_locator(plt.LogLocator(base=10, numticks=15))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
        fig.show()
        return


def plot(model, incl, save=None):
    global directory

    # Check that model exists in the right directory
    try:
        assert os.path.exists(model)
    except AssertionError:
        print("attention file not found")
        #raise FileNotFoundError(f"No such file in directory /{directory} : {model}")

    # Create the stellar visible grid
    star = Star(mod=model, incl=incl)

    # Get observed SED
    wave_obs, sed_obs = np.loadtxt('/home/rieutord/Ester/postprocessing/SED/SED',unpack=True)
                                   #skiprows=1, unpack=True)

    # compute theoretical SED
    wave_mod, sed_mod = star.SED()

    # Plot theoretical SED of model (line) alongside observed SED (points)
    fig, ax = plt.subplots(figsize=(19.20, 10.80))
    ax.loglog(wave_mod*1e4, sed_mod*1e-4)
    ax.loglog(wave_obs, sed_obs, 'o')

    fig.suptitle('SED', fontsize=24)
    ax.set_title('M = {:.5g} Msun, Z = {:.3f}, Xc = {:.5g}, i = {:.2f}'.format(star.mod.M/Star.Msun,star.mod.Z[0, 0],star.mod.Xc,star.incl), fontsize=18)
    ax.set_xlabel('$\lambda$ ($\mu$m)', fontsize=20)
    ax.set_ylabel('F$_\lambda$ (erg$\cdot$s$^{-1}\cdot$cm$^{-2}\cdot\mu$m$^{-1}$)', fontsize=20)
    ax.yaxis.set_minor_locator(plt.LogLocator(base=10, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(labelsize=16, direction='inout', which='both', length=5, width=1)

    # Save figure (filename must be given as save='filename' in function arguments)
    if save:
        fig.savefig(fname='/home/rieutord/Ester/postprocessing/SED/'+save, bbox_inches='tight')

    fig.show()
    return

def plot_all(incl, save=None):
    global models

    rows = 2
    cols = 2

    fig, axes = plt.subplots(nrows=rows, ncols=cols, sharex=True, sharey=True,
                             figsize=(19.20, 10.80))
    fig.subplots_adjust(wspace=0., hspace=0.)

    for i, model in enumerate(models):
        ax = axes.flatten()[i]
        star = Star(model, incl)

        wave_obs, sed_obs = np.loadtxt('/home/rieutord/Ester/postprocessing/SED/SED', unpack=True)
        wave_mod, sed_mod = star.SED()
        ax.loglog(wave_mod*1e4, sed_mod*1e-4)
        ax.loglog(wave_obs, sed_obs, 'o')
        ax.yaxis.set_minor_locator(plt.LogLocator(base=10, numticks=15))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
        ax.tick_params(labelsize=16, direction='inout', which='both', length=5, width=1)

        if ax.colNum != 0:
            ax.tick_params(left=False)
        else:
            ax.set_ylabel('F$_\lambda$ (erg$\cdot$s$^{-1}\cdot$cm$^{-2}\cdot\mu$m$^{-1}$)',
                          fontsize=20)

        if ax.rowNum != 1:
            ax.tick_params(bottom=False)
        else:
            ax.set_xlabel('$\lambda$ ($\mu$m)', fontsize=20)

        text = ('M  = {star.mod.M/Star.Msun:.2f} M$\odot$\nZ  = {star.mod.Z[0, 0]:.3f}\n'
                'Xc = {star.mod.Xc:.3f}')
        props = dict(boxstyle='round', facecolor='white', alpha=0.7)
        # place a text box in lower left in axes coords
        ax.text(0.75, 0.90, text, transform=ax.transAxes, verticalalignment='top',
                bbox=props, fontdict=dict(fontsize=14))

    if save:
        fig.savefig(fname='/home/rieutord/Ester/postprocessing/SED/'+save, bbox_inches='tight')
    fig.show()
    return

def plot_one_for_all(incl, save=None, res=False, error=0.1):

    global models

    fig, ax = plt.subplots(figsize=(10.80, 10.80))

    wave_obs, sed_obs = np.loadtxt('/home/rieutord/Ester/postprocessing/SED/SED', unpack=True)

    ax.loglog(wave_obs, sed_obs, 'o')

    ax.set_xlabel('$\lambda$ ($\mu$m)', fontsize=20)

    ax.set_ylabel('F$_\lambda$ (erg$\cdot$s$^{-1}\cdot$cm$^{-2}\cdot\mu$m$^{-1}$)', fontsize=20)
    ax.yaxis.set_minor_locator(plt.LogLocator(base=10, numticks=15))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(labelsize=16, direction='inout', which='both', length=5, width=1)


    linestyles = itertools.cycle([(0, (5, 5)), 'dashdot', (0, (3, 5, 1, 5, 1, 5)), 'dotted', '-'])

    # Print reduced Chi2 of fit with arbitrary uncertainty on observed data
    if res:
        print('Reduced chi2 ({error*100}% artificial error on data):\n')

    for i, model in enumerate(models):
        star = Star(mod=model, incl=incl)

        wave_mod, sed_mod = star.SED()
        ax.loglog(wave_mod*1e4, sed_mod*1e-4, ls=next(linestyles),
                  label=('M = {star.mod.M/Star.Msun:.2f} M$\odot$, Z = {star.mod.Z[0, 0]:.3f}, '
                         'Xc = {star.mod.Xc:.2f}, $\overline{{\mathrm{{T}}}} = '
                         '{int(star.Teff_mean())} $K'),
                         zorder=0)

        if res:
            wave_mod_chi2, sed_mod_chi2 = star.SED(wavelength=wave_obs*1e-4)
            sum_squared_res = sum(((sed_mod_chi2*1e-4 - sed_obs)/(error*sed_obs))**2)
            chi2_red = sum_squared_res / len(sed_mod_chi2)
            print('model {star.mod.M/Star.Msun:.2f}: {chi2_red}')

    ax.legend(fontsize=16, loc='lower left')

    if save:
        fig.savefig(fname='./'+save, bbox_inches='tight')
    fig.show()
    return

if __name__ == '__main__':
    path  = '/home/rieutord/Ester/postprocessing/SED/'
    directory = os.path.basename(os.path.normpath(path))

    model70 = path + '2d_m1.70_z0.0120_xc0.603_obk74.4.h5'
    model75 = path + '2d_m1.75_z0.0161_xc0.759_obk74.4.h5'
    model80 = path + '2d_m1.80_z0.0201_xc0.900_obk74.4.h5'

    models = [model70, model75, model80]



