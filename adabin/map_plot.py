
from adabin.adabin import adap_bin, make_seg_map
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import (MinMaxInterval, PercentileInterval, ImageNormalize, simple_norm, LinearStretch,SqrtStretch ,LogStretch)
from random import random
from matplotlib.colors import LinearSegmentedColormap

def _display_map(ax, signal, noise, min_SN,
                log=True, vrange=(-3, 0),
                cmap=plt.cm.RdYlBu_r, fontsize=20):
    img = np.flipud(np.where(signal / noise > min_SN, signal, np.nan))
    h, w = signal.shape
    
    if log:
        img = np.log10(img)

    if min_SN < 0.:
        im = ax.imshow(img, cmap=plt.cm.get_cmap('RdYlBu', int(np.ptp(img))))
    else:
        print("[display_map]: %d pixels with S/N > %d.\n" %(
              np.sum(signal / noise > min_SN), min_SN))
        Interval=PercentileInterval(99.8)
        norm1 = ImageNormalize(img,Interval,stretch=LinearStretch())
        im = ax.imshow(img, cmap=cmap, norm=norm1)
    
    return im


def recon_plot(signal_map, noise_map, min_SN, eline='SII6731', vrange=(-3, 0), fontsize=20):

    ticks = range(vrange[0], vrange[1] + 1)
    ticklabels = []
    for tick in ticks:
        ticklabels.append('10$^{%d}$' %(tick))

    signal = np.where(signal_map > 0., signal_map, 0.)
    noise = np.where(noise_map > 0., noise_map, np.mean(noise_map))

    rec_s, rec_n, maps = adap_bin(signal, noise, min_SN)
    
    fig = plt.figure(figsize=(24, 9))
    plt.subplots_adjust(left=.04, bottom=.08, right=.98, top=.8, wspace=.15)

    ax = plt.subplot(131)
    im = _display_map(ax, signal, noise, min_SN)
    cbar_ax = fig.add_axes([0.04, 0.92, 0.28, 0.05])
    cbar = plt.colorbar(im, orientation="horizontal", cax=cbar_ax, ticks=ticks)
    cbar.ax.set_xticklabels(ticklabels)
    cbar.set_label('Line flux ($10^{-16}$ erg)', size=20)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xlabel('x (arcsec)', fontsize=fontsize)
    ax.set_ylabel('y (arcsec)', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)

    ax = plt.subplot(132)

    s, n, im = make_seg_map(signal, noise, min_SN)
    Interval=PercentileInterval(95.8)
    norm1 = ImageNormalize(im,Interval,stretch=SqrtStretch())
    
    img = ax.imshow(im,origin='lower', cmap='RdYlBu', norm=norm1)
    
    cbar_ax = fig.add_axes([0.354, 0.92, 0.28, 0.05])
    cbar = plt.colorbar(img, orientation="horizontal",cax=cbar_ax)
    cbar.set_label('Bin number', size=20)
    cbar.ax.tick_params(labelsize=15)
    
    ax.set_xlabel('x (arcsec)', fontsize=fontsize)
    ax.set_ylabel('y (arcsec)', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)

    ax = plt.subplot(133)
    im = _display_map(ax, rec_s, rec_n, min_SN)
    cbar_ax = fig.add_axes([0.695, 0.92, 0.28, 0.05])
    cbar = plt.colorbar(im, orientation="horizontal", cax=cbar_ax, ticks=ticks)
    cbar.ax.set_xticklabels(ticklabels)
    cbar.set_label('Line flux ($10^{-16}$ erg)', size=20)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xlabel('x (arcsec)', fontsize=fontsize)
    ax.set_ylabel('y (arcsec)', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)

    plt.show()


