
# coding: utf-8

import numpy as np

import yt
from yt.utilities.physical_constants import mp, kb
from yt.units.yt_array import YTQuantity

from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import pylab
import matplotlib.pyplot as plt

from plutopy import create_fields

def mass_loss(ts, file_name='mass_loss.txt'):
    """
    Synopsis
    --------
    This function calculates the mass-loss from 
    the planet, star and system.

    Parameters
    ----------
    ts: list-like
        list of data sourses that forms a time-serise.

    filename: string-like
        name for output file.     

    Returns
    -------
    None: does not return a variable.
        Writes result to file 'outputs/'+filename

    TODO
    ----
    - Make general ont dependant on the planet always 
      being at 0.047 au.
    """

    # Initial list to stor mass-loss values.
    planet_mass_loss = []
    star_mass_loss = []
    system_mass_loss = []
    time = []

    Rsun = 6.955e10 # [cm]

    for i, ds in enumerate(ts):
        
        create_fields(ds)
        ds.periodicity = (True, True, True)

        # Planet
        cp = yt.YTArray([0.047, 0.0, 0.0], 'au').convert_to_units('cm')
        rp = YTQuantity(2.0, 'Rsun').in_cgs()
        sp = ds.sphere(cp, rp)
        pr = YTQuantity(0.75, "Rsun").in_cgs() 
        surfp = ds.surface(sp, 'radius_planet', pr)
        density_flux_p = surfp.calculate_flux("velocity_x", 
                                              "velocity_y", 
                                              "velocity_z", 
                                              "density"
                                              )
        density_flux_p *= Rsun**2

        #Star
        cs = yt.YTArray([0.0, 0.0, 0.0], 'au').convert_to_units('Rsun')
        rs = YTQuantity(5.0, 'Rsun').in_cgs()
        ss = ds.sphere(cs, rs)
        sr = YTQuantity(2.0, "Rsun").in_cgs()
        surfs = ds.surface(ss, 'radius', sr)
        density_flux_s = surfs.calculate_flux("velocity_x", 
                                              "velocity_y", 
                                              "velocity_z", 
                                              "density"
                                              )
        density_flux_s *= Rsun**2

        #System
        cs = yt.YTArray([0.0, 0.0, 0.0], 'au').convert_to_units('Rsun')
        rs = YTQuantity(32.0, 'Rsun').in_cgs()
        ss = ds.sphere(cs, rs)
        sr = YTQuantity(31.0, "Rsun").in_cgs()
        surfs = ds.surface(ss, 'radius', sr)
        density_flux_system = surfs.calculate_flux("velocity_x", 
                                                   "velocity_y", 
                                                   "velocity_z", 
                                                   "density"
                                                   )
        density_flux_system *= Rsun**2

        planet_mass_loss.append(density_flux_p)
        star_mass_loss.append(density_flux_s)
        system_mass_loss.append(density_flux_system)
        time.append(ds.current_time.convert_to_units('s'))

        print("file number:", i+1)
        print("Planet mass-loss =", density_flux_p)
        print("Star mass-loss   =", density_flux_s)
        print("System mass-loss =", density_flux_system)

        np.savetxt('outputs/'+file_name, (time, star_mass_loss, planet_mass_loss, system_mass_loss))


def mass_loss_plot(input_file='mass_loss.txt', output_file='mass_loss.png'):
    """
    Synopsis
    --------
    This function plots the results from the mass-loss 
    calculations in function 'mass_loss'.

    Parameters
    ----------
    input_file: string-like
        name on file which mass-losses are saved to.

    output_file: string-like
        name for output file.     

    Returns
    -------
    None : does not return a variable.
        Produces a image in png/pdf format.

    TODO
    ----
    - None :D
    """

    time, star_mass_loss, planet_mass_loss, system_mass_loss = np.loadtxt(
        'outputs/'+input_file)

    skip = 0
    fig = plt.figure()
    fig=plt.figure(figsize=(10, 8), dpi= 80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)

    ax.plot(time[planet_mass_loss > 0.0]/1.0e+3, 
        planet_mass_loss[planet_mass_loss > 0.0], 'g', label='Planet mass-loss lowUV')
    ax.plot(time/1.0e+3, star_mass_loss, 'k-.', label='Star mass-loss')
    ax.plot(time/1.0e+3, system_mass_loss, 'k--', label='System mass-loss')

    ax.set_xlabel('Time $[\mathrm{ks}]$')
    ax.set_ylabel('Mass-loss $[\mathrm{g/s}]$')
    ax.set_ylim(1.0e9, 1.0e13)
    ax.set_yscale("log", nonposy='clip')
    ax.legend()
    plt.savefig('plots/'+output_file)
    plt.close()




