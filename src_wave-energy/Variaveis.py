import xarray as xr
import numpy as np
import math as m
import pandas as pd
from scipy import signal
import os
import sys
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from xarray import DataArray
#import matplotlexitib.ticker as plticker
import datetime
from glob import glob

#from hovmoller_climatologia import dados_para_xarray


def wef(Hs, T, d, lat, dens, dir_wave, dir_coast):
    """
    This function calculates the vector Wave Energy Flux (WEF) of gravity waves
    that propagate in an air-water interface (such as the surface of the ocean).
    It assumes the Linear (Airy) Wave Theory for "intermediate waters". This is
    the general case that leads to the simplified cases "deep waters" and "shallow
    waters". Therefore, this function is valid for any water depths that are in
    agreement with the Linear (Airy) Wave Theory. Since the formulas are adjusted
    so that the significant wave height (Hs) is used as an input (instead of a
    sinusoidal height, H), this code is suitable for actual ocean waves. It does
    not consider any refraction or diffraction phenomena that may occur between
    the measurement point and a supposed coast. The component of the WEF vector
    that is perpendicular to the coast (output "Pper") is calculated via simple
    vector projection. All outputs are given in units kJ/(m*s) = kW/m.

    Usage from an external code at the same directory:
    from wef import wef
    dict = wef(Hs, T, d, lat, dens, dir_wave, dir_coast)
    (P, Pper, Px, Py) = (dict['P'], dict['Pper'], dict['Px'], dict['Py'])

    Inputs:
    1.        Hs: Wave significant height in meters.
    2.         T: Wave period in seconds.
    3.         d: Local depth in meters.
    4.       lat: Local latitude in degrees and decimal degrees (used to
                  calculate the local effective gravity due to Earth's
                  rotation). Positive: northern hemisphere. Negative: southern
                  hemisphere.
    5.      dens: Local water density in kg/m3.
    6.  dir_wave: Direction FROM WHICH the wave is propagating, measured in
                  degrees from the geographical north towards east (i.e.
                  clockwise). Example: a wave that propagates from east to west
                  has dir_wave = 90.
    7. dir_coast: Direction TO WHICH the coast is facing, measured in degrees
                  from the geographical north towards east (i.e. clockwise).
                  Example: a beach that faces east has dir_coast = 90.

    Outputs:
    1.    P: Magnitude of the vector Wave Energy Flux (WEF) in kW/m.
    2. Pper: Magnitude, in kW/m, of the component of the vector Wave Energy Flux
             (WEF) that is perpendicular to the coast. If positive: wave
             propagates towards the coast. If negative: wave propagates away
             from the coast.
    3.   Px: Zonal (east-west) component of the vector Wave Energy Flux (WEF) in
             kW/m. If positive: wave towards east. If negative: wave towards
             west.
    4.   Py: Meridional (north-south) component of the vector Wave Energy Flux
             (WEF) in kW/m. If positive: wave towards north. If negative: wave
             towards south.

       Author: Thalles A. A. Araujo (https://orcid.org/0000-0002-6204-6886)
      Contact: <thalles.araujo@alumni.usp.br>
     Creation: 2019-Feb-27
    Last update: 2022-Mar-26

    Author's note: This code was developed as part of the process of submitting
    the scientific article in <DOI>. The author encourages its free use while
    respectfully requests for its source to be properly cited. The author takes no
    responsibility for its use by others. Any constructive criticism is greatly
    appreciated. This same algorithm can be downloaded for both Python and
    MATLAB/Octave programming languages at <https://github.com/thallesaaa>, where
    other codes are freely available too.
    """

    # Importing required libraries:

    # Defining required physical constants:
    gconst  = (6.6743015)*(10.0)**(-11) # Gravitational constant in SI.
    emass   = (5.9736)*(10.0)**(24)     # Earth's mass in kg.
    eradius = (6.371)*(10.0)**(6)       # Earth's radius in m.
    omega   = 2.0*np.pi/86164.0905      # Earth's angular speed in rad/s.

    # Conversion of latitude from degrees to radians:
    lat = np.pi*lat/180

    # Effective gravity (considers centripetal force due to Earth's rotation):
    g = gconst*emass/(eradius**2) - (omega**2)*eradius*(np.cos(lat)**2)

    # Conversion of dir_wave and dir_coast due to the trigonometric convention:
    dir_wave_converted = -dir_wave - 90
    dir_coast_converted = -dir_coast + 90
    # This angle conversion is needed because an azimuth angle has its zero on
    # north and grows clockwise whereas a conventional angle in the trigonometric
    # circle has its zero on the right ("east") and grows anticlockwise.

    # Vectorial wave azimuth (describe direction TO which wave is propagating):
    dx_wave = np.cos(dir_wave_converted*np.pi/180)
    dy_wave = np.sin(dir_wave_converted*np.pi/180)

    # Vectorial coast azimuth (describe direction TO which the coast is facing):
    dx_coast = np.cos(dir_coast_converted*np.pi/180)
    dy_coast = np.sin(dir_coast_converted*np.pi/180)

    # Calculating the wavelength (L) via iteration of the dispersion relation:
    L1 = g*T**2/(2*np.pi) # Initial value for L (its deep water value).
    L2 = (g*T**2/(2*np.pi))*np.tanh(2*np.pi*d/L1) # First iteration.
  
    # All of the other iterations, until L differs to its real value by 1 mm:
    while (np.abs(L1-L2) > 0.001):
        L1 = (g*T**2/(2*np.pi))*np.tanh(2*np.pi*d/L2)
        L2 = (g*T**2/(2*np.pi))*np.tanh(2*np.pi*d/L1)
    L = L2 # Final value for the wavelength (L).

    # Calculating the wavenumber:
    K = 2*np.pi/L # Magnitude of the wavenumber vector.
    k = K*dx_wave # Zonal (east-west) component of the wavenumber vector.
    l = K*dy_wave # Meridional (north-south) component of the wavenumber vector.

    # Calculating the wave angular frequency:
    omega = 2*np.pi/T

    # Calculating the wave energy density:
    E = (dens*g*Hs**2)/16

    # Calculating the phase velocity:
    cx = omega*k/(K**2)
    cy = omega*l/(K**2)

    # Calculating vector Wave Energy Flux (WEF) in J/(m.s) = W/m:
    Px = E*(cx/2)*(1 + (2*K*d)/np.sinh(2*K*d))
    Py = E*(cy/2)*(1 + (2*K*d)/np.sinh(2*K*d))
    # Conversion to kJ/(m.s) = kW/m:
    Px = Px/1000 # Zonal (east-west) component of the WEF vector.
    Py = Py/1000 # Meridional (north-south) component of the WEF vector.

    # Magnitude of the vector Wave Energy Flux (WEF) in kJ/(m.s) = kW/m:
    P = np.sqrt(Px**2 + Py**2)

    # Magnitude, in kJ/(m.s) = kW/m, of the component of the WEF vector that is
    # perpendicular to the coast:
    Pper = -(Px*dx_coast + Py*dy_coast)
    Pperx = -Pper*dx_coast # Its zonal component (not an output).
    Ppery = -Pper*dy_coast # Its meridional component (not an output).
    Hsp = ((P/5)**0.5)
    return {'Hs': Hs, 'Tp': T, 'Dir': dir_wave,'P': P, 'Pper': Pper, 'Px': Px, 'Py': Py, 'Ppar':Ppery, 'Hsp':Hsp}


#tempo_modelo = slice('1993-01-01T00:00:00','2019-12-31T23:00:00') 
ds = xr.open_mfdataset('/p1-nemo/mbonjour/onda_dados/waverys_1993_2019.nc', parallel=True)
lat=ds['latitude']
lon=ds['longitude']

tan_ptos=pd.read_csv('/p1-nemo/rtecchio/teste_chico/pontos_branco_dir_perp.csv', sep=';', decimal=',')



# Constantes
dens = 1030
d=100

pontos = tan_ptos['PONTO'].unique()
for ponto in pontos: 
    dados =  tan_ptos.loc[tan_ptos['PONTO'] == ponto]
    lat_ponto, lon_ponto = float(dados.lat), float(dados.lon)
    dir_coast= float(dados.DIR_NORM)
    print('coletando dados para o ponto',ponto,'de lat/lon:',lat,lon)
    iwmo_ptos= ds.sel(latitude=lat_ponto, longitude=lon_ponto, method='nearest')
    # dados horarios de Hs e T e Dir
    tempo_escolhido= slice('2010-01-01T00:00:00','2019-12-31T23:00:00')
    tempo= iwmo_ptos['time'].sel(time=tempo_escolhido).squeeze()
    Hs = iwmo_ptos['VHM0'].squeeze()
    T = iwmo_ptos['VTPK'].squeeze()
    dir_wave = iwmo_ptos['VMDR'].squeeze()
    for t in tempo.time:
      time = t.values
      #print(time)
      energy_dict = wef(float(Hs.sel(time=t)), float(T.sel(time=t)), d, lat_ponto, dens, float(dir_wave.sel(time=t)), dir_coast)
      energy_list = []
      for key in energy_dict.keys():
        energy_list.append(float(energy_dict[key]))
      if t ==  tempo.time[0]:
        df = pd.DataFrame({'Hs':energy_list[0],'Tp':energy_list[1], 'Dir':energy_list[2],'P':energy_list[3],  'PPer':energy_list[4],'Px':energy_list[5],'Py':energy_list[6], 'Ppar':energy_list[7], 'Hsp':energy_list[8]},index=[time])
      else:
        tmp = pd.DataFrame({'Hs':energy_list[0],'Tp':energy_list[1], 'Dir':energy_list[2],'P':energy_list[3],  'PPer':energy_list[4],'Px':energy_list[5],'Py':energy_list[6], 'Ppar':energy_list[7], 'Hsp':energy_list[8]},index=[time])
        df = pd.concat([df,tmp])
    df['ponto'] = ponto
    
    if ponto == pontos[0]:
        energia_ptos = df
    else:
         energia_ptos = pd.concat([energia_ptos,df])
 
energia_ptos['time'] = energia_ptos.index
energia_ptos.index = range(len(energia_ptos))
energia_ptos.to_csv('teste_chico/variáveis_branco.csv')