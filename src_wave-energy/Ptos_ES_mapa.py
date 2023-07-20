import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs 
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
import cmocean as cmo

def map_features(ax):
    ax.add_feature(COASTLINE)
    ax.add_feature(BORDERS, edgecolor='#383838')
    return ax

def Brazil_states(ax):    
    states = NaturalEarthFeature(category='cultural', scale='50m',
                             facecolor=google_maps_colors[2], name='admin_1_states_provinces_lines')
    _ = ax.add_feature(states, edgecolor='#383838')
    
    cities = NaturalEarthFeature(category='cultural', scale='50m', facecolor=google_maps_colors[2],
                                  name='populated_places')
    _ #= ax.add_feature(cities, edgecolor='#383838')
    
def grid_labels_params(ax,i):
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    if i not in [0,3]:
        gl.left_labels = False
    gl.xlabel_style = {'size': 20, 'color': '#383838'}
    gl.ylabel_style = {'size': 20, 'color': '#383838'}
    ax.spines['geo'].set_edgecolor('#383838')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return ax

#Carregando a batimetria
ds = xr.open_dataset('/p1-nemo/rtecchio/Dados/GEBCO/gebco_2023_costa_s_se.nc')
prof2 = ds['elevation'][:]
prof_filtred = prof2.where(prof2 <= 0)
prof_filtred.min()
prof_filtred.max()
# Arquivo com a posição dos pontos na costa do ES
ptos = pd.read_csv('/p1-nemo/rtecchio/teste_chico/ptos_iso_100_BRANCO.csv',sep=",")
lat=ptos['lat']
lon=ptos['lon']

###############################################
#Plotando
google_maps_colors = ['#ebebeb', '#c7c7c7', '#b0b0b0', '#999999', '#7e7e7e', '#666666', '#4d4d4d', '#333333', '#1a1a1a', '#0c7b12']

# make figure with points 
data_min = prof_filtred.min()
data_max = prof_filtred.max()
interval = 50
levels = np.linspace(data_min,data_max,interval)

proj = ccrs.PlateCarree() 
plt.close('all')
fig = plt.figure(constrained_layout=False,figsize=(12,12))
ax = fig.add_subplot(111, projection=proj)
ax.set_extent([-42,-37,-22,-17.7])
grid_labels_params(ax,0)
Brazil_states(ax)
#cf1 = ax.contourf(iso100.lon,iso100.lat,iso100)
#cf1= ax.contourf(prof_filtred.lon,prof_filtred.lat,prof_filtred, cmap='cmo.deep_r',)
cf1= ax.contourf(prof_filtred.lon,prof_filtred.lat,prof_filtred,levels=levels, cmap='cmo.deep_r',)
ax.plot(lon, lat, 'o',color='k',  )
ax.text(lon[0]+0.2,lat[0]-0.2, '01',  color='black', fontsize=14, transform=proj)
ax.text(lon[1]+0.2,lat[1]-0.2, '02',  color='black', fontsize=14, transform=proj)
ax.text(lon[2]+0.2,lat[2]-0.2, '03',  color='black', fontsize=14, transform=proj)
# cb=plt.colorbar(cf1, ax=ax, shrink=0.8, aspect=15)
# cb.ax.tick_params(labelsize=12)
# cb.set_label('Profundidade (m)',fontsize=12)
map_features(ax)
#plt.show()
plt.savefig('Pontos_branco.png', dpi=300)