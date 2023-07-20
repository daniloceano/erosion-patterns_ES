import xarray as xr
import numpy as np
import pandas as pd
from geopy.distance import geodesic
import math


    
#extraindo os pontos da iso de 100
ds = xr.open_dataset('/p1-nemo/rtecchio/Dados/GEBCO/gebco_2023_costa_s_se.nc')
prof2 = ds['elevation'][:]
prof_filtred = prof2.where(prof2 <= 0)
prof_filtred.min()
prof_filtred.max()
#extraindo lat e lon da iso de 100 metros
lats = ds.variables['lat'][:]
lons = ds.variables['lon'][:]
prof = ds['elevation'][:]
iso100 = prof.where(prof == -100)
iso100_stacked = iso100.stack(x=['lat','lon'])
# Selecionando pontos e espaçando com a mesma distância 

iso100_coords = iso100_stacked[iso100_stacked.notnull()]
first_lat = float(iso100_coords['lat'][0])
first_lon = float(iso100_coords['lon'][0])

# itere sobre os próximos pontos, adicionando-os ao xarray com um espaçamento de 80 km
lats, lons = [], []
for lon in iso100.lon:
    for lat in iso100.lat:
        #print(float(lat),float(lon))
        if iso100.sel(lat=lat,lon=lon).values == -100:
            lats.append(float(lat)), lons.append(float(lat))

iso=pd.DataFrame([lats,lons]).T
iso = pd.DataFrame({'lat': lats, 'lon': lons})
iso.to_csv('teste_chico/pontos_iso_s_se.csv')

# selecione as coordenadas do ponto atual
list_of_coords = pd.read_csv('/p1-nemo/rtecchio/teste_chico/pontos_iso_s_se.csv')
isobata_100 = list_of_coords.sort_values(by='lat')
isobata_filtrada_1 = isobata_100[~(isobata_100['lat'] > -17) & (isobata_100['lon'] < -40)]
isobata_100_coords = isobata_filtrada_1[~(isobata_filtrada_1['lon'] < -52.5)]
#isobata_100_coords = isobata_filtrada_2[~(isobata_filtrada_2['lat'] > -30) & (isobata_filtrada_2['lon'] < -49.5)]
isobata_100_coords.reset_index(inplace=True)

#coords = pd.DataFrame(columns=['lat','lon'])
# first_lat = float(iso100_coords['lat'][0])
#first_lon = float(iso100_coords['lon'][0])
first_lat = float(isobata_100_coords['lat'][0])
first_lon = float(isobata_100_coords['lon'][0])

final_lat = []
final_lon =[]



for i in range(len(isobata_100_coords)):
    #current_lat = float(lats[i].values)
    #current_lon = float(lons[i].values)
    current_lat = float(isobata_100_coords['lat'][i])
    current_lon = float(isobata_100_coords['lon'][i])
    #coords = pd.concat([coords, pd.DataFrame([[current_lat,current_lon]], columns=['lat','lon'])])
    #print(current_lat,current_lon) 
    
    # calcule a distância entre o primeiro ponto e o ponto atual
    R = 6371  # raio da terra em km
    dLat = math.radians(current_lat - first_lat)
    dLon = math.radians(current_lon - first_lon)
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(math.radians(first_lat)) * math.cos(math.radians(current_lat)) * math.sin(dLon/2) * math.sin(dLon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    distance = R * c
    print (distance)
    # se a distância for maior que 80 km, adicione o ponto atual ao xarray
    # if  80<= distance <=85:
        # print (distance)
        # #ds = iso100.sel(lat=current_lat, lon=current_lon, method="nearest")
        # #ds = iso100.sel[{'lat': current_lat, 'lon': current_lon}]
        # # atualize o primeiro ponto
        # final_lat.append(current_lat)
        # final_lon.append(current_lon)
    if  200 <= distance <= 210:
       # print (distance)
        #ds = iso100.sel(lat=current_lat, lon=current_lon, method="nearest")
        #ds = iso100.sel[{'lat': current_lat, 'lon': current_lon}]
        # atualize o primeiro ponto
        final_lat.append(current_lat)
        final_lon.append(current_lon)
        first_lat = current_lat
        first_lon = current_lon
    if distance > 210:
        first_lat = current_lat
        first_lon = current_lon

#precisei criar um arquivo CSV pois não consegui pontos no ES aí criei um manualmente
ptos_iwmo = pd.DataFrame([final_lat,final_lon]).T
ptos_iwmo = pd.DataFrame({'lat': final_lat, 'lon': final_lon})
ptos_iwmo.to_csv('teste_chico/pontos_iwmo.csv')

ptos = pd.read_csv('/p1-nemo/rtecchio/teste_chico/pontos_iwmo.csv',sep=",")
lat=ptos['lat']
lon=ptos['lon']

#Plota os pontos na iso de 100.
import matplotlib
import matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs 
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, COASTLINE
from cartopy.feature import BORDERS
import cmocean as cmo

google_maps_colors = ['#ebebeb', '#c7c7c7', '#b0b0b0', '#999999', '#7e7e7e', '#666666', '#4d4d4d', '#333333', '#1a1a1a', '#0c7b12']
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
   
# make figure with points 

data_min = prof_filtred.min()
data_max = prof_filtred.max()
interval = 50
levels = np.linspace(data_min,data_max,interval)

proj = ccrs.PlateCarree() 
plt.close('all')
fig = plt.figure(constrained_layout=False,figsize=(12,12))
ax = fig.add_subplot(111, projection=proj)
ax.set_extent([-60,-37,-34,-17.5])
grid_labels_params(ax,0)
Brazil_states(ax)
#cf1 = ax.contourf(iso100.lon,iso100.lat,iso100)
#cf1= ax.contourf(prof_filtred.lon,prof_filtred.lat,prof_filtred, cmap='cmo.deep_r',)
cf1= ax.contourf(prof_filtred.lon,prof_filtred.lat,prof_filtred,levels=levels, cmap='cmo.deep_r',)
ax.plot(lon, lat, 'o',color='k',  )
ax.text(lon[1]+0.2,lat[1]-0.2, '02',  color='black', fontsize=14, transform=proj)
ax.text(lon[3]+0.2,lat[3]-0.2, '04',  color='black', fontsize=14, transform=proj)
ax.text(lon[5]+0.2,lat[5]-0.2, '06',  color='black', fontsize=14, transform=proj)
ax.text(lon[7]+0.2,lat[7]-0.2, '08',  color='black', fontsize=14, transform=proj)
ax.text(lon[9]+0.2,lat[9]-0.2, '10',  color='black', fontsize=14, transform=proj)
# cb=plt.colorbar(cf1, ax=ax, shrink=0.8, aspect=15)
# cb.ax.tick_params(labelsize=12)
# cb.set_label('Profundidade (m)',fontsize=12)
map_features(ax)
#plt.show()
plt.savefig('Pontos_IWMO.png', dpi=500)