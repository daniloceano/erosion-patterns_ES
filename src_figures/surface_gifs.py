# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    surface_gifs.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/13 16:29:58 by Danilo            #+#    #+#              #
#    Updated: 2023/07/13 17:32:31 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.animation as animation

path_era_files = '/p1-nemo/mbonjour/novos_wt_dados/'
df = pd.read_csv('../database/dates-limits.csv', index_col=0, parse_dates=[0])

def plot_map_hgt_winds(lon, lat, hgt, u, v, date, subsampling_factor=10):
    skip_coords = (slice(None, None, subsampling_factor))
    skip_vars = (slice(None, None, subsampling_factor), slice(None, None, subsampling_factor))

    # Create a new figure and axes with a specific projection
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    # Plot geopotential height (hgt)
    cf = ax.contourf(lon, lat, hgt, cmap='coolwarm')
    ax.coastlines()

    # Plot wind vectors (u and v)
    ax.quiver(lon[skip_coords], lat[skip_coords], u[skip_vars], v[skip_vars],
              transform=ccrs.PlateCarree(), scale=200,
              width=0.003, headwidth=4, headlength=5, headaxislength=6,
              alpha=0.8, zorder=10)

    # Set plot title and labels
    ax.set_title(f'{date:%Y-%m-%d %HZ}')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Add a colorbar
    cbar = plt.colorbar(cf, ax=ax, orientation='horizontal',
                        pad=0.05, label='Geopotential Height')

    return fig, ax

def animate_frames(frames, filename, interval=200):
    ani = animation.ArtistAnimation(fig, frames, interval=interval, blit=True)
    ani.save(filename, writer='pillow', fps=4, dpi=150)

for date_index, row in df.iterrows():
    start_date = row['start']
    end_date = row['end']
    year = date_index.year

    u_file = os.path.join(path_era_files, f'u_{year}_maior.nc')
    v_file = os.path.join(path_era_files, f'v_{year}_maior.nc')
    hgt_file = os.path.join(path_era_files, f'hgt_{year}_maior.nc')

    u_data = xr.open_dataset(u_file).sel(time=slice(start_date, end_date))
    v_data = xr.open_dataset(v_file).sel(time=slice(start_date, end_date))
    hgt_data = xr.open_dataset(hgt_file).sel(time=slice(start_date, end_date))

    frames = []
    for t in range(len(u_data.time)):
        # Extract the variables for the current time step
        u = u_data.isel(time=t)['u']
        v = v_data.isel(time=t)['v']
        hgt = hgt_data.isel(time=t)['z']
        date = pd.to_datetime(u_data.time[t].values)

        # Get longitude and latitude coordinates
        lon = u.longitude
        lat = u.latitude

        # Plot map of hgt and wind vectors
        fig, ax = plot_map_hgt_winds(lon, lat, hgt, u, v, date)
        frames.append([ax])

        plt.close()

    # Create animation for the row and save as GIF
    filename = f'animation_row_{date_index.strftime("%Y-%m-%d")}.gif'
    animate_frames(frames, filename)

    # Close the opened files
    u_data.close()
    v_data.close()
    hgt_data.close()