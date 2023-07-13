# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    surface_gifs.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/13 16:29:58 by Danilo            #+#    #+#              #
#    Updated: 2023/07/13 19:00:41 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from celluloid import Camera

PATH_ERA_FILES = '/p1-nemo/mbonjour/novos_wt_dados/'
DF_FILENAME = '../database/dates-limits.csv'
OUTPUT_DIR = 'animation_output/'

def load_dataframe(filename):
    df = pd.read_csv(filename, index_col=0, parse_dates=[0])
    return df

def load_data(year, start_date, end_date):
    u_file = os.path.join(PATH_ERA_FILES, f'u_{year}_maior.nc')
    v_file = os.path.join(PATH_ERA_FILES, f'v_{year}_maior.nc')
    hgt_file = os.path.join(PATH_ERA_FILES, f'hgt_{year}_maior.nc')

    u_data = xr.open_dataset(u_file).sel(time=slice(start_date, end_date))
    v_data = xr.open_dataset(v_file).sel(time=slice(start_date, end_date))
    hgt_data = xr.open_dataset(hgt_file).sel(time=slice(start_date, end_date))

    return u_data, v_data, hgt_data

def plot_map_hgt_winds(ax, lon, lat, hgt, u, v, date, subsampling_factor=10):
    skip_coords = (slice(None, None, subsampling_factor))
    skip_vars = (slice(None, None, subsampling_factor), slice(None, None, subsampling_factor))

    ax.clear()

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

    # Add a colorbar with appropriate vmin and vmax values
    cbar = plt.colorbar(cf, ax=ax, orientation='horizontal',
                        pad=0.05, label='Geopotential Height', extend='both')

    return cbar

def create_animation(df, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for date_index, row in df.iterrows():
        start_date = row['start']
        end_date = row['end']
        year = date_index.year

        u_data, v_data, hgt_data = load_data(year, start_date, end_date)

        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
        camera = Camera(fig)

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
            cbar = plot_map_hgt_winds(ax, lon, lat, hgt, u, v, date)

            # Capture the current frame
            camera.snap()

        animation = camera.animate()

        # Create animation for the row and save as GIF
        filename = os.path.join(output_dir, f'animation_row_{date_index.strftime("%Y-%m-%d")}.gif')
        animation.save(filename, writer='pillow', fps=4, dpi=150)

        # Close the opened files
        u_data.close()
        v_data.close()
        hgt_data.close()

def main():
    df = load_dataframe(DF_FILENAME)
    create_animation(df, OUTPUT_DIR)

if __name__ == '__main__':
    main()
