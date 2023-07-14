# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    surface_gifs.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo  <danilo.oceano@gmail.com>          +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/13 16:29:58 by Danilo            #+#    #+#              #
#    Updated: 2023/07/14 12:55:44 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PATH_ERA_FILES = '/p1-nemo/mbonjour/novos_wt_dados/'
DF_FILENAME = '../database/dates-limits.csv'
OUTPUT_DIR = '../figures_surface-patterns'

def load_dataframe(filename):
    df = pd.read_csv(filename, index_col=0, parse_dates=[0])
    return df

def load_data(year, start_date, end_date):
    u_file = os.path.join(PATH_ERA_FILES, f'u_{year}_maior.nc')
    v_file = os.path.join(PATH_ERA_FILES, f'v_{year}_maior.nc')
    hgt_file = os.path.join(PATH_ERA_FILES, f'hgt_{year}_maior.nc')

    if not all(os.path.exists(file) for file in [u_file, v_file, hgt_file]):
        print(f"Error: Some files do not exist for year {year}")
        return None

    u_data = xr.open_dataset(u_file).sel(time=slice(start_date, end_date))
    v_data = xr.open_dataset(v_file).sel(time=slice(start_date, end_date))
    hgt_data = xr.open_dataset(hgt_file).sel(time=slice(start_date, end_date))

    return u_data, v_data, hgt_data

def plot_map_hgt_winds(ax, lon, lat, hgt, u, v, date, subsampling_factor=10):
    skip_coords = (slice(None, None, subsampling_factor))
    skip_vars = (slice(None, None, subsampling_factor), slice(None, None, subsampling_factor))

    # Plot geopotential height (hgt)
    cf = ax.contourf(lon, lat, hgt, cmap='coolwarm', levels=20, extend='both')
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

    return cf

def create_panels(df, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for date_index, row in df.iterrows():
        start_date = row['start']
        end_date = row['end']
        year = date_index.year

        data = load_data(year, start_date, end_date)

        if data is None:
            continue

        u_data, v_data, hgt_data = data

        times = pd.to_datetime(u_data.time.values)
        times_12h = times[times.hour % 12 == 0]

        num_timesteps = len(times_12h)
        num_cols = min(num_timesteps, 4)
        num_rows = (num_timesteps + num_cols - 1) // num_cols
        num_subplots = num_rows * num_cols

        print(f'Starting date: {start_date}, Ending date: {end_date}')
        print(f'Number of rows: {num_rows}, number of columns: {num_cols}')
        print(f'Total number of figures: {num_subplots}')
        print(f'Length of time variable: {len(times)}')
        print(f'Number of times subsampled for 12h: {len(times_12h)}')

        fig = plt.figure(figsize=(12, 9))
        gs = GridSpec(num_rows, num_cols, figure=fig)
        cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.03])  # Colorbar axes

        cbar = None  # Initialize colorbar variable

        for t, time_12h in enumerate(times_12h):
            # Extract the variables for the current time step
            u = u_data.sel(time=time_12h)['u']
            v = v_data.sel(time=time_12h)['v']
            hgt = hgt_data.sel(time=time_12h)['z']
            date = time_12h

            # Get longitude and latitude coordinates
            lon = u.longitude
            lat = u.latitude

            # Create subplots within the grid
            print(f'Plotting panel for time: {date}')
            ax = fig.add_subplot(gs[t // num_cols, t % num_cols], projection=ccrs.PlateCarree())

            # Plot map of hgt and wind vectors
            cf = plot_map_hgt_winds(ax, lon, lat, hgt, u, v, date)

            # Normalize colorbar
            if cbar is None:
                cbar = plt.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                                    pad=0.05, label='Geopotential Height', extend='both')
            else:
                cf.norm.vmin = cbar.norm.vmin
                cf.norm.vmax = cbar.norm.vmax

        # Adjust the layout and spacing between subplots
        plt.tight_layout(pad=3)

        # Save the figure
        filename = os.path.join(output_dir, f'panels_row_{date_index.strftime("%Y-%m-%d")}.png')
        plt.savefig(filename, dpi=150)
        print(f'Panel saved: {filename}')

        # Close the opened files
        u_data.close()
        v_data.close()
        hgt_data.close()

def main():
    df = load_dataframe(DF_FILENAME)
    create_panels(df, OUTPUT_DIR)

if __name__ == '__main__':
    main()
