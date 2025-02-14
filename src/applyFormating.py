import netCDF4 as nc
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy.stats import linregress
import matplotlib
matplotlib.use("Qt5Agg")
import glob
import warnings
import sys
from datetime import datetime
import socket
# -----------------------------
# 0. Set the work directory
#      and file path.
# -----------------------------
# work directory
work_dir = r'/home/lzh/ework/data_ZLa'
os.chdir(work_dir)

# input and output
dat_in = 'dat_in'
dat_out_csv = 'dat_out_csv'
dat_out_nc = 'dat_out_nc'
if not os.path.exists(dat_in):
       warnings.warn("Input folder no found!", RuntimeWarning)
       sys.exit(1)
if not os.path.exists(dat_out_csv):
       os.mkdir(dat_out_csv)
if not os.path.exists(dat_out_nc):
       os.mkdir(dat_out_nc)

# info for files.
files = glob.glob(os.path.join(dat_in, '*.DAT'))
dates = pd.date_range('2020-09-01', '2020-09-30')       # modify based on the dat_in

# -----------------------------
# 1. search for the file/files
# ------------------------------
for date in dates:
       print('Start Processing... ', date)
       # format to 'yyyy_mm_dd'
       date_formated = date.strftime('%Y_%m_%d')

       #  Multiple files may exist per day due to short period failure of the equipments,
       #      reconnection produce new files.
       matching_files = [file for file in files if date_formated in os.path.basename(file)]

       # skip date without data.
       if len(matching_files) == 0:
           continue

       # -----------------------------
       # 2. format data as DataFrame
       #        and save as *.csv
       # ------------------------------

       # site based.
       columns = ['RECORD', 'Ux_B4', 'Uy_B4', 'Uz_B4', 'Ts_B4',
                  'diag_sonic_B4', 'Ux_B3', 'Uy_B3', 'Uz_B3', 'Ts_B3', 'diag_sonic_B3',
                  'Ux_B2', 'Uy_B2', 'Uz_B2', 'Ts_B2', 'diag_sonic_B2', 'Ux_B1', 'Uy_B1',
                  'Uz_B1', 'Ts_B1', 'diag_sonic_B1', 'co2_B4', 'h2o_B4', 'sig_irga_B4',
                  'Press_irga_B4', 'co2_B2', 'h2o_B2', 'sig_irga_B2', 'Press_irga_B2',
                  'FW_B2', 'Tair_B4', 'RH_B4', 'Tair_B2', 'RH_B2']

       # Create formated dataframe, filled with nan.
       steps = 24 * 3600 * 10 # 24 h * 3600 s /h * 10Hz
       day_time = pd.date_range(start=date, periods=steps+1, freq='0.1s', inclusive='right')
       df_day = pd.DataFrame(np.nan, dtype='object', index=day_time, columns=columns)

       # Use df_day.loc[index,:] = df to update the formatted dataframe with real data.
       for idx, file in enumerate(matching_files):
              print(f'Process... {date_formated}', f'{len(matching_files)}: {idx+1}')
              df = pd.read_csv(file, delimiter=',', skiprows=[0, 2, 3], header=0, index_col=0, dtype='object')
              df.index = pd.to_datetime(df.index, format='mixed')

              # df_hour.index → Selects rows in df_day that match df_hour.
              df_day.loc[df.index, :] = df
              del df

       # save as csv
       file_csv = os.path.join(dat_out_csv, f'{date_formated}.csv')
       df_day.to_csv(file_csv)

       # ------------------------------
       # 3. convert data from csv to
       #        netcdf.
       # ------------------------------
       file_nc = os.path.join(dat_out_nc, f'{date_formated}.nc')
       if os.path.exists(file_nc):
           os.remove(file_nc)
       ncfile = nc.Dataset(file_nc, mode='w', format='NETCDF4')

       # Create dimensions
       time_dim = ncfile.createDimension('time', steps)
       height_dim = ncfile.createDimension('height', 4)

       # Create variables
       Ux = ncfile.createVariable('Ux', np.float32, ('time', 'height'))
       Ux.units = 'm/s'
       Ux.long_name = 'Zonal wind'

       Uy = ncfile.createVariable('Uy', np.float32, ('time', 'height'))
       Uy.units = 'm/s'
       Uy.long_name = 'Meridional wind'

       Uz = ncfile.createVariable('Uz', np.float32, ('time', 'height'))
       Uz.units = 'm/s'
       Uz.long_name = 'Vertical wind'

       Ts = ncfile.createVariable('Ts', np.float32, ('time', 'height'))
       Ts.units = 'deg C'
       Ts.long_name = 'Air temperature'

       time_var = ncfile.createVariable("time", "f8", ("time",))
       time_var.units = f"seconds since {date}"
       time_var.calendar = "standard"

       # Fill the data
       time_var[:] = nc.date2num(df_day.index.to_pydatetime(), units=time_var.units, calendar=time_var.calendar)

       Ux[:, 0] = df_day['Ux_B4']
       Ux[:, 1] = df_day['Ux_B3']
       Ux[:, 2] = df_day['Ux_B2']
       Ux[:, 3] = df_day['Ux_B1']
       
       Uy[:, 0] = df_day['Uy_B4']
       Uy[:, 1] = df_day['Uy_B3']
       Uy[:, 2] = df_day['Uy_B2']
       Uy[:, 3] = df_day['Uy_B1']
       
       Uz[:, 0] = df_day['Uz_B4']
       Uz[:, 1] = df_day['Uz_B3']
       Uz[:, 2] = df_day['Uz_B2']
       Uz[:, 3] = df_day['Uz_B1']
        
       Ts[:, 0] = df_day['Ts_B4']
       Ts[:, 1] = df_day['Ts_B3']
       Ts[:, 2] = df_day['Ts_B2']
       Ts[:, 3] = df_day['Ts_B1']

        # add global attributes
       ncfile.setncattr('creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
       ncfile.setncattr('created_by', socket.gethostname())
       ncfile.setncattr('description', '10Hz data from EC, contact: zhiheng.lan@wsu.edu')
       ncfile.close()


       del df_day
