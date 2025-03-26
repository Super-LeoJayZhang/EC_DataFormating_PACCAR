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
dates = pd.date_range('2020-09-01', '2020-09-30')

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

    # Create formated dataframe, filled with nan.
    steps = 24 * 3600 * 10  # 24 h * 3600 s /h * 10Hz
    day_time = pd.date_range(start=date, periods=steps + 1, freq='0.1s', inclusive='right')

    # Use df_day.loc[index,:] = df to update the formatted dataframe with real data.
    for idx, file in enumerate(matching_files):
        # print(idx)
        if idx == 0:
            print(f'Process... {date_formated}', f'{len(matching_files)}: {idx + 1}')
            df = pd.read_csv(file, delimiter=',', skiprows=[0, 2, 3], header=0, index_col=0, dtype='object')
            df.index = pd.to_datetime(df.index, format='mixed')
            df_day = pd.DataFrame(np.nan, dtype='object', index=day_time, columns=df.columns)
            # df_hour.index → Selects rows in df_day that match df_hour.
            try:
                df_day.loc[df.index, :] = df
            except KeyError:
                pass

            del df
        else:
            print(f'Process... {date_formated}', f'{len(matching_files)}: {idx + 1}')
            df = pd.read_csv(file, delimiter=',', skiprows=[0, 2, 3], header=0, index_col=0, dtype='object')
            df.index = pd.to_datetime(df.index, format='mixed')
            print(max(df.index))
            # check whether columns are same
            if len(df_day.columns) == len(df.columns):
                try:
                    df_day.loc[df.index, :] = df
                except KeyError:
                    pass
            else:
                print('measurements changed')
            del df

    # ------------------------------
    # 3. convert data from csv to
    #        netcdf.
    # ------------------------------
    file_nc = os.path.join(dat_out_nc, f'{date_formated}.nc')
    if os.path.exists(file_nc):
        os.remove(file_nc)
    ncfile = nc.Dataset(file_nc, mode='w', format='NETCDF4')

    levels = np.array([9, 10, 11, 12])

    # Create dimensions
    time_dim = ncfile.createDimension('time', steps)
    height_dim = ncfile.createDimension('height', 4)
    base_time_dim = ncfile.createDimension('base_time', 1)
    level_dim = ncfile.createDimension('level', len(levels))

    # Create variables
    level_var = ncfile.createVariable('level', 'f8', ('level', ))
    level_var[:] = levels
    level_var.units = 'index'
    level_var.long_name = 'vertical layer index'
    level_var.positive = 'up'
    level_var.axis = 'Z'

    base_time = ncfile.createVariable('base_time', 'i8')
    base_time.unit = 'second since 1970-1-1 0:00:00'
    base_time.ancillary_variables = 'time_offset'
    base_time.string = f'{date}'
    base_time[:] = int(pd.Timestamp(date).timestamp())


    time_offset = ncfile.createVariable('time_offset', 'f8', ('time'))
    time_offset.long_name = 'Time offset from base_time'
    time_offset.unit = f'seconds since {date}'
    time_offset[:] = np.arange(0, 86400, 0.1)

    time_var = ncfile.createVariable("time", "f8", ("time",))
    time_var.units = f"seconds since {date}"
    time_var.calendar = "standard"
    time_var[:] = np.arange(0, 86400, 0.1)

    Ux = ncfile.createVariable('Ux', np.float32, ('time', 'level'))
    Ux.units = 'm/s'
    Ux.long_name = 'Zonal wind'

    Uy = ncfile.createVariable('Uy', np.float32, ('time', 'level'))
    Uy.units = 'm/s'
    Uy.long_name = 'Meridional wind'

    Uz = ncfile.createVariable('Uz', np.float32, ('time', 'level'))
    Uz.units = 'm/s'
    Uz.long_name = 'Vertical wind'

    Ts = ncfile.createVariable('Ts', np.float32, ('time', 'level'))
    Ts.units = 'deg C'
    Ts.long_name = 'Sonic temperature'

    diag_sonic = ncfile.createVariable('diag_sonic', 'f4', ('time', 'level'))
    diag_sonic.units = '#'
    diag_sonic.long_name = 'Sonic diagnostic'

    h2o = ncfile.createVariable('h2o', np.float32, ('time', 'level'))
    h2o.units = 'g/m^3'
    h2o.long_name = 'water vapor density'

    co2 = ncfile.createVariable('co2', np.float32, ('time', 'level'))
    co2.units = 'mg/m^3'
    co2.long_name = 'CO2 density'

    Press = ncfile.createVariable('Press', np.float32, ('time', 'level'))
    Press.units = 'kPa'
    Press.long_name = 'Air Pressure'

    sig_irga = ncfile.createVariable('sig_irga', 'f4', ('time', 'level'))
    sig_irga.units = '#'
    sig_irga.long_name = 'IRGA signal strength'

    FW = ncfile.createVariable('FW', np.float32, ('time', 'level'))
    FW.units = 'degree C'
    FW.long_name = 'FW05 temperature'

    Tair = ncfile.createVariable('Tair', np.float32, ('time', 'level'))
    Tair.units = 'degree C'
    Tair.long_name = '1 hz air temperature'

    RH = ncfile.createVariable('RH', np.float32, ('time', 'level'))
    RH.units = '%'
    RH.long_name = '1 hz relative humidity'

    height_var = ncfile.createVariable("height", "f8", ("level",))
    height_var.units = 'm'
    height_var.long_name = 'layer height'
    height_var[:] = np.array([30.3, 40.2, 50.6, 60.5])

    sonic_azimuth = ncfile.createVariable("sonic_azimuth", "f8", ("level",))
    sonic_azimuth.units = 'degree'
    sonic_azimuth.long_name = 'Sonic azimuth'
    sonic_azimuth[:] = np.array([338.8, 337.2, 337.5, 340.3])

    # Fill the data
    for column in df_day.columns:
        for num in ['1', '2', '3', '4']:
            if num in column:
                column_num = int(num) - 1
                if 'Ux' in column:
                    Ux[:, column_num] = df_day[column]
                if 'Uy' in column:
                    Uy[:, column_num] = df_day[column]
                if 'Uz' in column:
                    Uz[:, column_num] = df_day[column]
                if 'Ts' in column:
                    Ts[:, column_num] = df_day[column]
                if 'diag_sonic' in column:
                    diag_sonic[:, column_num] = df_day[column]
                if 'h2o' in column:
                    h2o[:, column_num] = df_day[column]
                if 'co2' in column:
                    co2[:, column_num] = df_day[column]
                if 'Press' in column:
                    Press[:, column_num] = df_day[column]
                if 'sig_irga' in column:
                    sig_irga[:, column_num] = df_day[column]
                if 'FW' in column:
                    FW[:, column_num] = df_day[column]
                if 'Tair' in column:
                    Tair[:, column_num] = df_day[column]
                if 'RH' in column:
                    RH[:, column_num] = df_day[column]

    # add global attributes
    ncfile.setncattr('creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    ncfile.setncattr('created_by', socket.gethostname())
    ncfile.setncattr('description', 'multi-layer 10Hz data from EC (group A), contact: zhiheng.lan@wsu.edu')
    ncfile.setncattr('layer height (m)', '30.3 40.2 50.6 60.5')
    ncfile.setncattr('Sonic azimuth (degree)', '338.8 337.2 337.5 340.3')
    ncfile.close()




