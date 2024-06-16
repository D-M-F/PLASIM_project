#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# In[2]:


# nice figures
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.dpi']= 150
mpl.rcParams['font.size'] = 12


# In[6]:


# dim0=time, dim1=lat, dim2=lon
ds = xr.open_dataset('prt_1971_2000_ERA5_monthly_T42.nc')
# unit conversion (read the ERA5 dataset description for details)
prt_ERA5 = ds.var228 / (3600*24)

ds1 = xr.open_dataset('prc_1971_2000_PLASIM_monthly_T42.nc')
ds2 = xr.open_dataset('prl_1971_2000_PLASIM_monthly_T42.nc')
# by definition, prt = prc+prl (read the ERA5 dataset description for details)
prt_PLASIM = ds1.prc+ds2.prl


# In[10]:


# global template
def mapplot(var,tt,clabel,res='T42'):
    fig = plt.figure(figsize=(15,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(lw=0.5)
    
    # x_extent = [-180, -120, -60, 0, 60, 120, 180]
    # y_extent = [-90,-60, -30, 0, 30, 60,90]
    lonmin = -180
    lonmax = 180
    latmin = -90
    latmax = 90
    x_extent = np.linspace(lonmin,lonmax,7)
    y_extent = np.linspace(latmin,latmax,7)
    ax.set_xticks(x_extent,crs=ccrs.PlateCarree())
    ax.set_yticks(y_extent,crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    # T42 (T21 to be updated)
    lonnum = 128
    latnum = 64
    lon = np.arange(-180,180.1,360/lonnum)
    lat = np.array([ 87.863799,  85.096527,  82.312913,  79.525607,  76.7369  ,  73.947515,
        71.157752,  68.367756,  65.577607,  62.787352,  59.99702 ,  57.206632,
        54.4162  ,  51.625734,  48.835241,  46.044727,  43.254195,  40.463648,
        37.67309 ,  34.882521,  32.091944,  29.30136 ,  26.510769,  23.720174,
        20.929574,  18.138971,  15.348365,  12.557756,   9.767146,   6.976534,
         4.185921,   1.395307,  -1.395307,  -4.185921,  -6.976534,  -9.767146,
       -12.557756, -15.348365, -18.138971, -20.929574, -23.720174, -26.510769,
       -29.30136 , -32.091944, -34.882521, -37.67309 , -40.463648, -43.254195,
       -46.044727, -48.835241, -51.625734, -54.4162  , -57.206632, -59.99702 ,
       -62.787352, -65.577607, -68.367756, -71.157752, -73.947515, -76.7369  ,
       -79.525607, -82.312913, -85.096527, -87.863799])
    # Convert lat/lon to map coordinates
    Lon, Lat = np.meshgrid(lon, lat)

    # Plot the data
    var = np.roll(np.mean(var,0),int(lonnum/2),1)
    # fill the blank 180E (only useful for global plots)
    np.append(var,var[:,0].reshape(latnum,1),axis=1)
    pm = ax.pcolormesh(Lon, Lat, var, cmap=plt.cm.jet, vmin=0, vmax=2e-7)

    # Add a colorbar
    cb = plt.colorbar(pm, label=r"%s"%clabel, orientation="vertical")

    ax.set_title(r"%s"%tt)
    ax.set_aspect('auto')
    # plt.savefig("figs/%s"%tt, bbox_inches='tight')
    plt.show()

mapplot(prt_ERA5,'mean total precipitation in ERA5 (1971-2000)','total precipitation [m/s]')
mapplot(prt_PLASIM,'mean total precipitation in PLASIM (1971-2000)','total precipitation [m/s]')


# In[20]:


# regional template
def mapplot(var,tt,clabel,res='T42'):
    fig = plt.figure(figsize=(15,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(lw=0.5)
    
    # x_extent = [-180, -120, -60, 0, 60, 120, 180]
    # y_extent = [-90,-60, -30, 0, 30, 60,90]
    lonmin = -10
    lonmax = 22
    latmin = 35
    latmax = 48
    x_extent = np.linspace(lonmin,lonmax,7)
    y_extent = np.linspace(latmin,latmax,7)
    ax.set_xticks(x_extent,crs=ccrs.PlateCarree())
    ax.set_yticks(y_extent,crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    # T42 (T21 to be updated)
    lonnum = 128
    latnum = 64
    lon = np.arange(-180,180.1,360/lonnum)
    lat = np.array([ 87.863799,  85.096527,  82.312913,  79.525607,  76.7369  ,  73.947515,
        71.157752,  68.367756,  65.577607,  62.787352,  59.99702 ,  57.206632,
        54.4162  ,  51.625734,  48.835241,  46.044727,  43.254195,  40.463648,
        37.67309 ,  34.882521,  32.091944,  29.30136 ,  26.510769,  23.720174,
        20.929574,  18.138971,  15.348365,  12.557756,   9.767146,   6.976534,
         4.185921,   1.395307,  -1.395307,  -4.185921,  -6.976534,  -9.767146,
       -12.557756, -15.348365, -18.138971, -20.929574, -23.720174, -26.510769,
       -29.30136 , -32.091944, -34.882521, -37.67309 , -40.463648, -43.254195,
       -46.044727, -48.835241, -51.625734, -54.4162  , -57.206632, -59.99702 ,
       -62.787352, -65.577607, -68.367756, -71.157752, -73.947515, -76.7369  ,
       -79.525607, -82.312913, -85.096527, -87.863799])
    lonidx = np.where((lon>=lonmin) & (lon<=lonmax))[0]
    latidx = np.where((lat>=latmin) & (lat<=latmax))[0]
    lon = lon[lonidx]
    lat = lat[latidx]
    # Convert lat/lon to map coordinates
    Lon, Lat = np.meshgrid(lon, lat)

    # Plot the data
    var = np.roll(np.mean(var,0),int(lonnum/2),1)
    # # fill the blank 180E (only useful for global plots)
    # np.append(var,var[:,0].reshape(latnum,1),axis=1)
    var = var[latidx,:][:,lonidx]
    pm = ax.pcolormesh(Lon, Lat, var, cmap=plt.cm.jet, vmin=0, vmax=2e-7)

    # Add a colorbar
    cb = plt.colorbar(pm, label=r"%s"%clabel, orientation="vertical")
    
    plt.xlim(lonmin,lonmax)
    plt.ylim(latmin,latmax)

    ax.set_title(r"%s"%tt)
    ax.set_aspect('auto')
    # plt.savefig("figs/%s"%tt, bbox_inches='tight')
    plt.show()

mapplot(prt_ERA5,'mean total precipitation in ERA5 (1971-2000)','total precipitation [m/s]')
mapplot(prt_PLASIM,'mean total precipitation in PLASIM (1971-2000)','total precipitation [m/s]')

