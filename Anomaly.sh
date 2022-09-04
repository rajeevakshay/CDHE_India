#!/bin/bash
## Script for converting daily Tmax and Soil moisture to pentad and to estimate their standardised anomalies
## Written by Akshay Rajeev
ddir='/data/';
cdo timselmean,5 -del29feb $ddir/ERA5_Tmax.nc $ddir/Tx_pent.nc;
cdo timselmean,5 -del29feb $ddir/ERA5_SM.nc $ddir/SM_pent.nc;

cdo ydaymean $ddir/Tx_pent.nc $ddir/Tx_clim.nc;
cdo ydaymean $ddir/SM_pent.nc $ddir/SM_clim.nc;

cdo sub $ddir/Tx_pent.nc $ddir/Tx_clim.nc $ddir/Tx_an.nc;
cdo sub $ddir/SM_pent.nc $ddir/SM_clim.nc $ddir/SM_an.nc;

cdo ydaystd $ddir/Tx_an.nc $ddir/Tx_std.nc;
cdo ydaystd $ddir/SM_an.nc $ddir/SM_std.nc;

cdo div $ddir/Tx_an.nc $ddir/Tx_std.nc $ddir/STA.nc;
cdo div $ddir/SM_an.nc $ddir/SM_std.nc $ddir/SSI.nc;