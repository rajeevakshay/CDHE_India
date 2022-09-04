%% Script for identifying drought years based on monsoon season SSI
%% Written by Akshay Rajeev
close all; clear all; clc;
%% Initializing data
ifile1 = '\Data\ERA5_Pr.nc';
ifile2 = '\Data\ERA5_SM.nc';
ifile3 = '\Data\ERA5_Tmax.nc';
lon = ncread(ifile2,'longitude');
lat = ncread(ifile2,'latitude');
pr = ncread(ifile1,'Prec');
sm = ncread(ifile2,'SM');
tmx = ncread(ifile3,'Tmax');
time = ncread(ifile2,'time');
t = datetime(1950,01,01)+days(time);
[y,m,d] = ymd(t);


pr(pr==-9999)=NaN;
sm(sm==-9999)=NaN;
tmx(tmx==-9999)=NaN;

ti = [];si = [];pi=[];
for yy =1950:2020
    
        tidx = find(y == yy & m == 6):find(y == yy & m == 9,1,'last'); %% Selecting only monsoon season data

        pi = [pi;yy nanmean(nanmean(sum(pr(:,:,tidx),3)))];
        si = [si;yy nanmean(nanmean(mean(sm(:,:,tidx),3)))];
        ti = [ti;yy nanmean(nanmean(mean(tmx(:,:,tidx),3)))];
   
end
%% Estimating whole india anomaly of monsoon precipitation, temperature and standardised anomaly of monsoon soil moisture (Also used in Figure.1)
pan = pi(:,2)-mean(pi(:,2));
panp = (pan/mean(pi(:,2)))*100;
tan = ti(:,2)-mean(ti(:,2));
san = si(:,2)-mean(si(:,2));
sstan = san/std(san);

b1 = Sen_Slope(pan);
b2 = Sen_Slope(tan);
b3 = Sen_Slope(sstan);

[h1,p1] = Mann_Kendall(pan,0.05);
[h2,p2] = Mann_Kendall(tan,0.05);
[h3, p3] = Mann_Kendall(sstan,0.05);

outdata=[panp,tan,sstan];

dry = ti(find(sstan<=-1),1); %% Estimation of drougth years based on SSI


