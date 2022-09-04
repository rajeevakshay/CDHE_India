%% Script to estimate the distribution of drought/non-drought CDHE Intensity.
%% Written by Akshay Rajeev
close all; clear all; clc;
%% Intialising dataset
ifile2 = '\Data\ERA5_SM_pent.nc';
ifile3 = '\Data\ERA5_Tx_pent.nc';
lon = ncread(ifile2,'longitude');
lat = ncread(ifile2,'latitude');
sm = ncread(ifile2,'SM');
tmx = ncread(ifile3,'Tmax');
time = ncread(ifile2,'time');
t = datetime(1950,01,03)+days(time);
[y,m,d] = ymd(t);

thresh = 1.0;

sm(sm==-9999)=NaN;
tmx(tmx==-9999)=NaN;

dry = [1960;1965;1968;1972;1974;1979;1982;1986;1987;2002;2014];
yrs=1950:2020;
ndry = setdiff(yrs,dry);
ti = [];si = [];tt=[];ttt=[];
for yy =1950:2020
    
        tidx = find(y == yy & m == 6):find(y == yy & m == 9,1,'last'); %% Selecting only monsoon season data
        if length(tidx)<25
            tidx = [tidx,tidx(end)+1];
        end
        tt = [tt;t(tidx)];
        ttt = [ttt t(tidx)];
        si = [si, reshape(nanmean(nanmean(sm(:,:,tidx))),length(tidx),1)];
        ti = [ti, reshape(nanmean(nanmean(tmx(:,:,tidx))),length(tidx),1)];
   
end
%% Spatially averaging for entire India
sclim = mean(si,2);
tclim = mean(ti,2);
yr = ymd(tt);
sta=[];ssi=[];
for ii = 1:length(1950:2020)
    tan = ti(:,ii)-tclim;
    tstan(:,ii) = tan/std(tan);
    sta = [sta;tan/std(tan)];
    san = si(:,ii)-sclim;
    sstan(:,ii) = san/std(san);
    ssi=[ssi;san/std(san)];
end
Cint = tstan-sstan; %% Estimating CDHE intensity
Int = [];
for k = 1:length(yrs)
    id = find(tstan(:,k)>=thresh & sstan(:,k)<=-thresh);
    int1 = mean(Cint(id,k));
    if isnan(int1)
        int1=0;
    end
    Int=[Int;int1];
end
%% Estmating CDHE intensity in drought and non-drought periods.
[val1,pos1]=intersect(yrs,dry);
[val2,pos2]=intersect(yrs,ndry);
dint=Int;
nint = Int;
nint(pos1)=0;
dint(pos2)=0;
di = dint(dint~=0);
ni = nint(nint~=0);

figure;
histogram(di)
hold on
histogram(ni)
figure;
ksdensity(di)
hold on; 
ksdensity(ni)
[fn,xn] = ksdensity(ni);
[fd,xd] = ksdensity(di);
mn=[mean(di),mean(ni)];
[h1,p1] = kstest2(di,ni,0.05);