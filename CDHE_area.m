%% Script to estimate the distribution of area (%) affected by CDHEs in India.
%% Written by Akshay Rajeev
close all; clear all; clc;
%% Intialising dataset
ifile2 = '\data\ERA5_SM_pent.nc';
ifile3 = '\data\ERA5_Tx_pent.nc';
ifile4 = '\data\ERA5_SSI.nc';
ifile5 = '\data\ERA5_STA.nc';
lon = ncread(ifile2,'longitude');
lat = ncread(ifile2,'latitude');
sm = ncread(ifile2,'SM');
tmx = ncread(ifile3,'Tmax');
SSI = ncread(ifile4,'SM');
STA = ncread(ifile5,'Tmax');
time = ncread(ifile2,'time');
t = datetime(1950,01,03)+days(time);
[y,m,d] = ymd(t);
time2 = ncread(ifile4,'time');
t2 = datetime(1950,01,03)+days(time2);
[y2,m2,d2] = ymd(t2);

tot_area = 31*31*sum(sum(~isnan(SSI(:,:,1)),2)); %% Estimating total area

thresh = 1.0;

sm(sm==-9999)=NaN;
tmx(tmx==-9999)=NaN;
SSI(SSI==-9999)=NaN;
STA(STA==-9999)=NaN;
zz = SSI;
zz(~isnan(SSI))=0; %% Creating new meshgrid as same size as SSI
zz(STA>=thresh & SSI<=(-thresh))=1; %% Estimating CDHEs occurrence at each grid at each timestep

dry = [1960;1965;1968;1972;1974;1979;1982;1986;1987;2002;2014];
yrs=(1950:2020)';
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
sta=[];ssi=[];
for ii = 1:length(1950:2020)
    tan = ti(:,ii)-tclim;
    tstan(:,ii) = tan/std(tan);
    sta = [sta;tan/std(tan)];
    san = si(:,ii)-sclim;
    sstan(:,ii) = san/std(san);
    ssi=[ssi;san/std(san)];
end

Ar = [];
for k = 1:length(yrs)
    id = find(tstan(:,k)>=thresh & sstan(:,k)<=-thresh); %% Estimating CDHE frequency based on All-India averaged values

    area=[];
    for kk = 1:length(id)
        tp = find(t2==ttt(id(kk),k));
        ar=[];
        ar = ((31*31*nansum(nansum(zz(:,:,tp))))./tot_area)*100; %% Estimating percentage area for each CDHE event
        area = [area;ar];
    end
    Ar = [Ar;mean(area)];
    
end
Ar(isnan(Ar))=0;

%% Estmating CDHE area in drought and non-drought periods.
[val1,pos1]=intersect(yrs,dry);
[val2,pos2]=intersect(yrs,ndry);
darea=Ar;
narea= Ar;
narea(pos1)=0;
darea(pos2)=0;
dar = darea(darea~=0);
nar = narea(narea~=0);


figure;
histogram(dar)
hold on
histogram(nar)
figure;
ksdensity(dar)
hold on;
ksdensity(nar)
xlim([0 50])
[fn,xn] = ksdensity(nar);
[fd,xd] = ksdensity(dar);
[h1,p1] = kstest2(dar,nar,0.05);
mn=[mean(dar),mean(nar)];
