%% Script to estimate Monsoon Precip anomaly, Monsoon Temp anomaly, CDHE Frequency, dry spell frequency and warm spell frequency
%% Written by Akshay Rajeev
close all; clear all; clc;

%% Initializing data
ifile1 = '\dir\Precipitation\ERA5_Pr.nc';
ifile2 = '\dir\Soil_Moisture\ERA5_SM.nc';
ifile3 = '\dir\Temperature\ERA5_Tmax.nc';
ifile4 = '\dir\Soil_Moisture\ERA5_SM_pent.nc';
ifile5 = '\dir\Temperature\ERA5_Tx_pent.nc';

lon = ncread(ifile2,'longitude');
lat = ncread(ifile2,'latitude');
pr = ncread(ifile1,'Prec'); %%  Daily Precipitation raw data
sm = ncread(ifile2,'SM'); %% Daily Soil Moisture raw data
tmx = ncread(ifile3,'Tmax'); %% Daily Temperature raw data
smp = ncread(ifile4,'SM'); %% Pentad Soil Moisture
tmp = ncread(ifile5,'Tmax'); %% Pentad Temperature

time = ncread(ifile2,'time');
t = datetime(1950,01,01)+days(time);
[y,m,d] = ymd(t);
time2 = ncread(ifile4,'time');
t2 = datetime(1950,01,03)+days(time2);
[y2,m2,d2] = ymd(t2);

pr(pr==-9999)=NaN;
sm(sm==-9999)=NaN;
tmx(tmx==-9999)=NaN;
smp(smp==-9999)=NaN;
tmp(tmp==-9999)=NaN;
yrs = (1950:2020)';
thresh = 1;

ddir = '\data\ISMR\';
B1 = dir(fullfile(ddir,'*'));
fname = setdiff({B1.name},{'.','..'});

mat1 = zeros([25 71 5]);
FF=[];dd=[];
for ii = 1:length(fname)%% Looping through All-India and each Homogenous rainfall region
    ii
    LL = readtable([ddir,fname{ii}]);
    LON = LL.Field1;
    LAT = LL.Field2;
    p=[];tx=[];s=[];sp=[];tp=[];
    ['Loop 1 Started']
    for i = 1:length(LL.FID)
        id1 = find(lat==LAT(i));
        id2 = find(lon==LON(i));
        p = [p,reshape(pr(id2,id1,:),length(t),1)];
        s = [s,reshape(sm(id2,id1,:),length(t),1)];
        tx = [tx,reshape(tmx(id2,id1,:),length(t),1)];
        sp = [sp,reshape(smp(id2,id1,:),length(t2),1)];
        tp = [tp,reshape(tmp(id2,id1,:),length(t2),1)];
    end
    ['Loop 1 Ended']
	%% Spatial Averaging
    p=mean(p,2);
    s=mean(s,2);
    tx=mean(tx,2);
    sp=mean(sp,2);
    tp=mean(tp,2);
    
    ['Loop 2 Started']
    td=[];Pi=[];Si=[];Ti=[];pi=[];si=[];ti=[];sip=[];tip=[];
    for yy =1950:2020
        tidx = (find(y == yy & m == 6):find(y == yy & m == 9,1,'last'))';
        tidx2 = (find(y2 == yy & m2 == 6):find(y2 == yy & m2 == 9,1,'last'))';
            if length(tidx2)<25
                tidx2 = [tidx2;tidx2(end)+1];
            end
       td = [td;t(tidx)];
       Pi = [Pi; sum(p(tidx))];
       Si = [Si; mean(s(tidx))];
       Ti = [Ti; mean(tx(tidx))];
       pi = [pi, (p(tidx))];
       ti = [ti, (tx(tidx))];
       
       sip = [sip, (sp(tidx2))];
       tip = [tip, tp(tidx2)];
    end
    ['Loop 2 Ended']
    Pclim = mean(Pi);
    Tclim = mean(Ti);
    Sclim = mean(Si);
    Pan = (Pi-Pclim);
    Panp = ((Pi-Pclim)./Pclim)*100; % Precipitation Anomaly during summer monsoon
    Tan = Ti-Tclim; % Temperature Anomaly during summer monsoon
    San = Si-Sclim;
    Ssi = San./std(San); %% Annual monsoon standardized soil moisture anomaly  to estimate drought/non-drought year
    pclim = mean(pi,2);
    tclim = mean(ti,2);
        
        dyrs = yrs(find(Ssi<=-1)); %Estimation of drought years
        ['Loop 2 Started']
        sta=[];spi=[];
        for ij = 1:length(1950:2020)
            tan = ti(:,ij)-tclim;
            sta = [sta;tan/std(tan)];
            pan = pi(:,ij)-pclim;
            spi=[spi;pan/std(pan)];
            
            tpan = tip(:,ij)-mean(tip,2);
            tstan(:,ij)= tpan./std(tpan);
            span = sip(:,ij)-mean(sip,2);
            ssi(:,ij)= span./std(span);           
        end
        ['Loop 3 Ended']
		%%% Estimation of dry, wet and warm spells
        yr = ymd(td);
        dr=[];wt=[];
        for j=1950:2020
            idy = spi(find(yr==j));
            wet=0;dry=0;
            for jj=2:length(idy)-1
                win = idy(jj-1:jj+1);
                if all(win <= -thresh)
                   dry = dry+1;
                elseif all(win >= thresh)
                   wet=wet+1;
                end
            end
            dr=[dr;dry];
            wt=[wt;wet];
        end
        prc = prctile(sta,90);
        ht=[];
        for j2=1950:2020
            idy2 = sta(find(yr==j2));
            hot=0;
            for jj2=2:length(idy2)-1
                win2 = idy2(jj2-1:jj2+1);
                if all(win2 >= prc)
                   hot = hot+1;

                end
            end
            ht=[ht;hot];
        end
        
    %% Estimation of CDHE frequency
    for k = 1:length(1950:2020)
        id = find(tstan(:,k)>=thresh & ssi(:,k)<=-thresh);
%         mat1(id,k,ii)=1;
        frq(k,1) = length(id);
    end
    ['Loop 4 Ended']
	%% Estimation of normal and non-drought years
    norm = yrs(find(frq==0));
    ndry = setdiff(setdiff(yrs,dyrs),norm);
    nn = setdiff(yrs,dyrs);
    [val1,pos1]=intersect(yrs,dyrs);
    [val2,pos2]=intersect(yrs,ndry);
    [val3,pos3]=intersect(yrs,norm);
    ff = yrs;
    ff(pos1)=-1;
    ff(pos2)=1;
    ff(pos3)=0;
    FF = [FF,ff];
    F = [yrs,FF];
	%% Fraction of drought and non-drought CDHEs
    drt = sum(frq(pos1))/(length(dyrs)*25);
    ndrt = sum(frq(pos2))/(length(nn)*25);
    Cint = tstan-ssi; %% CDHE Intensity
    ['Loop 5 Started']
Int = [];
for kk = 1:length(yrs)
    id3 = find(tstan(:,kk)>=thresh & ssi(:,kk)<=-thresh);
    int1 = mean(Cint(id3,kk));
    if isnan(int1)
        int1=0;
    end
    Int=[Int;int1];
end
['Loop 5 Ended']

outdata = [Panp,Tan,frq,dr,ht]; %% Monsoon Precip anomaly, Monsoon Temp anomaly, CDHE Frequency, dry spell frequency, warm spell frequency]
dlmwrite(['\outdir\ISMR_data\',char(fname(ii))],outdata,' ');
    
end
    