%% Script for estimating CDHE frequency and drought/non-drought composite STA and SSI (For Fig.2)
%% Written by Akshay Rajeev
close all; clear all; clc;

ifile2 = '\Data\ERA5_SM_pent.nc';
ifile3 = '\Data\ERA5_Tx_pent.nc';
ifile4 = '\Data\Processed\ERA5_SSI.nc';
ifile5 = '\Data\ERA5_STA.nc';
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

tot_area = 31*31*sum(sum(~isnan(SSI(:,:,1)),2));

thresh = 1.0;

sm(sm==-9999)=NaN;
zz = SSI;
zz(~isnan(SSI))=0; %% Creating new meshgrid as same size as SSI
zz(STA>=thresh & SSI<=(-thresh))=1; %% Estimating CDHEs occurrence at each grid at each timestep

dry = [1960;1965;1968;1972;1974;1979;1982;1986;1987;2002;2014];
yrs=(1950:2020)';
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
mat1 = zeros([25 71]);

%% Estimating All-India frequency of CDHEs
dss=[];dst=[];ss=[];st=[];
for k = 1:length(yrs)
    id = find(tstan(:,k)>=thresh & sstan(:,k)<=-thresh); %% Estimating CDHEs for each year
    for kk = 1:length(id)
    tp = find(t2==ttt(id(kk),k));
    if ismember(yrs(k),dry)
        dss = cat(3,dss,SSI(:,:,tp));
        dst = cat(3,dst,STA(:,:,tp));
    else
        ss = cat(3,ss,SSI(:,:,tp));
        st = cat(3,st,STA(:,:,tp));
    end 
    end
%     area=[];
%     for kk = 1:length(id)
%         tp = find(t2==ttt(id(kk),k));
%         ar=[];
%         ar = ((31*31*nansum(nansum(zz(:,:,tp))))./tot_area)*100;
%        area = [area;ar];
%     end
%      frq(k,1) = sum(area>=20);
     mat1(id,k)=1;
     frq2(k,1) = length(id);
    
end
    norm = yrs(find(frq2==0));
    ndry = setdiff(setdiff(yrs,dry),norm);
    [val1,pos1]=intersect(yrs,dry);
    [val2,pos2]=intersect(yrs,ndry);
    [val3,pos3]=intersect(yrs,norm);
    ff = yrs;
    ff(pos1)=-1;
    ff(pos2)=1;
    ff(pos3)=0;
    %FF = [FF,ff];
    F = [yrs,ff];
s1 = Sen_Slope(frq2);
[h1,p1] = Mann_Kendall(frq2,0.05); %% Trend in CDHE frequency

scomp = mean(ss,3); %% Composite SSI during non-drought
dscomp = mean(dss,3); %% Composite SSI during drought
tcomp = mean(st,3); %% Composite STA during non-drought
dtcomp = mean(dst,3); %% Composite STA during drought

