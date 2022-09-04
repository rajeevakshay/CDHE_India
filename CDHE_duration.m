%% Script to estimate the distribution of drought/non-drought CDHE duration.
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
%% Estimation of CDHE durations
dur=[];
for j=1950:2020
    idy = find(yr==j);
id = find(sta(idy)>=thresh & ssi(idy)<=-thresh);
if length(id)==1
    dur = [dur;1];
elseif length(id)>1
    days=1;ds=[];
    for jj=2:length(id)
        df = id(jj)-id(jj-1);
        if df == 1
            days = days+1;
        elseif df>1
            ds = [ds;days];
            days=1;            
            
        end
        
        
    end
    ds = [ds;days];
    dur=[dur;mean(ds)];
else
    dur = [dur;0];
end
end
%% Estmating CDHE duration in drought and non-drought periods.
[val1,pos1]=intersect(yrs,dry);
[val2,pos2]=intersect(yrs,ndry);
ddur=dur;
ndur = dur;
ndur(pos1)=0;
ddur(pos2)=0;
drd = ddur(ddur~=0);
nrd = ndur(ndur~=0);


figure;
ksdensity(drd)
hold on; 
ksdensity(nrd)

[fn,xn] = ksdensity(nrd);
[fd,xd] = ksdensity(drd);
mn=[mean(drd),mean(nrd)];
[h1,p1] = kstest2(drd,nrd,0.05);
s1 = Sen_Slope(dur);
%[h1,p1] = Mann_Kendall(dur,0.05);