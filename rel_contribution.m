close all; clear all;clc;
%% Script for estimating relative contribution of each driver
%% Written by Akshay Rajeev
data = readtable('\data\outdominance.txt');

tot = data.Pan+data.Tanm+data.Dsp+data.Wsp;
contP = ((data.Pan)./tot)*100;
contT = ((data.Tanm)./tot)*100;
contDs = ((data.Dsp)./tot)*100;
contWs = ((data.Wsp)./tot)*100;

out = [contP,contT,contDs,contWs]; %%  Relative contributions of Monsoon Precip anomaly, Monsoon Temp anomaly, dry spell frequency and warm spell frequency

