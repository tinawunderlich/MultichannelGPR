clear all
close all
clc

% Use make_Timeslices.m first and then use this script for plotting! Tsl
% only! (Select folder with timeslices)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires Tsl_slider_plot in folder Plotting and Subfunctions


% depth slices
dsl = 0;
%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!

if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Select folder with timeslices');
        else
            pfad=uigetdir([],'Select folder with timeslices');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Select folder with timeslices'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Select folder with timeslices');
        else
            pfad=uigetdir([],'Select folder with timeslices');
        end
    else
        pfad=uigetdir([],'Select folder with timeslices'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'),fullfile(curFold,'Subfunctions'));


load(fullfile(pfad,'t_startende.mat'));
if exist(fullfile(pfad,'coordtrans.mat'),'file')
    load(fullfile(pfad,'coordtrans.mat'));
else
    coordtrans=[1 1 1 1; 2 2 2 2];
end
if exist(fullfile(pfad,'tsl_interp.mat'),'file')
    load(fullfile(pfad,'tsl_interp.mat'));
    load(fullfile(pfad,'topo_interp.mat'));
    load(fullfile(pfad,'xgrid_interp.mat'));
    load(fullfile(pfad,'ygrid_interp.mat'));
    if exist(fullfile(pfad,'depth.mat'),'file')
        temp=load(fullfile(pfad,'depth.mat'));
        dsl=1;
        maxElev=temp.maxElevation;
    else
        maxElev=[];
        dsl=0;
    end
    Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,pfad,dsl,maxElev,coordtrans);
elseif exist(fullfile(pfad,'tsl.mat'),'file')
    load(fullfile(pfad,'tsl.mat'));
    load(fullfile(pfad,'topo.mat'));
    load(fullfile(pfad,'xgrid.mat'));
    load(fullfile(pfad,'ygrid.mat'));
    if exist(fullfile(pfad,'depth.mat'),'file')
        temp=load(fullfile(pfad,'depth.mat'));
        dsl=1;
        maxElev=temp.maxElevation;
    else
        maxElev=[];
        dsl=0;
    end
    Tsl_slider_plot(xgrid,ygrid,tsl,topo,t_startende,pfad,dsl,maxElev,coordtrans);
end

waitfor(gcf);

% set original path
path(oldpath);