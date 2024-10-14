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
dsl = 1; % if =1: depthscices, if =0: timeslices
%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
if ~ispc; menu('Select folder with timeslices','OK'); end
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Select folder with timeslices');
        else
            pfad=uigetdir([],'Select folder with timeslices');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Select folder with timeslices'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
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
        if isfield(temp,'followTopo')
            followTopo=temp.followTopo;
        else
            beep
            followTopo=input('Depthslices: Give followTopo=1 or 0: ');
        end
    else
        maxElev=[];
        dsl=0;
        followTopo=0;
    end
    Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,pfad,dsl,maxElev,followTopo,coordtrans);
elseif exist(fullfile(pfad,'tsl.mat'),'file')
    load(fullfile(pfad,'tsl.mat'));
    load(fullfile(pfad,'topo.mat'));
    load(fullfile(pfad,'xgrid.mat'));
    load(fullfile(pfad,'ygrid.mat'));
    if exist(fullfile(pfad,'depth.mat'),'file')
        temp=load(fullfile(pfad,'depth.mat'));
        dsl=1;
        maxElev=temp.maxElevation;
        if isfield(temp,'followTopo')
            followTopo=temp.followTopo;
        else
            followTopo=input('Depthslices: Give followTopo=1 or 0: ');
        end
    else
        maxElev=[];
        dsl=0;
        followTopo=0;
    end
    Tsl_slider_plot(xgrid,ygrid,tsl,topo,t_startende,pfad,dsl,maxElev,followTopo,coordtrans);
end

waitfor(gcf);

% set original path
path(oldpath);