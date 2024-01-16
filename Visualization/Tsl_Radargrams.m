clear all
close all
clc


% Display Mala-GPR-data as Timeslices and Radargrams with connected plots
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Requires Timeslices from script make_Timeslices (choose either Timeslices-folder 
% or interpolated timeslices) and processed radargrams
% in another folder (no processing of radargrams is provided here, so use
% processed radargrams). Radargrams could be from script make_Radargrams,
% Bins2radargrams, Mala3D or own creation (folder should contain radargrams.mat
% (cells with radargrams), global_coords.mat (cells with starting end
% ending coordinates of profile (for straight profiles) or coordinates for 
% each trace(for curved profiles)), t.mat (time vector)).
%
% You have to give the locations of both folders!
%
% requires content of Plotting and Subfunctions folder

timedepth=1; % 1=time, 2=depth

% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name - TSL
if ~ispc; menu('Choose folder with timeslices','OK'); end
if ispc
    if exist('tsltemp.temp') % read last opened folder from temp.temp
        fid=fopen('tsltemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_tsl=uigetdir(fn{1}{1},'Choose folder with timeslices');
        else
            pfad_tsl=uigetdir([],'Choose folder with timeslices');
        end
        fileattrib('tsltemp.temp','-h');
        fid=fopen('tsltemp.temp','wt');
        fprintf(fid,'%s',pfad_tsl);
        fclose(fid);
        fileattrib('tsltemp.temp','+h');
    else
        pfad_tsl=uigetdir([],'Choose folder with timeslices'); % path to radargram-folder

        fid=fopen('tsltemp.temp','wt');
        fprintf(fid,'%s',pfad_tsl);
        fclose(fid);
        fileattrib('tsltemp.temp','+h');
    end
else
    if exist('.tsltemp.temp') % read last opened folder from temp.temp
        fid=fopen('.tsltemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_tsl=uigetdir(fn{1}{1},'Choose folder with timeslices');
        else
            pfad_tsl=uigetdir([],'Choose folder with timeslices');
        end
    else
        pfad_tsl=uigetdir([],'Choose folder with timeslices'); % path to radargram-folder
    end

    fid=fopen('.tsltemp.temp','wt');
    fprintf(fid,'%s',pfad_tsl);
    fclose(fid);
end


% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with radargrams','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
        fileattrib('radtemp.temp','-h');
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    end
else
    if exist('.radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('.radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end

    fid=fopen('.radtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'),fullfile(curFold,'Subfunctions'));


%%% TSL
% read times for Tsl
t=load(fullfile(pfad_tsl,'t_startende.mat'));
t=t.t_startende;
% read coordinates and Tsl
if exist(fullfile(pfad_tsl,'xgrid_interp.mat'))
    temp=load(fullfile(pfad_tsl,'xgrid_interp.mat'));
    x=temp.xgrid_interp;
else
    temp=load(fullfile(pfad_tsl,'xgrid.mat'));
    x=temp.xgrid;
end
if exist(fullfile(pfad_tsl,'ygrid_interp.mat'))
    temp=load(fullfile(pfad_tsl,'ygrid_interp.mat'));
    y=temp.ygrid_interp;
else
    temp=load(fullfile(pfad_tsl,'ygrid.mat'));
    y=temp.ygrid;
end
if exist(fullfile(pfad_tsl,'mask_interp.mat'))
    temp=load(fullfile(pfad_tsl,'mask_interp.mat'));
    mask=temp.mask_interp;
else
    temp=load(fullfile(pfad_tsl,'mask.mat'));
    mask=temp.mask;
end
if exist(fullfile(pfad_tsl,'tsl_interp.mat'))
    temp=load(fullfile(pfad_tsl,'tsl_interp.mat'));
    tsl=temp.tsl_interp;
else
    temp=load(fullfile(pfad_tsl,'tsl.mat'));
    tsl=temp.tsl;
end
dx=x(1,2)-x(1,1);   % dx of Tsl

%%% Radargrams
temp=load(fullfile(pfad_rad,'global_coords.mat'));
coords=temp.global_coords;
temp=load(fullfile(pfad_rad,'radargrams.mat'));
radar=temp.radargrams;
temp=load(fullfile(pfad_rad,'t.mat'));
tvec=temp.t;

% List of profile names
profilliste=int2str([1:length(coords)]');


% plot in GUI
plot_TslRadargram(tsl,x,y,t,radar,tvec,coords,profilliste,timedepth);

waitfor(gcf);

% set original path
path(oldpath);