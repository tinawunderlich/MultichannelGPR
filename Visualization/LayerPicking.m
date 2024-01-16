clear all
close all
clc


% Picking of layers in radargrams
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires radargrams.mat, x.mat, global_coords.mat and t.mat/z.mat

% time or depth?
tzflag=1;   % time=1, depth=2 (only important for ylabel of radargrams and label for saving of picks)

% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

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
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'));

% Load data
load(fullfile(pfad,'radargrams.mat'));
load(fullfile(pfad,'x.mat'));
load(fullfile(pfad,'global_coords.mat'));
load(fullfile(pfad,'t.mat')); % migth also contain z

% Start GUI
plot_layerpicking(radargrams,global_coords,x,t,tzflag,pfad);

waitfor(gcf);

% set original path
path(oldpath);