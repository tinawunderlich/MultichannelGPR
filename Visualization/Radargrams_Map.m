clear all
close all
clc


% Display Radargrams and corresponding map with connected plots (Useful for
% larger profile spacing when timeslices are not available)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Requires processed radargrams: Radargrams could be from script make_Radargrams,
% Bins2radargrams, Mala3D or own creation (folder should contain radargrams.mat
% (cells with radargrams), global_coords.mat (cells with starting end
% ending coordinates of profile (for straight profiles) or coordinates for 
% each trace(for curved profiles)), t.mat (time vector)).
%
% You have to give the location of the folder containing radargrams.mat and
% corresponding files.
%
% requires content of Plotting and Subfunctions folder

timedepth=1; % 1=time, 2=depth

% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

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

%%% load Radargrams
temp=load(fullfile(pfad_rad,'global_coords.mat'));
coords=temp.global_coords;
temp=load(fullfile(pfad_rad,'radargrams.mat'));
radar=temp.radargrams;
temp=load(fullfile(pfad_rad,'t.mat'));
tvec=temp.t;

% List of profile names
profilliste=int2str([1:length(coords)]');

% plot in GUI
plot_MapRadargram(radar,tvec,coords,profilliste,timedepth);

waitfor(gcf);

% set original path
path(oldpath);