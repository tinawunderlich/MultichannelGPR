clear all
close all
clc

% Script for interactive testing of different options for t0-corretion,
% using Mala Mira raw data in rSlicer folder or other data in the correct
% mat-format
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires folders Plotting, Export_Import, Processing
%
% in case of DZT-Data the utm zone needs to be specified & Kann das weg????
UTMzone = 32; 

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
% get folder name

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
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Select rSlicer-folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder');
        end
    else
        foldername=uigetdir([],'Select rSlicer-folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Processing'),fullfile(curFold,'Subfunctions'),fullfile(curFold,'Export_Import'),fullfile(curFold,'GUIs'));


%% in case of Mira Mala Data
% get all available profile numbers
out=dir(fullfile(foldername,'/*.rad'));
if(~isempty(out))
    anz=1;
    for i=1:length(out)
        if ~startsWith(out(i).name,'.')
            n=out(i).name;  % name of file
            profilelist(anz,1)=str2num(n(end-6:end-4)); % profile number
            anz=anz+1;
        end
    end
    
    tempname=strsplit(n,'_'); % Name of data files without '_???.rd3'
    name=[tempname{1}];
    for i=2:length(tempname)-1
        name=[name,'_',tempname{i}];
    end
    
    % plot with gui options
    Check_t0_slider_plot(foldername,name,profilelist);
end

%% in case of GSSI-Data
out=dir(fullfile(foldername,'/*.DZT'));
if(~isempty(out))
    anz=1;
    for i=1:length(out)
        if ~startsWith(out(i).name,'.')
            n=out(i).name;  % name of file
            profilelist(anz,1)=str2num(n(end-6:end-4)); % profile number
            anz=anz+1;
        end
    end
    tempname=strsplit(n,'_'); % Name of data files without '_???.rd3'
    name=[tempname{1}];
    for i=2:length(tempname)-1
        name=[name,'_',tempname{i}];
    end
    
    % plot with gui options
    Check_t0_slider_plot_gssi(foldername,name,profilelist,UTMzone);
end

% check if figure is still open
fig=findobj('type','figure','name','Traces');
uiwait(fig); % wait until figure is closed

% set original path
path(oldpath);
