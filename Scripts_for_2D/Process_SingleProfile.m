clear all
close all
clc


% Read single profile of Mala Mira data (prepared for rSlicer), Spidar raw data or
% radargrams.mat (=import into MultichannelGPR from various systems) and do processing individually
% Only for testing of processing options, no saving of data, just saving of settings-file!
%
% Dr. Tina Wunderlich, CAU Kiel 2020-2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires folders Export_Import, Processing, Subfunctions, Migration,
% Plotting


% -------------------------------------------------------------------------
% Do not change the following part! 

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
            folder=uigetdir(fn{1}{1},'Select data-folder');
        else
            folder=uigetdir([],'Select data-folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',folder);
        fclose(fid);
    else
        folder=uigetdir([],'Select data-folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',folder);
        fclose(fid);
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            folder=uigetdir(fn{1}{1},'Select data-folder');
        else
            folder=uigetdir([],'Select data-folder');
        end
    else
        folder=uigetdir([],'Select data-folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',folder);
    fclose(fid);
end



% temporarily set path to required scripts
oldpath=path;
addpath('../Processing/','../Export_Import/','../Subfunctions/','../Migration/','../GUIs/');

%% MALA
temp=dir(fullfile(folder,'/*.rad')); % get list of all Mala Mira files
if (~isempty(temp)) % if Mala mira data available
    disp('Mala MIRA data found. Please wait!')
    anz=1;
    for i=1:length(temp)
        if ~startsWith(temp(i).name,'.')
            profilelist(anz)=str2num(temp(i).name(end-6:end-4)); % profile number
            anz=anz+1;
        end
    end

    tempname=strsplit(temp(end).name,'_'); % Name of data files without '_???.rd3'
    name=[tempname{1}];
    for i=2:length(tempname)-1
        name=[name,'_',tempname{i}];
    end

    processingTestGUI(folder,name,profilelist,[]);
end

%% SPIDAR
temp=dir(fullfile(folder,'/*.DT1')); % get list of all spidar files
if (~isempty(temp)) % if spidar data available
    disp('Spidar data found. Please wait!')
    chnum=[];
    pnum=[];
    for i=1:length(temp)
        if ~startsWith(temp(i).name,'.')
            temp1=strsplit(temp(i).name,'_'); % split in two parts (e.g. NIC01 and Line001.DT1)
            chnum=[chnum; str2double(temp1{1}(end-1:end))]; % channelNumber
            temp2=strsplit(temp1{2},'.'); % split in tow parts: e.g. Line001 and DT1
            pnum=[pnum; str2double(temp2{1}(end-2:end))]; % profile number
            if ~exist('name','var')
                name{1}=temp1{1}(1:end-2); % project name (e.g. NIC)
                name{2}=temp2{1}(1:end-3); % e.g. Line
            end
        end
    end
    profilelist=unique(pnum);
    channellist=unique(chnum);

    processingTestGUI(folder,name,profilelist,channellist);
end


%% MultichannelGPR format
temp=dir(fullfile(folder,'/*.mat')); % get list of all mat-files
if (~isempty(temp)) % if mat data available
    % check if all necessary files are available
    if any(cellfun(@(x) strcmp(x,'radargrams.mat'),{temp.name})) && ...
            any(cellfun(@(x) strcmp(x,'x.mat'),{temp.name})) && ...
            any(cellfun(@(x) strcmp(x,'t.mat'),{temp.name})) && ...
            any(cellfun(@(x) strcmp(x,'global_coords.mat'),{temp.name}))
        disp('*.mat-files found. Please wait!')
        load(fullfile(folder,'radargrams.mat'));
        load(fullfile(folder,'x.mat'));
        load(fullfile(folder,'t.mat'));
        load(fullfile(folder,'global_coords.mat'));

        processingTestGUI_mat(radargrams,t,x,global_coords,folder);
    else
        disp('mat-files found, but data is not complete. Please check again your folder.')
    end
end

waitfor(gcf);

% restore original path
path(oldpath);