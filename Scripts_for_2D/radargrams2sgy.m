clear all
close all
clc

% Script for exporting of processed radargrams (from Mala3D) or
% radargrams.mat to sgy files (e.g. for Kingdom Suite)
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% requires folder Export_Import
%
% You can either select the rSlicer-folder (=original processed files will be used, Mala3D had to be run before!) or a folder containing
% radargrams.mat and corresponding files (=any radargrams, could also be from DZT or RDT import).

constoff=0; % if =1: reduce coordinates by constant offset and save in mm accuracy, otherwise full coordinates and cm accuracy

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!

% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder or folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose rSlicer folder or folder with radargrams.mat');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Choose rSlicer folder or folder with radargrams.mat'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder or folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose rSlicer folder or folder with radargrams.mat');
        end
    else
        foldername=uigetdir([],'Choose rSlicer folder or folder with radargrams.mat'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'));

% get all available profile numbers
cd(foldername);
out=dir('*.rad');
out=out(~startsWith({out.name}, '.'));
if ~isempty(out) % -> rSlicer folder -> original processed profiles
    fname=out(1).name; % complete name and number.rad
    name=fname(1:end-8); % only file name without _000.rad
    for i=1:length(out)
        n=out(i).name;  % name of file
        profilelist(i,1)=str2num(n(end-6:end-4)); % profile number
    end
    
    % load data
    temp=load(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'));
    profileinfo=temp.profileinfo;
    
    for i=1:length(profilelist) % for all profiles
        disp(['Exporting profile ',int2str(profilelist(i,1)),'...'])
        m=matfile(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(profilelist(i)),'.mat']));
        temp=load(fullfile(foldername,'profiles2mat',[name,'_',int2str(profilelist(i)),'_info.mat']));
        info=temp.info;
        for j=1:profileinfo(i,4)    % for all channels
            traces=m.traces(:,(j-1)*profileinfo(profileinfo(:,1)==profilelist(i),5)+1:j*profileinfo(profileinfo(:,1)==profilelist(i),5));
            dt=profileinfo(i,2);
            ns=length(traces(:,1));
            x=info(4,info(3,:)==j)';
            y=info(5,info(3,:)==j)';
            z=info(6,info(3,:)==j)';
            t=0:dt:(ns-1)*dt;

            if ~exist(fullfile(foldername,'profiles2mat','proc','sgy'),'dir')
                mkdir(fullfile(foldername,'profiles2mat','proc','sgy'));
            end
            
            disp(['   Channel ',int2str(j)])
            export2sgy2D(traces,dt,x,y,fullfile(foldername,'profiles2mat','proc','sgy',[fname,'_',int2str(profilelist(i)),'_ch',int2str(j),'.sgy']),z,constoff);
        end
    end
    
else  % folder with radargrams.mat, ...
    if exist(fullfile(foldername,'radargrams.txt'),'file')
        fid=fopen(fullfile(foldername,'radargrams.txt'),'r');
        % test for which radargrams.txt-file (DZT_Convert or Mala)
        test=fscanf(fid,'%s',1);
        if strcmp(test,'Radargrams.mat') % Mala
            temp=textscan(fid,'%f','Headerlines',1);
            temp2=textscan(fid,'%f','Headerlines',1);
            chanlist=temp{1}';
            for i=1:length(temp2)
                profilelist=temp2{1}';
            end
            name='dummyname';
        else  % DZT_Convert
            temp=textscan(fid,'%f%s%f','Headerlines',1);
            chanlist=unique(temp{3}');
            plist=temp{2};
            for i=1:length(temp{2})
                c=strsplit(plist{i},'_','CollapseDelimiters',1);
                profilelist(i)=str2num(c{end});
            end
            name=c{1};
        end
        fclose(fid);
    else
        load(fullfile(foldername,'x.mat'));
        profilelist=[1:length(x)];
        chanlist=1;
        name='dummyname';
    end
 
    % load and export data
    radargrams=load(fullfile(foldername,'radargrams.mat'));
    radargrams=radargrams.radargrams;
    t=load(fullfile(foldername,'t.mat'));
    t=t.t;
    x=load(fullfile(foldername,'x.mat'));
    x=x.x;
    global_coords=load(fullfile(foldername,'global_coords.mat'));
    global_coords=global_coords.global_coords;
    for i=1:length(profilelist) % for all profiles
        disp(['Exporting profile ',int2str(profilelist(i)),'...'])
        for j=1:length(chanlist)    % for all channels
            traces=radargrams{(i-1)*length(chanlist)+j};
            dt=t(2)-t(1);
            x=global_coords{(i-1)*length(chanlist)+j}(:,1);
            y=global_coords{(i-1)*length(chanlist)+j}(:,2);
            if length(global_coords{(i-1)*length(chanlist)+j}(1,:))==3
                z=global_coords{(i-1)*length(chanlist)+j}(:,3);
            else
                z=zeros(size(x));
            end

            if ~exist(fullfile(foldername,'sgy'),'dir')
                mkdir(fullfile(foldername,'sgy'));
            end
            
            disp(['   Channel ',int2str(chanlist(j))])
            export2sgy2D(traces,dt,x,y,fullfile(foldername,'sgy',[name,'_',int2str(profilelist(i)),'_ch',int2str(chanlist(j)),'.sgy']),z,constoff);
        end
    end

end

% set original path
path(oldpath);
