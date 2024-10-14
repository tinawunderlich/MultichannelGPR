% Script for reading Mala-datafiles (prepared for rSlicer = all channels in
% one file) and exporting them to sgy
%
% Dr. Tina Wunderlich, CAU Kiel, 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files:
% *.rd3: Radar data
% *.rad: header file
% *.pos: position data files

clear all
close all
clc

% Select number of profiles:
profile_min=0;  % minimum profile number
profile_max=5;  % maximum profile number
% number of channels for this dataset
channels=16; % number of channels

changeDir=0; % if =1: change the sign of the y-antenna-GPS-offset, if =0: use offsets as written in file

add_Yoffset=0;    % add a constant offset to the y-antenna-GPS-offset (e.g. due to non-vertical GPS-stick), set =0 if not applicable
                    % negative if GPS is tilted towards front, positive if
                    % tilted towards back

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');


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
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            foldername=uigetdir([],'Choose rSlicer folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            foldername=uigetdir([],'Choose rSlicer folder');
        end
    else
        foldername=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% get name
temp=dir(fullfile(foldername,'/*.rad'));
tempname=strsplit(temp(end).name,'_'); % Name of data files without '_???.rd3'
name=[tempname{1}];
for i=2:length(tempname)-1
    name=[name,'_',tempname{i}];
end


% file name parts
[pathstr,fname,ext]=fileparts(name);    % divide name into parts


%% Load all profiles and save as sgy files

% Profile numbers
numbers=profile_min:profile_max;


disp('Read original profiles and save as sgy -> folder raw2syg')
% make folder to save raw data
if ~exist(fullfile(foldername,'raw2sgy'),'dir')
    mkdir(fullfile(foldername,'raw2sgy'));
end
% read profile data and coordinates
not=zeros(1,length(numbers));
lnum=length(numbers);

for i=1:lnum
    disp(['File ',int2str(i),'/',int2str(lnum)])
    % load data and coordinates
    [traces,dt,ns,tempx{i},tempy{i},tempz{i},channels]=readmala4parfor(foldername,name,numbers(i),changeDir,add_Yoffset);
    % delete traces with NaN-coordinates
    del=find(isnan(tempx{i}));
    traces(:,del)=[];
    tempx{i}(del)=[];
    tempy{i}(del)=[];
    tempz{i}(del)=[];
    if ~isempty(dt)
        % resort into matrices (without cells)
        numtr=length(tempx{i}); % number of all traces per profile
        numtrch=numtr/channels; % number of traces per channel
        info=zeros(9,numtr);   % profilenum, tracenum per channel, channelnum, x, y, z, dt, ns (per trace), tracenum per profile
        info(1,:)=[zeros(1,length(tempx{i}))+numbers(i)]; % profilenumber
        info(4:6,:)=[tempx{i}; tempy{i}; tempz{i}]; % x,y,z
        info(9,:)=[1:numtr]; % tracenumber per profile
        for ii=1:channels
            info(2:3,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[1:numtrch; zeros(1,numtrch)+ii]; % trace number per channel, channelnumber
            info(7:8,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[zeros(1,numtrch)+dt; zeros(1,numtrch)+ns]; % dt, ns
        end
        % save profile (raw data)
        for j=1:channels % one file for each channel and each profile
            if numbers(i)<10
                sgyname=fullfile(foldername,'raw2sgy',[name,'_00',int2str(numbers(i)),'_C',int2str(j),'.sgy']);
            elseif numbers(i)<100
                sgyname=fullfile(foldername,'raw2sgy',[name,'_0',int2str(numbers(i)),'_C',int2str(j),'.sgy']);
            else
                sgyname=fullfile(foldername,'raw2sgy',[name,'_',int2str(numbers(i)),'_C',int2str(j),'.sgy']);
            end
            x=info(4,info(3,:)==j);
            y=info(5,info(3,:)==j);
            z=info(6,info(3,:)==j);
            export2sgy2D(traces(:,info(3,:)==j),info(7,1),x,y,sgyname,z,0);
        end
    end
end

disp('Finished!')