clear all
close all
clc

%------------------------ Proceq_Convert -------------------------------------
%
% Converts sgy-files of Proceq equipment to mat-files (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, April 2025, tina.wunderlich@ifg.uni-kiel.de, CAU
% Kiel
%
%--------------------------------------------------------------------------


dataplot=1; % plot radargram for controlling? 1=yes, 0=no

% Coordinate settings:
filter_coords=1; % if =1: coordinates are filtered and traces at
% approx. the same position will be removed (e.g. at the beginning or
% end of the profile): a local difference
% over m coordinates for x and y is calculated separately. If this
% absolute local difference is smaller than nstd*std(diff) it will be removed.
m=15; % number of samples for local difference calculation
nstd=1;  % nstd*standard deviation as valid range
show_filter_plot=1; % if =1: show plot for quality control of filtering
coords_smooth=0; % if=1: smooth coordinates along profile
nsmooth=15;  % number of traces for smoothing of coords
utmzone=32; % give UTM-zone for conversion from WGS84 to UTM

offsetGPS_X=0; % Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=0; % Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)
offsetGPS_Z=0.4; % Offset between GPS and ground surface (always positive) (use only if not corrected inside the GPS acquisition)

% Options for calculating inline coordinates for each trace:
coords_opt=2;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
% =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)



%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
%
% LOAD DATA

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));


% get file names
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
            pfad=uigetdir(fn{1}{1},'Choose folder with sgy-folder(s)');
        else
            pfad=uigetdir([],'Choose folder with sgy-folder(s)');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with sgy-folder(s)'); % path to radargram-folder

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
            pfad=uigetdir(fn{1}{1},'Choose folder with sgy-folder(s)');
        else
            pfad=uigetdir([],'Choose folder with sgy-folder(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with sgy-folder(s)'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*/*.sgy'));
% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end

% go through list and extract infos
starttime=zeros(length(list),1);
for i=1:length(list)
    parts=strsplit(list(i).name,'_');
    list(i).givenname=parts{1};
    starttime(i)=str2num(parts{3});
    list(i).LF_HF=extractBefore(parts{6},'.sgy');
end

% sort list
[a,b]=sort(starttime);

LF_List=[];
HF_List=[];
for i=1:length(b)
    if strcmp(list(b(i)).LF_HF,'LF')
        LF_List=[LF_List; list(b(i))];
    else
        HF_List=[HF_List; list(b(i))];
    end
end

%% LF:
% go through list and read data
for i=1:length(LF_List)
    disp(['File ',int2str(i),' of ',int2str(length(LF_List))])

    filename=LF_List(i).name;
    [~,name{i},ext]=fileparts(filename);

    [traces,headers{i},coords]=readProceq(fullfile(LF_List(i).folder,filename));

    if ~isempty(traces)
        flag_ok(i)=1;
        % set in variables:
        radargrams{i}=traces;
        global_coords{i}=coords;
        xtemp=coords(:,1);
        ytemp=coords(:,2);
        % calculate profile coordinates
        if coords_opt==1 % difference to beginning
            x{i}=[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)'];
        elseif coords_opt==2 % cumulative sum
            x{i}=[0; reshape(cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)'),length(xtemp)-1,1)];
        end
        t=double(0:headers{i}.dt:headers{i}.dt*(headers{i}.ns-1)).*1e-3; % in ns

        if dataplot==1
            figure
            imagesc(x{i},t,radargrams{i})
            title(name{i})
            xlabel('x [m]')
            ylabel('t [ns]')
            colormap(flipud(gray))
        end
    else
        flag_ok(i)=0; % empty traces because of bad coordinates
    end
end

% write data
% delete empty cells
ok=~cellfun('isempty',radargrams);
radargrams=radargrams(ok);
x=x(ok);
global_coords=global_coords(ok);
headers=headers(ok);

disp('Start saving of LF mat-files...')


% make new folder
if ~exist(fullfile(pfad,'mat_LF'),'dir')
    mkdir(fullfile(pfad,'mat_LF'));
end
save(fullfile(pfad,'mat_LF','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(pfad,'mat_LF','t.mat'),'t','-v7.3');
save(fullfile(pfad,'mat_LF','x.mat'),'x','-v7.3');
save(fullfile(pfad,'mat_LF','global_coords.mat'),'global_coords','-v7.3');
save(fullfile(pfad,'mat_LF','h.mat'),'headers','-v7.3');

fid=fopen(fullfile(pfad,'mat_LF','radargrams.txt'),'wt');
fprintf(fid,'Nr.\tName\n');
anz=1;
for i=1:length(name)
    if flag_ok(i)==1
        fprintf(fid,'%d\t',anz);
        fprintf(fid,'%s\n',[name{i},'.sgy']);
        anz=anz+1;
    end
end
fclose(fid);

%% HF:
% go through list and read data
for i=1:length(HF_List)
    disp(['File ',int2str(i),' of ',int2str(length(HF_List))])

    filename=HF_List(i).name;
    [~,name{i},ext]=fileparts(filename);

    [traces,headers{i},coords]=readProceq(fullfile(HF_List(i).folder,filename));

    if ~isempty(traces)
        flag_ok(i)=1;
        % set in variables:
        radargrams{i}=traces;
        global_coords{i}=coords;
        xtemp=coords(:,1);
        ytemp=coords(:,2);
        % calculate profile coordinates
        if coords_opt==1 % difference to beginning
            x{i}=[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)'];
        elseif coords_opt==2 % cumulative sum
            x{i}=[0; reshape(cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)'),length(xtemp)-1,1)];
        end
        t=double(0:headers{i}.dt:headers{i}.dt*(headers{i}.ns-1)).*1e-3; % in ns

        if dataplot==1
            figure
            imagesc(x{i},t,radargrams{i})
            title(name{i})
            xlabel('x [m]')
            ylabel('t [ns]')
            colormap(flipud(gray))
        end
    else
        flag_ok(i)=0; % empty traces because of bad coordinates
    end
end

% write data
% delete empty cells
ok=~cellfun('isempty',radargrams);
radargrams=radargrams(ok);
x=x(ok);
global_coords=global_coords(ok);
headers=headers(ok);

disp('Start saving of HF mat-files...')


% make new folder
if ~exist(fullfile(pfad,'mat_HF'),'dir')
    mkdir(fullfile(pfad,'mat_HF'));
end
save(fullfile(pfad,'mat_HF','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(pfad,'mat_HF','t.mat'),'t','-v7.3');
save(fullfile(pfad,'mat_HF','x.mat'),'x','-v7.3');
save(fullfile(pfad,'mat_HF','global_coords.mat'),'global_coords','-v7.3');
save(fullfile(pfad,'mat_HF','h.mat'),'headers','-v7.3');

fid=fopen(fullfile(pfad,'mat_HF','radargrams.txt'),'wt');
fprintf(fid,'Nr.\tName\n');
anz=1;
for i=1:length(name)
    if flag_ok(i)==1
        fprintf(fid,'%d\t',anz);
        fprintf(fid,'%s\n',[name{i},'.sgy']);
        anz=anz+1;
    end
end
fclose(fid);

disp('   Finished!')

% set original path
path(oldpath);
