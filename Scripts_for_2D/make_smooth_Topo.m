clear all
close all
clc


% creates a smooth topography from profile data -> is used for
% correcting the height in global_coords for each profile (e.g. for
% topomigration)
% All heights from profiles are interpolated on larger area and smoothed.
% Then the heights for each profile are extracted again and put back into
% global_coords.mat.
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires profile data in radargrams.mat format

num_folders=1; % how many folders with radargrams.mat do you want to combine? 
% (e.g. 4 means that you can choose 4 different folders with data)

plot_flag=1; % if =1: plot some figures for checking

dx=0.05; % grid spacing of area grid

% Interpolation of topography
griding=1;  % 1: Griddata (linear interpolation)
% 2: Inverse distance weighting =IDW (use also the following
% parameters: radius and power)
radius=6;   % Radius for IDW and masking of interpolated timeslice (in bins)
power=10;    % Power of IDW

% smoothing area with 2D median filter
msize=55; % filter size in pixel

% smooth extracted profile heights with 'smooth'
mprof=55; % filter size in pixel
prof_num=[1:5:26]; % show these profile numbers for checking correctness (only if plot_flag==1)

%% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name
for i=1:num_folders
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
                pfad{i}=uigetdir(fn{1}{1},'Choose radargrams.mat folder');
            else
                pfad{i}=uigetdir([],'Choose radargrams.mat folder');
            end
            fid=fopen('temp.temp','wt');
            fprintf(fid,'%s',pfad{i});
            fclose(fid);
        else
            pfad{i}=uigetdir([],'Choose radargrams.mat folder'); % path to radargram-folder

            fid=fopen('temp.temp','wt');
            fprintf(fid,'%s',pfad{i});
            fclose(fid);
        end
    else
        if exist('.temp.temp') % read last opened folder from temp.temp
            fid=fopen('.temp.temp','r');
            fn=textscan(fid,'%s');
            fclose(fid);
            if ~isempty(fn{1})
                pfad{i}=uigetdir(fn{1}{1},'Choose radargrams.mat folder');
            else
                pfad{i}=uigetdir([],'Choose radargrams.mat folder');
            end
        else
            pfad{i}=uigetdir([],'Choose radargrams.mat folder'); % path to radargram-folder
        end

        fid=fopen('.temp.temp','wt');
        fprintf(fid,'%s',pfad{i});
        fclose(fid);
    end
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'),fullfile(curFold,'Processing'));


% read coordinates:
coords=[];
for j=1:length(pfad) % for each folder:
    % check if there is global_coords.mat or several global_coords_*.mat
    if exist(fullfile(pfad{j},'global_coords.mat'),'file')
        load(fullfile(pfad{j},'global_coords.mat'));
        for i=1:length(global_coords)
            coords=[coords; global_coords{i}];
        end
        flag=0;
    else
        % use datastore:
        readFcn = @(f) load(f).global_coords;
        fds = fileDatastore(fullfile(pfad{j},'global_coords_*.mat'), "ReadFcn", readFcn);
        while hasdata(fds)
            global_coords=read(fds);
            for i=1:length(global_coords)
                coords=[coords; global_coords{i}];
            end
        end
        flag=1;
    end
end

% binning:
[topobin,xgrid,ygrid]=bindata2(coords(:,3),coords(:,1),coords(:,2),min(coords(:,1))+dx/2:dx:max(coords(:,1))-dx/2,min(coords(:,2))+dx/2:dx:max(coords(:,2))-dx/2);

if plot_flag==1
    figure
    imagesc(xgrid(1,:),ygrid(:,1),topobin)
    xlabel('x')
    ylabel('y')
    title('Binned raw topography')
end

% Interpolation
disp('Interpolation... Please wait.');
if griding==1
    topo=griddata(xgrid(~isnan(topobin)),ygrid(~isnan(topobin)),topobin(~isnan(topobin)),xgrid,ygrid);
elseif griding==2
    topo=idw2d_tsl(topobin,xgrid,ygrid,eucmap,radius,power);
end

if plot_flag==1
    figure
    imagesc(xgrid(1,:),ygrid(:,1),topo)
    xlabel('x')
    ylabel('y')
    title('Interpolated topography')
end

% smoothing
disp('Smoothing with 2D median filter...')
topo=medianfilt2(topo,[msize msize]);
% mask:
mask=zeros(size(xgrid));
mask(~isnan(topobin))=1;
temp=ones(size(mask));
temp(mask==1)=0;
eucmap=chamfer_DT(temp);
mask_interp=ones(size(eucmap));
mask_interp(eucmap.*dx>dx*3)=0;

topo=topo.*mask_interp;

if plot_flag==1
    figure
    imagesc(xgrid(1,:),ygrid(:,1),topo)
    xlabel('x')
    ylabel('y')
    title('Interpolated & smoothed topography')
end


% creating new z-coordinates from smoothed topography
disp('Creating new height coordinates from smoothed topography')
F=scatteredInterpolant(xgrid(~isnan(topo)),ygrid(~isnan(topo)),topo(~isnan(topo)));
for j=1:num_folders
    disp(['  Folder ',int2str(j),'/',int2str(num_folders)])
    if ~exist(fullfile(pfad{j},'original_coords'),'dir')
        mkdir(fullfile(pfad{j},'original_coords'));
    end
    if flag==0
        copyfile(fullfile(pfad{j},'global_coords.mat'),fullfile(pfad{j},'original_coords','global_coords.mat')); % save original coords file as backup
        load(fullfile(pfad{j},'global_coords.mat'));
        for i=1:length(global_coords)
            if any(i==prof_num) && plot_flag==1
                figure
                plot(global_coords{i}(:,3),'Linewidth',2,'DisplayName','original heights')
                hold on
                temp=F(global_coords{i}(:,1),global_coords{i}(:,2)); % get new heights
                plot(temp,'Linewidth',2,'DisplayName','new heights')
                global_coords{i}(:,3)=smooth(temp,mprof); % smoothing
                plot(global_coords{i}(:,3),'Linewidth',2,'DisplayName','smoothed new heights')
                legend
                title(['Folder ',int2str(j),'/Profile ',int2str(prof_num(i==prof_num))])
            else
                global_coords{i}(:,3)=smooth(F(global_coords{i}(:,1),global_coords{i}(:,2)),mprof);
            end
        end
        save(fullfile(pfad{j},'global_coords.mat'),'global_coords','-v7.3');
    else
        % datastore...
        list=dir(fullfile(pfad{j},'global_coords_*.mat'));
        for k=1:length(list)
            disp(['    - File ',int2str(k),'/',int2str(length(list))])
            copyfile(fullfile(pfad{j},list(k).name),fullfile(pfad{j},'original_coords',list(k).name)); % save original coords file as backup
            global_coords=load(fullfile(pfad{j},list(k).name)).global_coords;
            for i=1:length(global_coords)
                if any(i==prof_num) && plot_flag==1
                    temp=F(global_coords{i}(:,1),global_coords{i}(:,2)); % get new heights
                    global_coords{i}(:,3)=smooth(temp,mprof); % smoothing

                    figure
                    plot(global_coords{i}(:,3),'Linewidth',2,'DisplayName','original heights')
                    hold on
                    plot(temp,'Linewidth',2,'DisplayName','new heights')
                    plot(global_coords{i}(:,3),'Linewidth',2,'DisplayName','smoothed new heights')
                    legend
                    title(['Folder ',int2str(j),'/Profile ',int2str(prof_num(i==prof_num))])
                else
                    global_coords{i}(:,3)=smooth(F(global_coords{i}(:,1),global_coords{i}(:,2)),mprof);
                end
            end
            save(fullfile(pfad{j},list(k).name),'global_coords','-v7.3');
        end
    end
end


% set original path
path(oldpath);
