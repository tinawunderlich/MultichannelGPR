% Script for reading radargrams.mat and associated files and
% binning them onto a rectangular grid (for each profile individually to
% get balanced channel energies and less stripes in timeslices)
% Select folder with radargrams.mat!
%
% Dr. Tina Wunderlich, CAU Kiel 2025-2026, tina.wunderlich@ifg.uni-kiel.de
% OPTIMIZED VERSION (with help of Claude.ai free version, manually checked!)
%
% requires MATLAB-files in following folders (path will be temporarily
% set):  Subfunctions

clear all
close all
clc

% Bin size of grid
dx=0.05; % [m]

radius=0.3; % radius in m for valid interpolation (-> mask)

% Automatic rotation of measurement area for minimum memory size
rotate_area=1;  % 1=yes (recommended), 0=no

% time or depth?
tz_flag=2;  % 1: time -> [ns], 2: depth -> [m]

% if depth: follow Topography or horizontal slices?
followTopo=0; % =1: yes; =0: horizontal slices
constV=0.1; % constant velocity [m/ns], only used if tz_flag=1 and followTopo=1
%%% Explanation:
% If tz_flag=1 and followTopo=0: ignore the topography and make horizontal
% slices parallel to t=0ns.
% If tz_flag=1 and followTopo=1: normalize data for each time step and profile, 
% use the constant velocity to convert time to depth,  and then use the
% topography to bend the slices parallel to topography. Then cut again
% horizontal slice out of the 3D cube (cutting the topography).
% If tz_flag=2 and followTopo=0: Input is topography-migrated data.
% Normalization of depth samples is done parallel to topography, but slices
% are cut horizontally (cutting the topography).
% If tz_flag=2 and followTopo=1: Input is topography-migrated data.
% Normalization of depth samples is done parallel to topography, and slices
% are also calculated parallel to topography.

% smoothing of topography along profiles (only if tz_flag=1. For tz_flag=2
% do smoothing before topo migration!)
n_1dTopo=155; % number of traces along profile for median filter

% smoothing of topography on map (2D) -> only for topo display
n_topo=15; % number of grid cells for 2d median filter

% Downsampling of data
downsampling=1; % if =1: yes (and use following settings)
downsampling_factor=20; % only take each downsampling-factor sample (e.g. only take every 2nd sample)
% cutting of range
cut_range=0; % if =1:yes
cut_time_depth=30; % choose time/depth for cutting (tz_flag=1&followTopo=0: [ns] or all other combinations: [m])

% save Sampleslices as geopng or geotif?
save_geopng=0; % 1=yes
save_geotif=1; % 1=yes
colperc=3; % Colorscale clipping in percent (if =0: autoscale min-max)
removeBorder=1; % =1: remove border artifacts from interpolation, =0: leave as it is
pix=6; % if removeBorder==1: how many pixels are removed from border around area
medianFilter=1; % do you want to apply a 2D-median filter (1=yes, 0=no)
msize=3; % filter size in pixel
% use squareroot of amplitudes for visualization?
sq=0; % 1=yes, 0=no
% for geotif: epsg code of CRS
epsg=25832;


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name
if ispc
    if exist('temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose folder with radargrams.mat');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    end
else
    if exist('.temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose folder with radargrams.mat');
        end
    else
        foldername=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end



% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'),fullfile(curFold,'Export_Import'));



%% check if correctly processed data available:
if exist(fullfile(foldername,'radargrams.mat'),'file') % all in one file
    % load profileinfo
    m=matfile(fullfile(foldername,'radargrams.mat'));
    load(fullfile(foldername,'x.mat'));
    load(fullfile(foldername,'global_coords.mat'));

    numtraces=cellfun(@(xx) numel(xx),x);
    numbers=1:numel(x); % profile numbers
    disp(['Found processed data of ',int2str(length(numbers)),' profiles.'])

    datastore_flag=0;
elseif exist(fullfile(foldername,'radargrams_1.mat'),'file') % several files -> datastore
    datastore_flag=1;

    % get all files (Radargrams)
    files = dir(fullfile(foldername, 'radargrams_*.mat'));
    % sort numerically:
    nums = arrayfun(@(f) sscanf(f.name, 'radargrams_%d.mat'), files);
    [~, sortIdx] = sort(nums);
    files = files(sortIdx);
    % with sorted list: 1, 2, 3, .... (instead of 1, 10, 11, ..., 19, 2, 20, ...)
    filePaths = fullfile(foldername, {files.name});

    % get all files (global_coords)
    files = dir(fullfile(foldername, 'global_coords_*.mat'));
    % sort numerically:
    nums = arrayfun(@(f) sscanf(f.name, 'global_coords_%d.mat'), files);
    [~, sortIdx] = sort(nums);
    files = files(sortIdx);
    % with sorted list: 1, 2, 3, .... (instead of 1, 10, 11, ..., 19, 2, 20, ...)
    filePaths_gc = fullfile(foldername, {files.name});

    radar_ds = fileDatastore(filePaths, "ReadFcn", @(f) readRadargrams(f));
    gc_ds = fileDatastore(filePaths_gc, "ReadFcn", @load);
    anz=numel(radar_ds.Files);
else
    disp('No data found.')
    return;
end
load(fullfile(foldername,'t.mat'));
dt=abs(t(2)-t(1));
ns=numel(t);

if tz_flag==2 % if depth input
    z_abs=t; % t contains absolute depth (=height) in m
    t=0:dt:(ns-1)*dt; % z-vector from 0 (top) to maxz (bottom)
end

% smoothing of topo along profiles:
if tz_flag==1 % only if time domain
    if datastore_flag==0
        global_coords = cellfun(@(c) [c(:,1:2), medfilt1_own(c(:,3), n_1dTopo)], global_coords, 'UniformOutput', false);
    end
end

%% Size of area
% read info files for coordinates
if datastore_flag==0
    nRows = cellfun(@(c) size(c, 1), global_coords); % number of traces per radragram
    idx = repelem(1:numel(global_coords), nRows)'; % profile number
    coords = vertcat(global_coords{:});
    temp=arrayfun(@(n) (1:n)', nRows, 'UniformOutput', false);
    rowIdx = vertcat(temp{:}); % tracenumber in profile
    % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
    xylist=zeros(numel(idx),7);
    xylist(:,1)=idx;
    xylist(:,2:4)=coords;
    xylist(:,5)=ones(size(xylist(:,5)));
    xylist(:,6)=rowIdx;
    xylist(:,7)=xylist(:,6);
else
    reset(gc_ds);
    i=1;
    offsetnum=0;
    xylist_temp=cell(anz,1);
    while hasdata(gc_ds)
        temp = read(gc_ds);
        global_coords=temp.global_coords;

        nRows = cellfun(@(c) size(c, 1), global_coords); % number of traces per radragram
        idx = repelem(1:numel(global_coords), nRows)'; % profile number in this file
        coords = vertcat(global_coords{:});
        temp=arrayfun(@(n) (1:n)', nRows, 'UniformOutput', false);
        rowIdx = vertcat(temp{:}); % tracenumber in profile
        % Number, x, y, z, channel of profile, tracenumber in profile,
        % tracenumber in channel, idx of file, number in file
        xylist_temp{i}=zeros(numel(idx),9);
        xylist_temp{i}(:,1)=idx+offsetnum;
        xylist_temp{i}(:,2:4)=coords;
        xylist_temp{i}(:,5)=ones(size(xylist_temp{i}(:,5)));
        xylist_temp{i}(:,6)=rowIdx;
        xylist_temp{i}(:,7)=xylist_temp{i}(:,6);
        xylist_temp{i}(:,8)=i; % file number
        xylist_temp{i}(:,9)=idx; % profile number in this file
        i=i+1; % increase file number
        offsetnum=offsetnum+max(idx); % increase offset number
    end
    xylist=vertcat(xylist_temp{:});
    fprintf('Found %d data files.\n',i-1);
end

profnum=unique(xylist(:,1)); % profile numbers

% optional: rotate area
disp('Reading coordinates for determining area size...')
if rotate_area==1
    [xylist(:,2:3),rotbest,shiftx,shifty,coordtrans]=rotatearea(xylist(:,2:3));

    fig1=figure('Visible','off');
    plot(xylist(:,2),xylist(:,3),'k.')
    hold on
    set(gca,'Dataaspectratio',[1 1 1])
    axis xy
    xlabel('x [m]')
    ylabel('y [m]')
end


%% create sampleslices and balance amplitudes for each channel
% make grids
[xgrid,ygrid]=meshgrid(min(xylist(:,2)):dx:max(xylist(:,2)),min(xylist(:,3)):dx:max(xylist(:,3)));
linearindex=reshape(1:numel(xgrid),size(xgrid)); % linear indices of grid

if tz_flag==1 && followTopo==0
    newfolder='SampleSlices_time_parallel2surface';
elseif tz_flag==1 && followTopo==1
    newfolder='SampleSlices_time_followTopo_cutHorizontally';
elseif tz_flag==2 && followTopo==0
    newfolder='SampleSlices_depth_horizontal';
elseif tz_flag==2 && followTopo==1
    newfolder='SampleSlices_depth_followTopo';
end


if ~exist(fullfile(foldername,newfolder),'dir')
    mkdir(fullfile(foldername,newfolder))
end
% save x/y-grids:
save(fullfile(foldername,newfolder,'xgrid.mat'),'xgrid');
save(fullfile(foldername,newfolder,'ygrid.mat'),'ygrid');

% bin edges:
xrg=min(xylist(:,2))-dx/2:dx:max(xylist(:,2))+dx/2;
yrg=min(xylist(:,3))-dx/2:dx:max(xylist(:,3))+dx/2;

% make topography bins:
topo=bindata2(xylist(:,4),xylist(:,2),xylist(:,3),xrg,yrg);

disp('Smooth and interpolate topography...')
% make mask for topo:
mask_topo=zeros(size(topo));
mask_topo(~isnan(topo))=1;
temp=ones(size(mask_topo));
temp(mask_topo==1)=0;
eucmap=chamfer_DT(temp);
mask_interp_topo=ones(size(eucmap));
mask_interp_topo(eucmap.*dx>radius)=0;
% smoothing with median filter:
topo=medianfilt2(topo,n_topo);
topo(mask_topo==0)=NaN;
F=scatteredInterpolant(xgrid(mask_topo>0),ygrid(mask_topo>0),topo(mask_topo>0)); % re-use this later, because same grid
topo_interp=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
topo_interp(mask_interp_topo==0)=NaN;

if tz_flag==2 
    disp('Get topography indices for each profile...')

    % get info from profiles
    row_ind_start=cell(numel(profnum),1);
    msg = '';
    if datastore_flag==0
        for i=1:numel(profnum) % for each profile:
            % alte Zeile löschen
            fprintf(repmat('\b', 1, length(msg)));
            % neue Zeile schreiben
            msg = sprintf('[%-*s]', length(profnum), repmat('.', 1, i));
            fprintf(msg);

            % load data of this profile
            if size(m.radargrams,1)>1
                temp=m.radargrams(i,1);
            else
                temp=m.radargrams(1,i); % -> traces (all channels)
            end
            traces=temp{1};

            % find first/last datapoint for each column:
            nan_mask = ~isnan(traces);

            [~, row_ind_start{i}] = max(nan_mask, [], 1);    % erste non-NaN
            [~, last_idx]         = max(flipud(nan_mask), [], 1);
            row_ind_end           = size(traces, 1) - last_idx + 1;  % letzte non-NaN

            thickness = abs(t(row_ind_end) - t(row_ind_start{i}));   % [1 x nCols]
            minz(i) = min(thickness);
        end
    else % datastore:
        i=1; % profile number overall
        j=1; % file number
        while hasdata(radar_ds)
            fprintf(' -- File %d (of %d):\n',j,anz);
            data=read(radar_ds);

            msg = '';
            for ii=1:numel(data.radargrams)
                fprintf(repmat('\b', 1, length(msg)));
                msg=sprintf('      Radargram %d (of %d)\n',ii,numel(data.radargrams));
                fprintf('%s', msg);

                % load data of this profile
                traces=data.radargrams{ii};

                % find first/last datapoint for each column:
                nan_mask = ~isnan(traces);

                [~, row_ind_start{i}] = max(nan_mask, [], 1);    % erste non-NaN
                [~, last_idx]         = max(flipud(nan_mask), [], 1);
                row_ind_end           = size(traces, 1) - last_idx + 1;  % letzte non-NaN

                thickness = abs(t(row_ind_end) - t(row_ind_start{i}));   % [1 x nCols]
                minz(i) = min(thickness);

                i=i+1;
            end
            fprintf('\n');
            j=j+1;
        end
        anz_all=j-1; % overall number of profiles over all groups
    end


    % one z-vector for all
    dz=abs(t(2)-t(1)); % the same for all profiles
    z_all=0:dz:min(minz);

    if followTopo==1
        % reduce samples at this point:
        % downsampling/cutting:
        if downsampling==1
            if cut_range==1
                timesamplenum=1:downsampling_factor:length(z_all(z_all<=cut_time_depth));
            else
                timesamplenum=1:downsampling_factor:length(z_all);
            end
        else
            if cut_range==1
                timesamplenum=1:length(z_all(z_all<=cut_time_depth));
            else
                timesamplenum=1:length(z_all);
            end
        end

        t=z_all(timesamplenum); % new depth vector (only named t for consistency!)
    elseif followTopo==0
        % no reduction of samples at this point -> done later!
        t=z_all;
        timesamplenum=1:length(z_all);
    end

    % save time vector and coordtrans:
    save(fullfile(foldername,newfolder,'t.mat'),'t');
    save(fullfile(foldername,newfolder,'coordtrans.mat'),'coordtrans');

    % create inital profnum & channum slices:
    slice_prof=NaN(size(xgrid));

    if followTopo==1 % parallel to topography:
        % read profile data:
        disp('-----------')
        disp(['Get data of profiles ',int2str(profnum(1)),'-',int2str(profnum(end)),'...'])
        fprintf('Profile\tData\tBinning\t\tSaved\tTime elapsed [s]\n')

        % extract and normalize slices:
        if datastore_flag==0
            for n=1:numel(x) % for each profile:
                tstart=tic;
                fprintf('%d\t',profnum(n));

                % load data of this profile
                if size(m.radargrams,1)>1
                    temp=m.radargrams(n,1);
                else
                    temp=m.radargrams(1,n); % -> traces (all channels)
                end
                traces=temp{1};

                % coords of this profile
                ctemp=xylist(xylist(:,1)==profnum(n),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

                fprintf('x\t');

                % get profile data for relevant time samples only, parallel to topography:
                row_indices = row_ind_start{n}(:)' + timesamplenum(:) - 1; % row indices for all columns
                col_indices = repmat(1:size(traces,2), numel(timesamplenum), 1);  % [nRows x nCols]
                linear_idx = sub2ind(size(traces), row_indices, col_indices);
                traces = traces(linear_idx);

                % --- normalization ---
                mu  = mean(traces(:,all(~isnan(traces),1)), 2);
                sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
                sg(sg==0) = 1; % avoid divide-by-zero
                traces = 100 .* (traces-mu) ./ sg;

                % bin data:
                dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
                % valid data points:
                validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
                % initialize variable for all data:
                nValid = numel(validdata);
                nSl    = size(dtemp, 3);
                profiledata  = zeros(nValid, nSl+1);
                profiledata(:,1) = validdata;

                % --- OPTIMIZATION: vectorized slice normalization ---
                % Extract valid pixels for all slices at once: nValid x nSl
                % Reshape dtemp to (ny*nx) x nSl, pick valid rows
                nY = size(dtemp,1); nX = size(dtemp,2);
                flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
                valid_flat = flat(validdata, :);      % nValid x nSl

                mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
                sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
                sg_sl(sg_sl==0) = 1;
                profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

                fprintf('x\t\t');

                % Save profiledata:
                save(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(n)),'.mat']),'profiledata');
                fprintf('x\t');

                fprintf('\t%.1f\n',toc(tstart));
            end
        else % datastore
            i=1; % profile number overall
            j=1; % file number
            reset(radar_ds);
            while hasdata(radar_ds)
                fprintf('  -- File %d (of %d) --\n',j,anz_all);
                data=read(radar_ds);
                j=j+1;

                for ii=1:numel(data.radargrams)
                    tstart=tic;
                    fprintf('%d\t',profnum(i));

                    % load data of this profile
                    traces=data.radargrams{ii};

                    % coords of this profile
                    ctemp=xylist(xylist(:,1)==profnum(i),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

                    fprintf('x\t');

                    % get profile data for relevant time samples only, parallel to topography:
                    row_indices = row_ind_start{i}(:)' + timesamplenum(:) - 1; % row indices for all columns
                    col_indices = repmat(1:size(traces,2), numel(timesamplenum), 1);  % [nRows x nCols]
                    linear_idx = sub2ind(size(traces), row_indices, col_indices);
                    traces = traces(linear_idx);

                    % --- normalization ---
                    mu  = mean(traces(:,all(~isnan(traces),1)), 2);
                    sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
                    sg(sg==0) = 1; % avoid divide-by-zero
                    traces = 100 .* (traces-mu) ./ sg;

                    % bin data:
                    dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
                    % valid data points:
                    validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
                    % initialize variable for all data:
                    nValid = numel(validdata);
                    nSl    = size(dtemp, 3);
                    profiledata  = zeros(nValid, nSl+1);
                    profiledata(:,1) = validdata;

                    % --- OPTIMIZATION: vectorized slice normalization ---
                    % Extract valid pixels for all slices at once: nValid x nSl
                    % Reshape dtemp to (ny*nx) x nSl, pick valid rows
                    nY = size(dtemp,1); nX = size(dtemp,2);
                    flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
                    valid_flat = flat(validdata, :);      % nValid x nSl

                    mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
                    sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
                    sg_sl(sg_sl==0) = 1;
                    profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

                    fprintf('x\t\t');

                    % Save profiledata:
                    save(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(i)),'.mat']),'profiledata');
                   
                    fprintf('x\t');

                    fprintf('\t%.1f\n',toc(tstart));

                    i=i+1;
                end
            end
        end

    elseif followTopo==0 % horizontal slices cutting the topography:
        % read profile data:
        disp('-----------')
        disp(['Get data of profiles ',int2str(profnum(1)),'-',int2str(profnum(end)),' and normalize parallel to topography...'])
        fprintf('Profile\tData\tNormalize\tBinning\tTime elapsed [s]\n')

        % extract and normalize slices:
        if datastore_flag==0
            for n=1:numel(x) % for each profile:
                tstart=tic;
                fprintf('%d\t',profnum(n));

                % load data of this profile
                if size(m.radargrams,1)>1
                    temp=m.radargrams(n,1);
                else
                    temp=m.radargrams(1,n); % -> traces (all channels)
                end
                traces=temp{1};
                [M,N]=size(traces); % original size

                % coords of this profile
                ctemp=xylist(xylist(:,1)==profnum(n),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

                fprintf('x\t');

                % get profile data for ALL time samples, parallel to topography:
                row_indices = row_ind_start{n}(:)' + timesamplenum(:) - 1; % row indices for all columns
                col_indices = repmat(1:size(traces,2), numel(timesamplenum), 1);  % [nRows x nCols]
                linear_idx = sub2ind(size(traces), row_indices, col_indices);
                traces = traces(linear_idx);

                % --- normalization ---
                mu  = mean(traces(:,all(~isnan(traces),1)), 2);
                sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
                sg(sg==0) = 1; % avoid divide-by-zero
                traces = 100 .* (traces-mu) ./ sg;

                fprintf('x\t');

                % put normalized data back to topography bended radargram:
                radargram=NaN(M,N);
                radargram(linear_idx)=traces;
                traces=radargram;
                clear radargram;

                % extract only relevant samples:
                % downsampling/cutting:
                if downsampling==1
                    if cut_range==1
                        tsnum=1:downsampling_factor:length(t(t<=cut_time_depth));
                    else
                        tsnum=1:downsampling_factor:length(t);
                    end
                else
                    if cut_range==1
                        tsnum=1:length(t(t<=cut_time_depth));
                    else
                        tsnum=1:length(t);
                    end
                end
                tnew=t(tsnum); % new depth vector (only named t for consistency!)
                traces=traces(tsnum,:);

                % bin data along horizontal slices:
                dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
                % each depth sluce for this profile has difefrent valid
                % values (due to topo):
                validdata = linearindex(any(~isnan(dtemp), 3)); % indices with at least one non-nan value at any depth
                % initialize variable for all data:
                nValid = numel(validdata);
                nSl    = size(dtemp, 3);
                profiledata  = zeros(nValid, nSl+1);
                profiledata(:,1) = validdata;

                % --- OPTIMIZATION: vectorized slice normalization ---
                % Extract valid pixels for all slices at once: nValid x nSl
                % Reshape dtemp to (ny*nx) x nSl, pick valid rows
                nY = size(dtemp,1); nX = size(dtemp,2);
                flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
                valid_flat = flat(validdata, :);      % nValid x nSl

                mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
                sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
                sg_sl(sg_sl==0) = 1;
                profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

                % Save profiledata:
                save(fullfile(foldername,newfolder,['profiledata_',int2str(numbers(n)),'.mat']),'profiledata');
                fprintf('x\t');

                fprintf('\t%.1f\n',toc(tstart));
            end
        else % datastore:
            i=1; % profile number overall
            j=1; % file number
            reset(radar_ds);
            while hasdata(radar_ds)
                fprintf('   -- File %d (of %d) --\n',j,anz_all);
                data=read(radar_ds);
                j=j+1;

                for ii=1:numel(data.radargrams)
                    % load data of this profile
                    traces=data.radargrams{ii};
                    [M,N]=size(traces); % original size

                    tstart=tic;
                    fprintf('%d\t',profnum(i));

                    % coords of this profile
                    ctemp=xylist(xylist(:,1)==profnum(i),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

                    fprintf('x\t');

                    % get profile data for ALL time samples, parallel to topography:
                    row_indices = row_ind_start{i}(:)' + timesamplenum(:) - 1; % row indices for all columns
                    col_indices = repmat(1:size(traces,2), numel(timesamplenum), 1);  % [nRows x nCols]
                    linear_idx = sub2ind(size(traces), row_indices, col_indices);
                    traces = traces(linear_idx);

                    % --- normalization ---
                    mu  = mean(traces(:,all(~isnan(traces),1)), 2);
                    sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
                    sg(sg==0) = 1; % avoid divide-by-zero
                    traces = 100 .* (traces-mu) ./ sg;

                    % put normalized data back to topography bended radargram:
                    radargram=NaN(M,N);
                    radargram(linear_idx)=traces;
                    traces=radargram;
                    clear radargram;

                    fprintf('x\t\t');

                    % extract only relevant samples:
                    % downsampling/cutting:
                    if downsampling==1
                        if cut_range==1
                            tsnum=1:downsampling_factor:length(t(t<=cut_time_depth));
                        else
                            tsnum=1:downsampling_factor:length(t);
                        end
                    else
                        if cut_range==1
                            tsnum=1:length(t(t<=cut_time_depth));
                        else
                            tsnum=1:length(t);
                        end
                    end
                    tnew=t(tsnum); % new depth vector (only named t for consistency!)
                    traces=traces(tsnum,:);

                    % bin data along horizontal slices:
                    dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
                    % each depth sluce for this profile has difefrent valid
                    % values (due to topo):
                    validdata = linearindex(any(~isnan(dtemp), 3)); % indices with at least one non-nan value at any depth
                    % initialize variable for all data:
                    nValid = numel(validdata);
                    nSl    = size(dtemp, 3);
                    profiledata  = zeros(nValid, nSl+1);
                    profiledata(:,1) = validdata;

                    % --- OPTIMIZATION: vectorized slice normalization ---
                    % Extract valid pixels for all slices at once: nValid x nSl
                    % Reshape dtemp to (ny*nx) x nSl, pick valid rows
                    nY = size(dtemp,1); nX = size(dtemp,2);
                    flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
                    valid_flat = flat(validdata, :);      % nValid x nSl

                    mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
                    sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
                    sg_sl(sg_sl==0) = 1;
                    profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

                    % Save profiledata:
                    save(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(i)),'.mat']),'profiledata');
                    fprintf('x\t');

                    fprintf('\t%.1f\n',toc(tstart));

                    i=i+1;
                end
            end
        end
        t=tnew;
        timesamplenum=tsnum;
    end

%%%%% TIME DOMAIN INPUT DATA
else %% if time domain input data:
    if followTopo==0
        % downsampling/cutting:
        if downsampling==1
            if cut_range==1
                timesamplenum=1:downsampling_factor:length(t(t<=cut_time_depth));
            else
                timesamplenum=1:downsampling_factor:length(t);
            end
        else
            if cut_range==1
                timesamplenum=1:length(t(t<=cut_time_depth));
            else
                timesamplenum=1:length(t);
            end
        end
    else
        % do not downsample at this point
        timesamplenum=1:numel(t);
    end
    % timesamplenum is a vector of time samples (indices) that are used for creation of
    % slices
    t=t(timesamplenum); % new time vector
    if tz_flag==1 && followTopo==1
        % convert time to depth:
        t=t./2*constV; % t is now depth in m!
    end

    % save time vector and coordtrans:
    save(fullfile(foldername,newfolder,'t.mat'),'t');
    save(fullfile(foldername,newfolder,'coordtrans.mat'),'coordtrans');

    % create inital profnum slice:
    slice_prof=NaN(size(xgrid));

    % read profile data:
    disp('-----------')
    disp(['Get data of profiles ',int2str(profnum(1)),'-',int2str(profnum(end)),'...'])
    fprintf('Profile\tData\tBinning\t\tSaved\tTime elapsed [s]\n')

    if datastore_flag==0
        for n=1:length(numbers) %  loop over profiles
            tstart=tic;
            fprintf('%d\t',profnum(n));

            % load data of this profile
            if size(m.radargrams,1)>1
                temp=m.radargrams(n,1);
            else
                temp=m.radargrams(1,n); % -> traces (all channels)
            end
            traces=temp{1};

            % coords of this profile
            ctemp=xylist(xylist(:,1)==profnum(n),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

            fprintf('x\t');

            % get profile data for relevant time samples only:
            traces=traces(timesamplenum,:);

            % ---  normalization ---
            mu  = mean(traces(:,all(~isnan(traces),1)), 2);
            sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
            sg(sg==0) = 1; % avoid divide-by-zero
            traces = 100 .* (traces-mu) ./ sg;

            % bin data:
            dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
            % valid data points:
            validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
            % initialize variable for all data:
            nValid = numel(validdata);
            nSl    = size(dtemp, 3);
            profiledata  = zeros(nValid, nSl+1);
            profiledata(:,1) = validdata;

            % --- OPTIMIZATION: vectorized slice normalization ---
            % Extract valid pixels for all slices at once: nValid x nSl
            % Reshape dtemp to (ny*nx) x nSl, pick valid rows
            nY = size(dtemp,1); nX = size(dtemp,2);
            flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
            valid_flat = flat(validdata, :);      % nValid x nSl

            mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
            sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
            sg_sl(sg_sl==0) = 1;
            profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

            fprintf('x\t\t');

            % Save profiledata:
            save(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(n)),'.mat']),'profiledata');
            fprintf('x\t');

            fprintf('\t%.1f\n',toc(tstart));
        end
    else % datastore
        i=1; % profile number overall
        j=1; % file number
        reset(radar_ds);
        while hasdata(radar_ds)
            fprintf('  -- File %d (of %d) --\n',j,numel(radar_ds.Files));
            data=read(radar_ds);
            j=j+1;

            for ii=1:numel(data.radargrams)
                % load data of this profile
                traces=data.radargrams{ii};

                tstart=tic;
                fprintf('%d\t',profnum(i));

                % coords of this profile
                ctemp=xylist(xylist(:,1)==profnum(i),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

                fprintf('x\t');

                % get profile data for relevant time samples only:
                traces=traces(timesamplenum,:);

                % ---  normalization ---
                mu  = mean(traces(:,all(~isnan(traces),1)), 2);
                sg  = std(traces(:,all(~isnan(traces),1)), 0, 2);
                sg(sg==0) = 1; % avoid divide-by-zero
                traces = 100 .* (traces-mu) ./ sg;

                % bin data:
                dtemp=bindata3_oneTracePerBin(traces,ctemp(1,:),ctemp(2,:),xrg,yrg);
                % valid data points:
                validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
                % initialize variable for all data:
                nValid = numel(validdata);
                nSl    = size(dtemp, 3);
                profiledata  = zeros(nValid, nSl+1);
                profiledata(:,1) = validdata;

                % --- OPTIMIZATION: vectorized slice normalization ---
                % Extract valid pixels for all slices at once: nValid x nSl
                % Reshape dtemp to (ny*nx) x nSl, pick valid rows
                nY = size(dtemp,1); nX = size(dtemp,2);
                flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
                valid_flat = flat(validdata, :);      % nValid x nSl

                mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
                sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
                sg_sl(sg_sl==0) = 1;
                profiledata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl; % slice normalization for this profile

                fprintf('x\t\t');

                % Save profiledata:
                save(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(i)),'.mat']),'profiledata');
                fprintf('x\t');

                fprintf('\t%.1f\n',toc(tstart));
                i=i+1;
            end
        end
    end
end

%% CREATING SAMPLE SLICES FROM DATA
disp('-----------')

for n=1:length(profnum)
    mp{n}=matfile(fullfile(foldername,newfolder,['profiledata_',int2str(profnum(n)),'.mat']));
end

% Pre-load linear indices and first-slice data for all profiles to avoid
% repeated matfile property access inside the time loop
all_linidx  = cell(length(profnum),1);
all_data_tt = cell(length(profnum),1); % will be filled per tt below
for n=1:length(profnum)
    try
        all_linidx{n} = mp{n}.profiledata(:,1); % first column is linear index in slice
    catch
        all_linidx{n} = [];
    end
end

if tz_flag==1 && followTopo==1
    fprintf('Create 3D data cube with topography\n')

    % create new 3D cube:
    maxZ=ceil(max(topo(:))*10)/10; % max topo
    minZ=floor(min(topo(:))*10)/10; % min topo
    z=fliplr(minZ-max(t):abs(t(2)-t(1)):maxZ); % new absolute depth vector (t is here in depth already), from top to bottom of cube

    nz=numel(z);

    % find 3rd dimension in cube for each topo point:
    topoIdx = interp1(z, 1:length(z), topo(:), 'nearest', 'extrap');
    topoIdx = reshape(topoIdx, size(topo));         % (ny x nx)-Matrix with index in 3rd dimension of cube

    % put slices of profiles in cube:
    msg = '';
    for n=1:length(profnum)
        % alte Zeile löschen
        fprintf(repmat('\b', 1, length(msg)));

        if isempty(all_linidx{n})
            msg = sprintf('[%-*s]', length(profnum), repmat('.', 1, n));
            fprintf(msg);
            continue;
        end

        msg = sprintf('[%-*s]', length(profnum), repmat('.', 1, n));
        fprintf(msg);

        indices{n}=zeros(numel(all_linidx{n}),ns);
        for i=1:numel(all_linidx{n}) % for each bin
            % get indices in 3D cube
            indices{n}(i,:)=topoIdx(all_linidx{n}(i)):topoIdx(all_linidx{n}(i))+ns-1;
        end

        slice_prof(all_linidx{n}) = profnum(n);
    end

    % do downsampling an cutting here:
    if downsampling==1
        if cut_range==1
            timesamplenum=1:downsampling_factor:length(z(z<=cut_time_depth));
        else
            timesamplenum=1:downsampling_factor:length(z);
        end
    else
        if cut_range==1
            timesamplenum=1:length(z(z<=cut_time_depth));
        else
            timesamplenum=1:length(z);
        end
    end
    z=z(timesamplenum);

    fprintf('\n');
    disp(['Creating ',int2str(numel(timesamplenum)),' sample slices'])
    fprintf('\nMask & Interpolation & Save\n')
    msg = '';
    % for each slice in cube:
    for i=1:numel(timesamplenum)
        % alte Zeile löschen
        fprintf(repmat('\b', 1, length(msg)));
        msg = sprintf('[%-*s]', numel(timesamplenum), repmat('.', 1, i));
        fprintf(msg);

        slice=NaN(size(xgrid));
        % fill slice with data:
        for n=1:length(profnum)
            % find indices for this slice:
            [in_r,in_c]=find(indices{n}==timesamplenum(i));
            temp=mp{n}.profiledata;
            for j=1:numel(in_c)
                slice(all_linidx{n}(in_r(j)))=temp(in_r(j),in_c(j)+1);
            end
        end

        if any(~isnan(slice(:)))
            % Mask:
            mask{i}=zeros(size(xgrid));
            mask{i}(~isnan(slice))=1;
            temp=ones(size(mask{i}));
            temp(mask{i}==1)=0;
            eucmap=chamfer_DT(temp);
            mask_interp{i}=ones(size(eucmap));
            mask_interp{i}(eucmap.*dx>radius)=0;

            % Interpolation:
            F=scatteredInterpolant(xgrid(mask{i}>0),ygrid(mask{i}>0),slice(mask{i}>0));
            slice=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
            slice(mask_interp{i}==0)=NaN;
        else
            % no data in this slice
            mask{i}=zeros(size(xgrid));
            mask_interp{i}=zeros(size(xgrid));
        end

        % save slice:
        save(fullfile(foldername,newfolder,['slice_',int2str(i),'.mat']),'slice');
    end


    

elseif (tz_flag==1 && followTopo==0) || (tz_flag==2 && followTopo==1)  % slices parallel to surface
    disp(['Creating ',int2str(numel(timesamplenum)),' sample slices'])

    fprintf('#\tData\tMask\tInterpolation\tSaved\tTime elapsed [s]\n')

    spinnerChars = {'|', '/', '-', '\'};
    for tt=1:length(timesamplenum) % for each time sample
        tstart=tic;
        fprintf('%3d\t\t',tt);

        slice=NaN(size(xgrid));

        % --- OPTIMIZATION: read column tt+1 from matfile once per profile ---
        for n=1:length(profnum)
            if isempty(all_linidx{n}), fprintf('\b%s', spinnerChars{mod(step-1, 4) + 1}); continue; end
            try
                fprintf('\b%s', spinnerChars{mod(n-1, 4) + 1});
                col_data = mp{n}.profiledata(:, tt+1);  % read single column
                slice(all_linidx{n}) = col_data;
                if tt==1
                    slice_prof(all_linidx{n}) = profnum(n);
                end
            catch
                bla=1;
            end
        end
        fprintf('\bx\t');

        % Mask: (the same for each slice)
        if tt==1
            mask{1}=zeros(size(slice));
            mask{1}(~isnan(slice))=1;
            temp=ones(size(mask{1}));
            temp(mask{1}==1)=0;
            eucmap=chamfer_DT(temp);
            mask_interp{1}=ones(size(eucmap));
            mask_interp{1}(eucmap.*dx>radius)=0;
        end
        fprintf('x\t');

        % Interpolation:
        F.Values=slice(mask{1}>0);
        slice=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
        slice(mask_interp{1}==0)=NaN;
        fprintf('x\t\t');

        % save slice:
        save(fullfile(foldername,newfolder,['slice_',int2str(tt),'.mat']),'slice');
        fprintf('x\t');

        fprintf('\t%.1f\n',toc(tstart));
    end

elseif tz_flag==2 && followTopo==0 % horizontal slices cutting through topography
    fprintf('Creating %d sampleslices:\n',numel(timesamplenum))
    fprintf('#\tData\tMask\tInterpolation\tSaved\tTime elapsed [s]\n')

    mask=cell(numel(timesamplenum),1);
    mask_interp=cell(numel(timesamplenum),1);

    spinnerChars = {'|', '/', '-', '\'};

    for tt=1:length(timesamplenum) % for each time sample
        tstart=tic;
        fprintf('%3d\t\t',tt);

        slice=NaN(size(xgrid));

        % --- OPTIMIZATION: read column tt+1 from matfile once per profile ---
        for n=1:length(profnum)
            if isempty(all_linidx{n}), fprintf('\b%s', spinnerChars{mod(n-1, 4) + 1}); continue; end
            try
                fprintf('\b%s', spinnerChars{mod(n-1, 4) + 1});
                col_data = mp{n}.profiledata(:, tt+1);  % read single column
                slice(all_linidx{n}) = col_data;
                if tt==1
                    slice_prof(all_linidx{n}) = profnum(n);
                end
            catch
                bla=1;
            end
        end
        fprintf('\bx\t');


        % Mask: (different for each slice)
        mask{tt}=zeros(size(slice));
        mask{tt}(~isnan(slice))=1;
        temp=ones(size(mask{tt}));
        temp(mask{tt}==1)=0;
        eucmap=chamfer_DT(temp);
        mask_interp{tt}=ones(size(eucmap));
        mask_interp{tt}(eucmap.*dx>radius)=0;
        fprintf('x\t');

        % Interpolation:
        if any(any(mask{tt}>0))==1 % if there is data in this slice
            F=scatteredInterpolant(xgrid(mask{tt}>0),ygrid(mask{tt}>0),slice(mask{tt}>0));
            slice=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
            slice(mask_interp{tt}==0)=NaN;
        end
        fprintf('x\t\t');

        % save slice:
        save(fullfile(foldername,newfolder,['slice_',int2str(tt),'.mat']),'slice');
        fprintf('x\t');

        fprintf('\t%.1f\n',toc(tstart));
    end

end
fprintf('\n')
disp('-----------')



%% geopng/geotif:
if save_geopng==1 || save_geotif==1
    disp('Saving sampleslices as georeferenced images...')

    clear slice;
    i=1;
    while exist(fullfile(foldername,newfolder,['slice_',int2str(i),'.mat']),'file')
        temp=load(fullfile(foldername,newfolder,['slice_',int2str(i),'.mat']));
        bla=struct2cell(temp);
        slice{i}=bla{1}.*mask_interp{min(numel(mask_interp),i)};
        i=i+1;
    end

    if removeBorder==1 % remove interpolation artifacts around area
        disp('Remove interpolation border around area...')
        for i=1:numel(slice)
            dist{i} = chamfer_DT(mask_interp{min(numel(mask_interp),i)});
        end
        dist{i+1}=chamfer_DT(mask_interp_topo); % extra mask for topo
    else
        dist=[];
    end

    % GEOPNG:
    if save_geopng==1
        disp('-- Geopng/pgw --')
        if tz_flag==1 && followTopo==0
            saveallslices_geopng(xgrid,ygrid,slice,topo_interp,t,fullfile(foldername,newfolder),colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix,1);
        elseif tz_flag==1 && followTopo==1
            saveallslices_geopng(xgrid,ygrid,slice,topo_interp,z,fullfile(foldername,newfolder),colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix,2);
        elseif tz_flag==2 && followTopo==1
            saveallslices_geopng(xgrid,ygrid,slice,topo_interp,t,fullfile(foldername,newfolder),colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix,2);
        elseif tz_flag==2 && followTopo==0
            saveallslices_geopng(xgrid,ygrid,slice,topo_interp,z_abs(timesamplenum),fullfile(foldername,newfolder),colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix,2);
        end
    end

    % GEOTIF:
    if save_geotif==1
        disp('-- Geotif --')
        if ~exist(fullfile(foldername,newfolder,'georef_tif'),'dir')
            mkdir(fullfile(foldername,newfolder,'georef_tif'));
        end
        % get UTM coords of corners auf area:
        corners=helmert([xgrid(end,1) ygrid(end,1); xgrid(end,end) ygrid(end,end); xgrid(1,1) ygrid(1,1); xgrid(1,end) ygrid(1,end)],coordtrans(:,1:2),coordtrans(:,3:4));
        % topo:
        if removeBorder==1 % remove interpolation artifacts around area
            topo_interp(dist{end}<=pix)=NaN;
        end
        if ~exist('minZ','var')
            maxZ=ceil(max(topo(:))*10)/10; % max topo
            minZ=floor(min(topo(:))*10)/10; % min topo
        end
        writeGeoTIFF(fullfile(foldername,newfolder,'georef_tif','Topo_interp.tif'), flipud(topo_interp), epsg, corners, 'jet', [minZ maxZ]);
        % sampleslices:
        if tz_flag==1 && followTopo==0
            saveallslices_geotif(xgrid,ygrid,corners,epsg,slice,t,fullfile(foldername,newfolder,'georef_tif'),colperc,sq,medianFilter,msize,removeBorder,dist,pix,1);
        elseif tz_flag==1 && followTopo==1
            saveallslices_geotif(xgrid,ygrid,corners,epsg,slice,z,fullfile(foldername,newfolder,'georef_tif'),colperc,sq,medianFilter,msize,removeBorder,dist,pix,2);
        elseif tz_flag==2 && followTopo==1
            saveallslices_geotif(xgrid,ygrid,corners,epsg,slice,t,fullfile(foldername,newfolder,'georef_tif'),colperc,sq,medianFilter,msize,removeBorder,dist,pix,2);
        elseif tz_flag==2 && followTopo==0
            saveallslices_geotif(xgrid,ygrid,corners,epsg,slice,z_abs(timesamplenum),fullfile(foldername,newfolder,'georef_tif'),colperc,sq,medianFilter,msize,removeBorder,dist,pix,2);
        end
    end
end


%%
disp('-----------')
disp('Saving additional infos...')
save(fullfile(foldername,newfolder,'slice_profilenum.mat'),'slice_prof');
save(fullfile(foldername,newfolder,'topo_interp.mat'),'topo_interp');
save(fullfile(foldername,newfolder,'mask_interp.mat'),'mask_interp');
save(fullfile(foldername,newfolder,'mask.mat'),'mask');
save(fullfile(foldername,newfolder,'mask_interp_topo.mat'),'mask_interp_topo');
if rotate_area==1
    saveas(fig1,fullfile(foldername,newfolder,'area.png'));
end

% write config file:
fid=fopen(fullfile(foldername,newfolder,'configuration.txt'),'wt');
fprintf(fid,['Bin size in m: ',num2str(dx),'\n']);
fprintf(fid,['Original number of samples: ',int2str(ns),'\n']);
fprintf(fid,['Original sampling interval: ',num2str(dt),' ns\n']);
fprintf(fid,['Original range: ',num2str((ns-1)*dt),' ns\n']);
if downsampling==1
    fprintf(fid,['Downsampling:\n  Number of samples: ',num2str(length(t)),'\n']);
    fprintf(fid,['  Sampling interval: ',num2str(t(2)-t(1)),' ns\n']);
    fprintf(fid,['  Range: ',num2str(max(t)),' ns\n']);
end
if rotate_area==1
    fprintf(fid,['Area rotated by ',int2str(rotbest),' degree.\n']);
    fprintf(fid,['Area shifted by ',num2str(shiftx),' m in x-direction and ',num2str(shifty),' m in y-direction.\n']);
end
fprintf(fid,['tz_flag: ',num2str(tz_flag,0),'\n']);
fprintf(fid,['followTopo: ',num2str(followTopo,0),'\n']);
fprintf(fid,['constV: ',num2str(constV,2),' m/ns\n']);
fclose(fid);


% set original path
path(oldpath);

% End of script.

%%--------------------------------------------------------------------------

function [xy,rotbest,shiftx,shifty,coordtrans]=rotatearea(xy)
%%% Rotate area for minimum memory
disp('Find optimum rotation angle...')
rot=-45:5:45;
area_sz = zeros(size(rot));

% Build all rotation matrices at once and apply vectorized
for r=1:length(rot)
    rmat=[cosd(rot(r)) -sind(rot(r)); sind(rot(r)) cosd(rot(r))];
    new = xy * rmat';   % (N x 2) * (2 x 2) 
    area_sz(r)=(max(new(:,1))-min(new(:,1)))*(max(new(:,2))-min(new(:,2)));
end
rotbest=rot(area_sz==min(area_sz));
disp(['Optimum rotation angle is ',num2str(rotbest),' degree. Area has been rotated. Saving coordtrans.mat for later transformation.'])

rmat=[cosd(rotbest) -sind(rotbest); sind(rotbest) cosd(rotbest)];
new = xy * rmat'; 

% move origin
shiftx=floor(min(new(:,1)));
shifty=floor(min(new(:,2)));
new(:,1)=new(:,1)-shiftx;
new(:,2)=new(:,2)-shifty;

% save coordinate pairs for later transformation
coordtrans=[new(new(:,1)==min(new(:,1)),:) xy(new(:,1)==min(new(:,1)),:);...
    new(new(:,1)==max(new(:,1)),:) xy(new(:,1)==max(new(:,1)),:);...
    new(new(:,2)==min(new(:,2)),:) xy(new(:,2)==min(new(:,2)),:);...
    new(new(:,2)==max(new(:,2)),:) xy(new(:,2)==max(new(:,2)),:)];
xy=new;
disp(['Area size is now ',int2str(round(max(new(:,1))-min(new(:,1)))),' x ',int2str(round(max(new(:,2))-min(new(:,2)))),' m (x/y).'])
end



function [xy]=apply_rotatearea(xy,rot,shiftx,shifty)
%%% Rotate area with given parameters 
rmat=[cosd(rot) -sind(rot); sind(rot) cosd(rot)];
xy = xy * rmat';
xy(:,1)=xy(:,1)-shiftx;
xy(:,2)=xy(:,2)-shifty;
end



function saveallslices_geotif(xgrid,ygrid,corners,epsg,slice,t,foldername,colperc,sq,medianFilter,msize,removeBorder,dist,pix,tz)

    for numtsl=1:length(slice)
        disp(['   ',int2str(numtsl),'/',int2str(length(slice))])
    
        if medianFilter==1
            slice{numtsl}=medianfilt2(slice{numtsl},[msize msize]);
        end
    
        if removeBorder==1 % remove interpolation artifacts around area
            slice{numtsl}(dist{numtsl}<=pix)=NaN;
        end
    
        % prepare colorscale:
        if sq==1
            cdata=sqrt(slice{numtsl});
        else
            cdata=slice{numtsl};
        end
        cmin=min(cdata(:));
        cmax=max(cdata(:));
    
        tslname = fullfile(foldername,make_fname(numtsl,'.tif',t,tz));
    
        if colperc==0
            writeGeoTIFF(tslname, flipud(cdata), epsg, corners, flipud(gray(16)), [cmin cmax]);
        else
            coldata=sort(cdata(~isnan(cdata)));
            if ~isempty(coldata) && length(coldata)>2
                cmin=coldata(round(length(coldata)/100*colperc));
                cmax=coldata(end-round(length(coldata)/100*colperc));
            end
            writeGeoTIFF(tslname, flipud(cdata), epsg, corners, flipud(gray(16)), [cmin cmax]);
        end
    end
end



function saveallslices_geopng(xgrid,ygrid,slice,topo,t,pfad,colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix,tz)
dx=abs(xgrid(1,1)-xgrid(1,2));
% topo:
if removeBorder==1 % remove interpolation artifacts around area
    topo(dist{end}<=pix)=NaN;
end
cdata=topo;
cmin=min(cdata(:));
cmax=max(cdata(:));

if ~exist(fullfile(pfad,'georef_png'),'dir')
    mkdir(fullfile(pfad,'georef_png'));
end
cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
cdata(isnan(cdata))=0;  % set nan to 0
imwrite(flipud(cdata).*256,parula(256),fullfile(pfad,'georef_png','Topo_interp.png'),'Transparency',0);
fid=fopen(fullfile(pfad,'georef_png','Topo_interp.pgw'),'wt');
fprintf(fid,[num2str(dx),'\n0\n0\n',num2str(-dx),'\n',num2str(min(xgrid(:))),'\n',num2str(max(ygrid(:)))]);
fclose(fid);

% slices
for numtsl=1:length(slice)
    disp(['   ',int2str(numtsl),'/',int2str(length(slice))])

    if medianFilter==1
        slice{numtsl}=medianfilt2(slice{numtsl},[msize msize]);
    end

    if removeBorder==1 % remove interpolation artifacts around area
        slice{numtsl}(dist{numtsl}<=pix)=NaN;
    end

    % Georeferenced png:
    if sq==1
        cdata=sqrt(slice{numtsl});
    else
        cdata=slice{numtsl};
    end

    cmin=min(cdata(:));
    cmax=max(cdata(:));

    tslname = fullfile(pfad,'georef_png',make_fname(numtsl,'.png',t,tz));

    if colperc==0
        cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
        cdata(isnan(cdata))=0;  % set nan to 0
        imwrite(flipud(cdata).*256,flipud(gray(256)),tslname,'Transparency',0);
    else
        coldata=sort(cdata(~isnan(cdata)));
        if ~isempty(coldata) && length(coldata)>2
            cmin=coldata(round(length(coldata)/100*colperc));
            cmax=coldata(end-round(length(coldata)/100*colperc));
            range=cmax-cmin;
            cdata=(cdata-cmin)/range;
            cdata(cdata<=0)=0;
            cdata(cdata>=1)=1;
            m=ones(size(cdata));
            m(isnan(cdata))=0;
        end
        im=cdata.*256;
        im(im<=2)=2;
        im(isnan(cdata))=0;  % set nan to 0
        imwrite(flipud(im),flipud(gray(256)),tslname,'Transparency',0);
    end


    % write pngw
    fname = make_fname(numtsl,'.pgw',t,tz);
    if ~exist('coordtrans','var')    % local
        fid=fopen(fullfile(pfad,'georef_png',fname),'wt');
        fprintf(fid,[num2str(dx),'\n0\n0\n',num2str(-dx),'\n',num2str(min(xgrid(:))),'\n',num2str(max(ygrid(:)))]);
        fclose(fid);
    else % global
        write_geoPNGW(xgrid,ygrid,coordtrans,fullfile(pfad,'georef_png',fname));
    end

end
end

%%
function fnameStr = make_fname(numtsl,extension,t,tz)
    % create filename
    if tz==2 % depth
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'m',extension];
    elseif tz==1 % time
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'ns',extension];
    end
end

%%
function []=write_geoPNGW(x,y,coordtrans,filename)

dx=abs(x(1,2)-x(1,1));
dy=abs(y(2,1)-y(1,1));

% determine global coords of upper left and upper right pixel
pix_ol=helmert([min(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
pix_or=helmert([max(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
alpha=atand(abs(pix_ol(2)-pix_or(2))/abs(pix_ol(1)-pix_or(1))); % angle against west
if pix_ol(2)>pix_or(2)
    alpha=-alpha;
end
% determine pixel-lengths in all directions
A=dx*cosd(alpha);
D=dx*sind(alpha);
E=dy*cosd(alpha);
B=dy*sind(alpha);
    
% write pngw
fid=fopen(filename,'wt');
fprintf(fid,[num2str(A),'\n',num2str(D),'\n',num2str(B),'\n',num2str(-E),'\n',num2str(pix_ol(1)),'\n',num2str(pix_ol(2))]);
fclose(fid);
end


function y = medfilt1_own(x, k)
    half = floor(k/2);
    % Ränder mit erstem/letztem Wert padden
    x_pad = [repmat(x(1), half, 1); x; repmat(x(end), half, 1)];
    n = length(x);
    y = zeros(size(x));
    
    for i = 1:n
        y(i) = median(x_pad(i : i+k-1));
    end
end