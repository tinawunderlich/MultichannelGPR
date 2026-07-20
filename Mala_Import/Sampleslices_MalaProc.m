% Script for reading processed Mala-datafiles in profiles2mat/proc and
% binning them onto a rectangular grid (for each channel individually to
% get balanced channel energies and less stripes in timeslices)
%
% ONLY WORKS IN TIME DOMAIN! IF YOU WANT TO USE DEPTH DOMAIN OR TOPOGRAPHY
% SETTINGS, PLEASE EXPORT PROFILES AS radargrams.mat AND USE
% SAMPLESLICES_FROMRADARGRAMS.m
%
% Dr. Tina Wunderlich, CAU Kiel 2025-2026, tina.wunderlich@ifg.uni-kiel.de
% OPTIMIZED VERSION (by help of Claude.ai free version, manually checked!)
%
% requires MATLAB-files in following folders (path will be temporarily
% set):  Subfunctions

clear all
close all
clc

% Bin size of grid
dx=0.05; % [m]

radius=0.3; % radius in m for valid interpolation (-> mask)

% virtual channels for interpolation:
virt_chan_num=3; % if =0: only virtual channels between real channels,
                    % if e.g. =3: extrapolate also to 3 channels before and
                    % after last real channel to the sides

% Automatic rotation of measurement area for minimum memory size
rotate_area=1;  % 1=yes (recommended), 0=no


% Downsampling of data
downsampling=1; % if =1: yes (and use following settings)
downsampling_factor=20; % only take each downsampling-factor sample (e.g. only take every 2nd sample)
% cutting of range
cut_range=0; % if =1:yes
cut_time=20; % choose time for cutting [ns]

% save Sampleslices as geopng?
save_geopng=1; % 1=yes
colperc=3; % Colorscale clipping in percent (if =0: autoscale min-max)
removeBorder=0; % =1: remove border artifacts from interpolation, =0: leave as it is
pix=6; % if removeBorder==1: how many pixels are removed from border around area
medianFilter=1; % do you want to apply a 2D-median filter (1=yes, 0=no)
msize=3; % filter size in pixel
% use squareroot of amplitudes for visualization?
sq=0; % 1=yes, 0=no


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
    if exist('.temp.temp','file') % read last opened folder from temp.temp
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
name_withoutGPS=name; % name of files when not using GPS
for i=2:length(tempname)-1
    name=[name,'_',tempname{i}];
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));



%% check if correctly processed data available:
if ~exist(fullfile(foldername,'profiles2mat','proc'),'dir')
    disp('No processed data for binning found (no data in profiles2mat/proc). Stopping!')
    return;
end
% load profileinfo
temp=load(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'));
profileinfo=temp.profileinfo; % profilnumber, dt, ns, channels, numtraces-per-channel
t=0:profileinfo(1,2):profileinfo(1,2)*(profileinfo(1,3)-1);
dt=profileinfo(1,2);
numbers=profileinfo(:,1); % profile numbers
disp(['Found processed data of ',int2str(length(numbers)),' profiles (Min: ',int2str(min(numbers)),'; Max: ',int2str(max(numbers)),').'])


%% Size of area
% read info files for coordinates
disp('Reading coordinates for determining area size...')
xylist=zeros(sum(profileinfo(:,4).*profileinfo(:,5)),7);
anz=0;
channel_num=cell(length(numbers),1);
numtraces=zeros(length(numbers),1);
for i=1:length(numbers)
    % read info for this profile:
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat'])); % variable info
    elseif exist(fullfile(foldername,'profiles2mat',[name_withoutGPS,'_',int2str(numbers(i)),'_info_proc.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name_withoutGPS,'_',int2str(numbers(i)),'_info_proc.mat']));
        name=name_withoutGPS; % set correct name
    else
        continue;
    end

    nrows = length(info(4,:));
    idx = anz+1:anz+nrows;
    xylist(idx,:) = [repmat(numbers(i),nrows,1) info(4,:)' info(5,:)' info(6,:)' info(3,:)' (1:nrows)' info(2,:)'];  % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
    anz = anz + nrows;

    if i==1
        channels=unique(xylist(xylist(:,1)==numbers(i),5));
        channels=channels(channels~=0);
        ns=profileinfo(1,3);
    end
    numtraces(i)=max(info(2,:)); % number of traces in this profile for each channel
end


% optional: rotate area
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

newfolder='SampleSlices_time_parallel2surface';

if ~exist(fullfile(foldername,newfolder),'dir')
    mkdir(fullfile(foldername,newfolder))
end
% save x/y-grids:
save(fullfile(foldername,newfolder,'xgrid.mat'),'xgrid','-v7.3');
save(fullfile(foldername,newfolder,'ygrid.mat'),'ygrid','-v7.3');

% bin edges:
xrg=min(xylist(:,2))-dx/2:dx:max(xylist(:,2))+dx/2;
yrg=min(xylist(:,3))-dx/2:dx:max(xylist(:,3))+dx/2;

% make topography bins:
topo=bindata2(xylist(:,4),xylist(:,2),xylist(:,3),xrg,yrg);

% downsampling/cutting:
if downsampling==1
    if cut_range==1
        timesamplenum=1:downsampling_factor:length(t(t<=cut_time));
    else
        timesamplenum=1:downsampling_factor:length(t);
    end
else
    if cut_range==1
        timesamplenum=1:length(t(t<=cut_time));
    else
        timesamplenum=1:length(t);
    end
end
% timesamplenum is a vector of time samples (indices) that are used for creation of
% slices
t=t(timesamplenum); % new time vector

% save time vector and coordtrans:
save(fullfile(foldername,newfolder,'t.mat'),'t','-v7.3');
save(fullfile(foldername,newfolder,'coordtrans.mat'),'coordtrans','-v7.3');

% create inital profnum & channum slices:
slice_chan=NaN(size(xgrid));
slice_prof=NaN(size(xgrid));

% read profile data:
disp('-----------')
disp(['Get data of profiles ',int2str(numbers(1)),'-',int2str(numbers(end)),'...'])
fprintf('Profile\tData\tVirtual channels\tBinning\t\t\tSaved\tTime elapsed [s]\n')

for n=1:length(numbers) %  loop over profiles
    tstart=tic;
    fprintf('%d\t',numbers(n));

    % load data of this profile
    load(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(numbers(n)),'.mat'])); % -> traces (all channels)

    % coords of this profile
    ctemp=xylist(xylist(:,1)==numbers(n),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

    fprintf('x\t');

    % get profile data for relevant time samples only:
    traces=traces(timesamplenum,:);

    % --- OPTIMIZATION: vectorized per-channel normalization ---
    for ch=1:length(channels)
        mask_ch = (ctemp(3,:)==channels(ch));
        blk = traces(:, mask_ch); % only data of one channel
        mu  = mean(blk(:,all(~isnan(blk),1)), 2);
        sg  = std(blk(:,all(~isnan(blk),1)), 0, 2);
        sg(sg==0) = 1; % avoid divide-by-zero
        traces(:, mask_ch) = 100 .* (blk-mu) ./ sg;
    end

    % --- OPTIMIZATION: vectorized virtual channel creation ---
    trnum=unique(ctemp(4,:)); % all trace numbers
    nTr = length(trnum);
    virtchan=[min(channels)-virt_chan_num*0.5:0.5:min(channels)-0.01, ...
              min(channels)+0.5:max(channels)-0.5, ...
              fliplr(max(channels)+0.5*virt_chan_num:-0.5:max(channels)+0.01)];
    nVirt = numel(virtchan);
    nSamp = numel(timesamplenum);

    virtchan_data = zeros(nSamp, nTr * nVirt);
    virtchan_xy   = zeros(4, nTr * nVirt);

    % Build lookup: for each trace number, gather channel data matrix
    % interp1 over channels dimension — vectorize across all traces at once
    % by building a 3-D array: (samples x channels x traces)
    nCh = length(channels);
    % Check all trace numbers have same channel count (typical for Mala)
    % and build 3-D arrays for batch interp
    allSame = all(accumarray(ctemp(4,:)', 1) == nCh);

    if allSame && nTr > 1
        % Fast path: reshape into 3D and interp once
        % Sort ctemp by (trnum, channel) for reliable reshape
        [~, sortIdx] = sortrows(ctemp([4,3],:)');
        traces_sorted = traces(:, sortIdx);   % nSamp x (nTr*nCh)
        xy_sorted     = ctemp(:, sortIdx);

        % Reshape: nSamp x nCh x nTr  and  4 x nCh x nTr
        D3 = reshape(traces_sorted, nSamp, nCh, nTr);   % data
        X3 = reshape(xy_sorted,     4,     nCh, nTr);   % coords

        % channels per trace (same for all): use channels vector
        c_vec = squeeze(X3(3, :, 1));  % nCh channel IDs for first trace

        % interp over channel dimension for all traces simultaneously
        % Result: nSamp x nVirt x nTr
        Dv = zeros(nSamp, nVirt, nTr);
        Xv = zeros(2,     nVirt, nTr);
        for vv = 1:nVirt
            % linear interp weight
            xi = virtchan(vv);
            % find bracket
            lo = find(c_vec <= xi, 1, 'last');
            hi = find(c_vec >= xi, 1, 'first');
            if isempty(lo) || isempty(hi)
                % extrapolate: use two nearest
                if isempty(lo), lo=1; hi=2; end
                if isempty(hi), hi=nCh; lo=nCh-1; end
            end
            if lo==hi
                w=1; % exact match
                Dv(:, vv, :) = D3(:, lo, :);
                Xv(:, vv, :) = X3(1:2, lo, :);
            else
                w = (xi - c_vec(lo)) / (c_vec(hi) - c_vec(lo));
                Dv(:, vv, :) = (1-w)*D3(:, lo, :) + w*D3(:, hi, :);
                Xv(:, vv, :) = (1-w)*X3(1:2, lo, :) + w*X3(1:2, hi, :);
            end
        end
        % Reshape back to 2D
        virtchan_data = reshape(Dv, nSamp, nVirt*nTr);
        virtchan_xy(1:2,:) = reshape(Xv, 2, nVirt*nTr);
        % channel and trace labels
        virtchan_xy(3,:) = repmat(virtchan(:), nTr, 1)';
        virtchan_xy(4,:) = repelem(trnum(:)', nVirt);
    else
        % Fallback: original per-trace loop (for unequal channel counts)
        for tr=1:nTr
            if ~mod(tr,round(length(trnum)/10))
                fprintf('.');
            end
            d=traces(:,ctemp(4,:)==trnum(tr));
            xytemp=ctemp(1:2,ctemp(4,:)==trnum(tr));
            c=ctemp(3,ctemp(4,:)==trnum(tr));
            virtchan_data(:,nVirt*(tr-1)+1:nVirt*tr)=interp1(c,d',virtchan,'linear','extrap')';
            virtchan_xy(1:2,nVirt*(tr-1)+1:nVirt*tr)=interp1(c,xytemp',virtchan,'linear','extrap')';
            virtchan_xy(3,nVirt*(tr-1)+1:nVirt*tr)=virtchan;
            virtchan_xy(4,nVirt*(tr-1)+1:nVirt*tr)=trnum(tr);
        end
    end
    fprintf('x\t\t');

    % add to real channel data:
    channels_all=[virtchan(:); channels(:)]; %  virtual channels first
    traces=[virtchan_data traces];
    ctemp=[virtchan_xy ctemp];

    % bin data:
    profiledata=[];
    chan_prof=[];
    for ch=1:length(channels_all)
        if ~mod(ch,2)
            fprintf('.');
        end
        chanInd=ctemp(3,:)==channels_all(ch);
        dtemp=bindata3_oneTracePerBin(traces(:,chanInd),ctemp(1,chanInd),ctemp(2,chanInd),xrg,yrg);
        % valid data points:
        validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
        % initialize variable for all data:
        nValid = numel(validdata);
        nSl    = size(dtemp, 3);
        pdata  = zeros(nValid, nSl+1);
        pdata(:,1) = validdata;

        % --- OPTIMIZATION: vectorized slice normalization ---
        % Extract valid pixels for all slices at once: nValid x nSl
        % Reshape dtemp to (ny*nx) x nSl, pick valid rows
        nY = size(dtemp,1); nX = size(dtemp,2);
        flat = reshape(dtemp, nY*nX, nSl);   % (ny*nx) x nSl
        valid_flat = flat(validdata, :);      % nValid x nSl

        mu_sl  = mean(valid_flat, 1, 'omitnan');   % 1 x nSl
        sg_sl  = std(valid_flat,  0, 1, 'omitnan'); % 1 x nSl
        sg_sl(sg_sl==0) = 1;
        pdata(:, 2:end) = 100 .* bsxfun(@minus, valid_flat, mu_sl) ./ sg_sl;

        profiledata=[profiledata; pdata];
        chan_prof=[chan_prof; zeros(nValid,1)+channels_all(ch)];
    end
    chan_prof(:,2)=numbers(n); % profile number
    fprintf('x\t');

    % Save profiledata:
    save(fullfile(foldername,newfolder,['profiledata_',int2str(numbers(n)),'.mat']),'profiledata','-v7.3');
    save(fullfile(foldername,newfolder,['chan_prof_',int2str(numbers(n)),'.mat']),'chan_prof','-v7.3');
    fprintf('x\t');

    fprintf('\t%.1f\n',toc(tstart));    
end

%%
disp('-----------')
disp(['Creating ',int2str(numel(timesamplenum)),' sample slices'])
fprintf('#\tData\t\t\tSaved\tMask\tTime elapsed [s]\n')

for n=1:length(numbers)
    m{n}=matfile(fullfile(foldername,newfolder,['profiledata_',int2str(numbers(n)),'.mat']));
    mcp(n)=load(fullfile(foldername,newfolder,['chan_prof_',int2str(numbers(n)),'.mat']));
end

% Pre-load linear indices and first-slice data for all profiles to avoid
% repeated matfile property access inside the time loop
disp('Pre-loading profile indices...')
all_linidx  = cell(length(numbers),1);
all_data_tt = cell(length(numbers),1); % will be filled per tt below
for n=1:length(numbers)
    try
        all_linidx{n} = m{n}.profiledata(:,1); % first column is linear index in slice
    catch
        all_linidx{n} = [];
    end
end

for tt=1:length(timesamplenum) % for each time sample
    tstart=tic;
    fprintf('%d\t',tt);

    slice=NaN(size(xgrid));

    % --- OPTIMIZATION: read column tt+1 from matfile once per profile ---
    for n=1:length(numbers)
        if isempty(all_linidx{n}), fprintf('.'); continue; end
        try
            fprintf('.');
            col_data = m{n}.profiledata(:, tt+1);  % read single column
            slice(all_linidx{n}) = col_data;
            if tt==1
                slice_chan(all_linidx{n}) = mcp(n).chan_prof(:,1);
                slice_prof(all_linidx{n}) = mcp(n).chan_prof(:,2);
            end
        catch
            bla=1;
        end
    end
    fprintf('x\t');

    % save slice:
    save(fullfile(foldername,newfolder,['slice_',int2str(tt),'.mat']),'slice','-v7.3');
    fprintf('x\t');

    fprintf('\t%.1f\n',toc(tstart));
end




disp('-----------')
disp('Creating mask')
% Mask, the same for all
mask=zeros(size(slice));
mask(~isnan(slice))=1;

temp=ones(size(mask));
temp(mask==1)=0;
eucmap=chamfer_DT(temp);
mask_interp=ones(size(eucmap));
mask_interp(eucmap.*dx>radius)=0;

disp('-----------')
disp('Interpolate topography...')
F=scatteredInterpolant(xgrid(mask>0),ygrid(mask>0),topo(mask>0));
topo_interp=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
topo_interp(mask_interp==0)=NaN;

% geopng:
if save_geopng==1
    disp('Saving sampleslices as geopng...')

    clear slice;
    i=1;
    while exist(fullfile(foldername,newfolder,['slice_',int2str(i),'.mat']),'file')
        temp=load(fullfile(foldername,newfolder,['slice_',int2str(i),'.mat']));
        bla=struct2cell(temp);
        slice{i}=bla{1}.*mask_interp;
        i=i+1;
    end

    if removeBorder==1 % remove interpolation artifacts around area
        disp('Remove interpolation border around area...')
        dist = chamfer_DT(mask_interp);
    else
        dist=[];
    end

    saveallslices(xgrid,ygrid,slice,topo_interp,t,fullfile(foldername,newfolder),colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix);
end



disp('-----------')
disp('Saving additional infos...')
save(fullfile(foldername,newfolder,'slice_channelnum.mat'),'slice_chan','-v7.3');
save(fullfile(foldername,newfolder,'slice_profilenum.mat'),'slice_prof','-v7.3');
save(fullfile(foldername,newfolder,'topo_interp.mat'),'topo_interp','-v7.3');
save(fullfile(foldername,newfolder,'mask_interp.mat'),'mask_interp','-v7.3');
save(fullfile(foldername,newfolder,'mask.mat'),'mask','-v7.3');
if rotate_area==1
    saveas(fig1,fullfile(foldername,newfolder,'area.png'));
end

% write config file:
fid=fopen(fullfile(foldername,newfolder,'configuration.txt'),'wt');
fprintf(fid,['Number of channels: ',int2str(length(channels)),'\n']);
fprintf(fid,['Bin size in m: ',num2str(dx),'\n']);
fprintf(fid,['Original number of samples: ',int2str(ns),'\n']);
fprintf(fid,['Original sampling interval: ',num2str(dt),' ns\n']);
fprintf(fid,['Original range: ',num2str((ns-1)*profileinfo(1,2)),' ns\n']);
if downsampling==1
    fprintf(fid,['Downsampling:\n  Number of samples: ',num2str(length(t)),'\n']);
    fprintf(fid,['  Sampling interval: ',num2str(t(2)-t(1)),' ns\n']);
    fprintf(fid,['  Range: ',num2str(max(t)),' ns\n']);
end
if rotate_area==1
    fprintf(fid,['Area rotated by ',int2str(rotbest),' degree.\n']);
    fprintf(fid,['Area shifted by ',num2str(shiftx),' m in x-direction and ',num2str(shifty),' m in y-direction.\n']);
end
fclose(fid);


% set original path
path(oldpath);

% End of script.

%--------------------------------------------------------------------------

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

function saveallslices(xgrid,ygrid,slice,topo,t,pfad,colperc,coordtrans,sq,medianFilter,msize,removeBorder,dist,pix)
dx=abs(xgrid(1,1)-xgrid(1,2));
for numtsl=1:length(slice)
    disp(['   ',int2str(numtsl),'/',int2str(length(slice))])

    if medianFilter==1
        slice{numtsl}=medianfilt2(slice{numtsl},[msize msize]);
    end

    if removeBorder==1 % remove interpolation artifacts around area
        slice{numtsl}(dist<=pix)=NaN;
    end

    % Georeferenced png:
    if sq==1
        cdata=sqrt(slice{numtsl});
    else
        cdata=slice{numtsl};
    end

    cmin=min(cdata(:));
    cmax=max(cdata(:));

    if ~exist(fullfile(pfad,'georef'),'dir')
        mkdir(fullfile(pfad,'georef'));
    end

    tslname = fullfile(pfad,'georef',make_fname(numtsl,'.png',t));

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
    fname = make_fname(numtsl,'.pgw',t);
    if ~exist('coordtrans','var')    % local
        fid=fopen(fullfile(pfad,fname),'wt');
        fprintf(fid,[num2str(dx),'\n0\n0\n',num2str(-dx),'\n',num2str(min(xgrid(:))),'\n',num2str(max(ygrid(:)))]);
        fclose(fid);
    else % global
        write_geoPNGW(xgrid,ygrid,coordtrans,fullfile(pfad,'georef',fname));
    end

end
end

%%
function fnameStr = make_fname(numtsl,extension,t)
    % create filename
    fnameStr = ['Tsl','_',num2str(numtsl,'%2d'),'_t',num2str(t(numtsl),2),'ns',extension];
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