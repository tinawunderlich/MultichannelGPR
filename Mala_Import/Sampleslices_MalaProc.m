% Script for reading processed Mala-datafiles in profiles2mat/proc and
% binning them onto a rectangular grid (for each channel individually to
% get balanced channel energies and less stripes in timeslices)
%
% Dr. Tina Wunderlich, CAU Kiel 2025, tina.wunderlich@ifg.uni-kiel.de
%
% requires MATLAB-files in following folders (path will be temporarily
% set):  Subfunctions


clear all
close all
clc

% Bin size of grid
dx=0.05; % [m]

radius=0.1; % radius in m for valid interpolation (-> mask)

% virtual channels for interpolation:
virt_chan_num=1; % if =0: only virtual channels between real channels,
                    % if e.g. =3: extrapolate also to 3 channels before and
                    % after last real channel to the sides

% Automatic rotation of measurement area for minimum memory size
rotate_area=1;  % 1=yes (recommended), 0=no

% time or depth?
tz_flag=1;  % 1: time -> [ns], 2: depth -> [m]

% if depth: follow Topography or horizontal slices?
followTopo=0; % =1: yes; =0: horizontal slices

% Downsampling of data
downsampling=1; % if =1: yes (and use following settings)
downsampling_factor=4; % only take each downsampling-factor sample (e.g. only take every 2nd sample)
% cutting of range
cut_range=1; % if =1:yes
cut_time=50; % choose time for cutting [ns]

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
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat'])); % variable info
        if i==1
            xylist(1:profileinfo(i,4)*profileinfo(i,5),:)=[zeros(length(info(4,:)),1)+numbers(i) info(4,:)' info(5,:)' info(6,:)' info(3,:)' [1:profileinfo(i,4)*profileinfo(i,5)]' info(2,:)']; % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
            anz=anz+profileinfo(i,4)*profileinfo(i,5);
        else
            xylist(anz+1:anz+length(info(4,:)),:)=[zeros(length(info(4,:)),1)+numbers(i) info(4,:)' info(5,:)' info(6,:)' info(3,:)' [1:length(info(3,:))]' info(2,:)']; % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
            anz=anz+length(info(4,:));
        end
    elseif exist(fullfile(foldername,'profiles2mat',[name_withoutGPS,'_',int2str(numbers(i)),'_info_proc.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name_withoutGPS,'_',int2str(numbers(i)),'_info_proc.mat']));
        name=name_withoutGPS; % set correct name
        if i==1
            xylist(1:profileinfo(i,4)*profileinfo(i,5),:)=[zeros(length(info(4,:)),1)+numbers(i) info(4,:)' info(5,:)' info(6,:)' info(3,:)' [1:profileinfo(i,4)*profileinfo(i,5)]' info(2,:)']; % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
            anz=anz+profileinfo(i,4)*profileinfo(i,5);
        else
            xylist(anz+1:anz+length(info(4,:)),:)=[zeros(length(info(4,:)),1)+numbers(i) info(4,:)' info(5,:)' info(6,:)' info(3,:)' [1:length(info(3,:))]' info(2,:)']; % Number, x, y, z, channel of profile, tracenumber in profile, tracenumber in channel
            anz=anz+length(info(4,:));
        end
    end

    if i==1
        channels=unique(xylist(xylist(:,1)==numbers(i),5));
        channels=channels(channels~=0);
        ns=profileinfo(1,3);
    end
    numtraces(i)=max(info(2,:)); % number of traces in this profile for each channel
end


% optional: rotate area
if rotate_area==1
    [xylist(:,2:3),rotbest,shiftx,shifty,coordtrans]=rotatearea(xylist(:,2:3)); % disp() included in function

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

if ~exist(fullfile(foldername,'SampleSlices'),'dir')
    mkdir(fullfile(foldername,'SampleSlices'))
end
% save x/y-grids:
save(fullfile(foldername,'SampleSlices','xgrid.mat'),'xgrid','-v7.3');
save(fullfile(foldername,'SampleSlices','ygrid.mat'),'ygrid','-v7.3');

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
save(fullfile(foldername,'SampleSlices','t.mat'),'t','-v7.3');
save(fullfile(foldername,'SampleSlices','coordtrans.mat'),'coordtrans','-v7.3');

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
    load(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(numbers(n)),'.mat'])); % traces

    % coords of this profile
    ctemp=xylist(xylist(:,1)==numbers(n),[2:3 5 7])'; % x/y-coordinates & channel number & trace number in channel

    fprintf('x\t');

    % get profile data for relevant time samples only:
    traces=traces(timesamplenum,:);
    % normalize that each channel has a mean of 100 for each time sample:
    for ch=1:channels
        traces(:,ctemp(3,:)==channels(ch))=100.*(traces(:,ctemp(3,:)==channels(ch))-mean(traces(:,ctemp(3,:)==channels(ch)),2))./std(traces(:,ctemp(3,:)==channels(ch)),0,2);
    end

    % create virtual channels in between real channels:
    trnum=unique(ctemp(4,:)); % all trace numbers
    virtchan=[min(channels)-virt_chan_num*0.5:0.5:min(channels)-0.01 min(channels)+0.5:max(channels)-0.5 fliplr(max(channels)+0.5*virt_chan_num:-0.5:max(channels)+0.01)]; % virtual channel numbers
    virtchan_data=zeros(numel(timesamplenum),max(trnum)*numel(virtchan)); % initialization
    virtchan_xy=zeros(4,max(trnum)*numel(virtchan));
    for tr=1:length(trnum)
        if ~mod(tr,round(length(trnum)/10))
            fprintf('.');
        end
        d=traces(:,ctemp(4,:)==trnum(tr)); % data for this trace number of all channels
        xytemp=ctemp(1:2,ctemp(4,:)==trnum(tr)); % xy coordinates for these data
        c=ctemp(3,ctemp(4,:)==trnum(tr)); % corresponding channel numbers
        % interpolate data and coordinates for this trace number:
        virtchan_data(:,numel(virtchan)*(tr-1)+1:numel(virtchan)*tr)=interp1(c,d',virtchan,'linear','extrap')'; % data of virtual channels
        virtchan_xy(1:2,numel(virtchan)*(tr-1)+1:numel(virtchan)*tr)=interp1(c,xytemp',virtchan,'linear','extrap')'; % xy of virtual channels
        % set trace number and channel number
        virtchan_xy(3,numel(virtchan)*(tr-1)+1:numel(virtchan)*tr)=virtchan;
        virtchan_xy(4,numel(virtchan)*(tr-1)+1:numel(virtchan)*tr)=trnum(tr);
    end
    % add to real channel data:
    channels_all=[virtchan(:); channels(:)]; %  virtual channels first
    traces=[virtchan_data traces];
    ctemp=[virtchan_xy ctemp];

    fprintf('x\t\t');

    % bin data:
    profiledata=[];
    chan_prof=[];
    for ch=1:length(channels_all)
        if ~mod(ch,2)
            fprintf('.');
        end
        dtemp=bindata3_oneTracePerBin(traces(:,ctemp(3,:)==channels_all(ch)),ctemp(1,ctemp(3,:)==channels_all(ch)),ctemp(2,ctemp(3,:)==channels_all(ch)),xrg,yrg);
        % valid data points:
        validdata=linearindex(~isnan(dtemp(:,:,1))); % indices with data
        % initialize variable for all data:
        pdata=zeros(numel(validdata),numel(timesamplenum));
        pdata(:,1)=validdata; 
        % normalize data that mean is 100 for each channel for each slice:
        for sl=1:size(dtemp,3)
            ddtemp=100.*(dtemp(:,:,sl)-mean(dtemp(:,:,sl),'all','omitnan'))/std(dtemp(:,:,sl),0,'all','omitnan');
            pdata(:,sl+1)=ddtemp(validdata);
        end
        profiledata=[profiledata; pdata];
        chan_prof=[chan_prof; zeros(size(validdata))+channels_all(ch)]; % channel number / profile number according to rows in profiledata
    end
    chan_prof(:,2)=numbers(n); % profile number
    fprintf('x\t');

    % Save profiledata:
    save(fullfile(foldername,'SampleSlices',['profiledata_',int2str(numbers(n)),'.mat']),'profiledata','-v7.3');
    save(fullfile(foldername,'SampleSlices',['chan_prof_',int2str(numbers(n)),'.mat']),'chan_prof','-v7.3');
    fprintf('x\t');

    fprintf('\t%.1f\n',toc(tstart));    
end

disp('-----------')
disp(['Creating ',int2str(numel(timesamplenum)),' sample slices'])
fprintf('#\tData\t\t\tSaved\tMask\tTime elapsed [s]\n')
for n=1:length(numbers)
    m{n}=matfile(fullfile(foldername,'SampleSlices',['profiledata_',int2str(numbers(n)),'.mat']));
    mcp(n)=load(fullfile(foldername,'SampleSlices',['chan_prof_',int2str(numbers(n)),'.mat']));
end

for tt=1:length(timesamplenum) % for each time sample
    tstart=tic;
    fprintf('%d\t',tt);

    slice=NaN(size(xgrid));
    
    % fill with data (in case that one bin has a virtual and real channel
    % data: real channel data overwrites the virtual channel)
    for n=1:length(numbers)
        fprintf('.');
        slice(m{n}.profiledata(:,1))=m{n}.profiledata(:,tt+1);
        if tt==1
            slice_chan(m{n}.profiledata(:,1))=mcp(n).chan_prof(:,1); % channel number
            slice_prof(m{n}.profiledata(:,1))=mcp(n).chan_prof(:,2); % profile number
        end
    end
    fprintf('x\t');

    % save slice:
    save(fullfile(foldername,'SampleSlices',['slice_',int2str(tt),'.mat']),'slice','-v7.3');
    fprintf('x\t');

    fprintf('\t%.1f\n',toc(tstart));
end


disp('-----------')
disp('Creating mask')
for tt=1:length(timesamplenum) % for each time sample
    % Mask:
    if tt==1 && tz_flag==1 % TIMEslice
        mask=zeros(size(slice));
        mask(~isnan(slice))=1;

        temp=ones(size(mask));
        temp(mask==1)=0;
        eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
        mask_interp=ones(size(eucmap)); % initialize new grid
        mask_interp(eucmap.*dx>radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
    elseif tz_flag==2 % DEPTHslices -> combine masks into maximum mask
        if tt==1
            mask=zeros(size(slice)); % initialize mask
        end
        mask(~isnan(slice))=1;

        temp=ones(size(mask));
        temp(mask>=1)=0;
        eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
        mask_interp=ones(size(eucmap)); % initialize new grid
        mask_interp(eucmap.*dx>radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
    end
end

disp('-----------')
disp('Interpolate topography...')
% interpolate topography:
F=scatteredInterpolant(xgrid(mask>0),ygrid(mask>0),topo(mask>0));
topo_interp=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
topo_interp(mask_interp==0)=NaN; % apply mask to topo

disp('-----------')
disp('Saving additional infos...')
% save slice_chan/slice_prof:
save(fullfile(foldername,'SampleSlices','slice_channelnum.mat'),'slice_chan','-v7.3');
save(fullfile(foldername,'SampleSlices','slice_profilenum.mat'),'slice_prof','-v7.3');
% save topo:
save(fullfile(foldername,'SampleSlices','topo_interp.mat'),'topo_interp','-v7.3');
% save masks:
save(fullfile(foldername,'SampleSlices','mask_interp.mat'),'mask_interp','-v7.3');
save(fullfile(foldername,'SampleSlices','mask.mat'),'mask','-v7.3');
% save figure
if rotate_area==1
    saveas(fig1,fullfile(foldername,'SampleSlices','area.png'));
end

% write config file:
fid=fopen(fullfile(foldername,'SampleSlices','configuration.txt'),'wt');
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
rot=[-45:5:45];
new=zeros(size(xy));
for r=1:length(rot)
    rmat=[cosd(rot(r)) -sind(rot(r)); sind(rot(r)) cosd(rot(r))]; % rotational matrix
    for rr=1:length(xy(:,1))
        new(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
    end
    area(r)=(max(new(:,1))-min(new(:,1)))*(max(new(:,2))-min(new(:,2)));    % area size of rotated coordinates
end
rotbest=rot(area==min(area)); % best rotation angle => smallest area
disp(['Optimum rotation angle is ',num2str(rotbest),' degree. Area has been rotated. Saving coordtrans.mat for later transformation.'])
rmat=[cosd(rotbest) -sind(rotbest); sind(rotbest) cosd(rotbest)]; % rotational matrix
for rr=1:length(xy(:,1))
    new(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
end
% move origin
shiftx=floor(min(new(:,1)));
shifty=floor(min(new(:,2)));
new(:,1)=new(:,1)-shiftx;
new(:,2)=new(:,2)-shifty;
% save coordinate pairs for later transformation
coordtrans=[new(new(:,1)==min(new(:,1)),:) xy(new(:,1)==min(new(:,1)),:);...
    new(new(:,1)==max(new(:,1)),:) xy(new(:,1)==max(new(:,1)),:);...
    new(new(:,2)==min(new(:,2)),:) xy(new(:,2)==min(new(:,2)),:);...
    new(new(:,2)==max(new(:,2)),:) xy(new(:,2)==max(new(:,2)),:)]; % [local x, local y, global x, global y]
% overwrite coordinates in position
xy=new;
disp(['Area size is now ',int2str(round(range(new(:,1)))),' x ',int2str(round(range(new(:,2)))),' m (x/y).'])
end

function [xy]=apply_rotatearea(xy,rot,shiftx,shifty)
%%% Rotate area with given parameters
rmat=[cosd(rot) -sind(rot); sind(rot) cosd(rot)]; % rotational matrix
for rr=1:length(xy(:,1))
    xy(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
end
% move origin
xy(:,1)=xy(:,1)-shiftx;
xy(:,2)=xy(:,2)-shifty;
end
