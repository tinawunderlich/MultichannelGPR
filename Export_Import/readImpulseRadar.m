function [traces,dt,ns,x,y,z,numchannels,flag]=readImpulseRadar(foldername,name,profile_num,add_Yoffset,GNSS_height,utmzone,filter_coords,smooth_coords)

% Read Impulse radar-data (Raptor)
% [traces,dt,ns,x,y,z,numchannels]=readImpulseRadar(foldername,name,profile_num,add_Yoffset,GNSS_height,utmzone,filter_coords,smooth_coords)
%
% Dr. Tina Wunderlich, CAU Kiel 2025, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% foldername: complete path
% name: Name of datafiles without '_...'
% profile_num: Number of profile
% add_Yoffset: additional offset for y-antenna-GPS-offset
% GNSS_height: height of GNSS antenna above ground (for tilting correction)
%
% Output:
% traces is matrix of size [ns,numchannels*numtraces_per_channel] 
% (order is: first all traces of 1. channels, all traces of 2. channel,....)
% x, y, z: vectors of length [1,numchannels*numtraces_per_channel] with
% coordinates for each trace
% dt: sampling interval in ns
% ns: number of samples per trace
% numchannels: number of channels
%
% requires helmert.m

if nargin==3
    changeDir=0;
    add_Yoffset=0;
    GNSS_height=0;
elseif nargin==4
    add_Yoffset=0;
    GNSS_height=0;
elseif nargin==5
    GNSS_height=0;
end


%%% Read header file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_A01.iprh']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_A01.iprh']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_A01.iprh']),'r');
end

flag=1;
if fid==-1
    traces=[];
    dt=[];
    ns=[];
    x=[];
    y=[];
    z=[];
    numchannels=[];
    flag=0; % -> bad file
    return;
end

i=1;
while ~feof(fid)
    temp{i}=fgetl(fid);
    
    if startsWith(temp{i},'SAMPLES:')
        ns=str2double(temp{i}(length('SAMPLES:')+1:end));
    elseif strfind(temp{i},'TIMEWINDOW:')
        range=str2double(temp{i}(length('TIMEWINDOW:')+1:end));
    elseif strfind(temp{i},'LAST TRACE:')
        num_traces=str2double(temp{i}(length('LAST TRACE:')+1:end));
    elseif strfind(temp{i},'CHANNELS:')
        numchannels=str2double(temp{i}(length('CHANNELS:')+1:end));
    elseif strfind(temp{i},'CH_X_OFFSET:')
        ch_x(1)=str2double(temp{i}(length('CH_X_OFFSET:')+1:end));
    elseif strfind(temp{i},'CH_Y_OFFSET:')
        ch_y(1)=str2double(temp{i}(length('CH_Y_OFFSET:')+1:end));
    end
    i=i+1;
end
fclose(fid);

% read header files of other channels:
for i=2:numchannels
    if profile_num<10
        if i<10
            fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_A0',int2str(i),'.iprh']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_A',int2str(i),'.iprh']),'r');
        end
    elseif profile_num>=10 && profile_num<100
        if i<10
            fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_A0',int2str(i),'.iprh']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_A',int2str(i),'.iprh']),'r');
        end
    else
        if i<10
            fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_A0',int2str(i),'.iprh']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_A',int2str(i),'.iprh']),'r');
        end
    end

    while ~feof(fid)
        temp=fgetl(fid);
        if strfind(temp,'CH_X_OFFSET:')
            ch_x(i)=str2double(temp(length('CH_X_OFFSET:')+1:end));
        elseif strfind(temp,'CH_Y_OFFSET:')
            ch_y(i)=str2double(temp(length('CH_Y_OFFSET:')+1:end));
        end
    end
    fclose(fid);
end
ch_y=ch_y+add_Yoffset; % add additional offset


%%% Read coordinate file
if utmzone~=0 % GNSS
    if profile_num<10
        fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.cor']),'r');
    elseif profile_num>=10 && profile_num<100
        fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.cor']),'r');
    else
        fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.cor']),'r');
    end
    if fid~=-1
        temp=textscan(fid,'%f%s%s%f%*c%f%*c%f%*c%f','Headerlines',1);
        fclose(fid);
        pos_orig=[temp{1} temp{4} temp{5} temp{6}];  % trace number, x, y, z
    end
else % Total station
    if profile_num<10
        fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.tsp']),'r');
    elseif profile_num>=10 && profile_num<100
        fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.tsp']),'r');
    else
        fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.tsp']),'r');
    end
    if fid~=-1
        temp=textscan(fid,'%f%f%f%f%s%f','Headerlines',5);
        fclose(fid);
        pos_orig=[temp{1} temp{2} temp{3} temp{4}];  % trace number, x, y, z
        [a,b]=unique(pos_orig(:,1));
        pos_orig=pos_orig(b,:); % only unique points
    end
end


if utmzone~=0
    % convert lat/lon to utm
    [pos_orig(:,2),pos_orig(:,3)]=wgs2utm(pos_orig(:,2),pos_orig(:,3),utmzone,'N');
end

if smooth_coords~=0
    % smooth coordinates with moving median filter
    hw=(smooth_coords-1)/2; % half width of window
    % pad coords
    temp=[repmat(pos_orig(1,2:4),[hw,1]); pos_orig(:,2:4); repmat(pos_orig(end,2:4),[hw,1])];
    % median filter
    for i=hw+1:size(pos_orig,1)-hw
        pos_orig(i-hw,2:4)=median(temp(i-hw:i+hw,:),1);
    end
end

if ~isempty(pos_orig) && num_traces/numchannels>=15  % if number of traces in file is too small -> omit file (=do not save in position-matrix)
    % interpolate between traces
    pos_new=[[pos_orig(1,1):pos_orig(end,1)]' interp1(pos_orig(:,1),pos_orig(:,2),[pos_orig(1,1):pos_orig(end,1)]') interp1(pos_orig(:,1),pos_orig(:,3),[pos_orig(1,1):pos_orig(end,1)]') interp1(pos_orig(:,1),pos_orig(:,4),[pos_orig(1,1):pos_orig(end,1)]')];  % trace number, x(Rw!), y(Hw!), z

    if filter_coords==1
        % filter coordinates (delete points at the same position)
        % Calculate rate of change between coordinates:
        rate=sqrt(diff(pos_new(:,2)).^2+diff(pos_new(:,3)).^2)./diff(pos_new(:,1));
        rate=[mean(rate); rate];
        del = (rate<mean(rate)/2); % Rate too small -> delete points
        trnum=pos_new(:,1);
        del=trnum(del); % del is now containing the trace numbers that need to be deleted later in the data
        if sum(del)<length(del) % if less than half of the points will be deleted
            pos_new(del,:)=[];
            % set new trace number:
            pos_new(:,1)=1:size(pos_new,1);
        end
    end

    if smooth_coords~=0
        % smooth coordinates with moving median filter
        hw=smooth_coords; % half width of window (increase smooth_coords*2)
        % pad coords
        temp=[repmat(pos_new(1,2:4),[hw,1]); pos_new(:,2:4); repmat(pos_new(end,2:4),[hw,1])];
        % median filter
        for i=hw+1:size(pos_new,1)-hw
            pos_new(i-hw,2:4)=median(temp(i-hw:i+hw,:),1);
        end
    end

    num_traces2=length(pos_new(:,1))*numchannels; % update number of traces (for all channels), if less valid coordinates are present

    % create positioning file
    pos=zeros(num_traces2,8);   % profile number, trace number, channel number, x, y, z, ns, dt
    pos(:,1)=profile_num;   % profile number

    tr=pos_new(1,1);
    for i=1:numchannels:num_traces2
        pos(i:i+numchannels-1,2)=tr;   % fill in trace number
        tr=tr+1;
    end

    % apply channel offsets and correct for tilted GNSS antenna due to topography
    for jj=1:numchannels   % for each channel
        trh(jj).x=pos_new(:,2);
        trh(jj).y=pos_new(:,3);
        trh(jj).z=pos_new(:,4);
        trh(jj).tracenum=1:size(pos_new,1);
        trh(jj).channum=zeros(size(pos_new(:,1)))+jj;
        [~,trh(jj),deltra(:,jj)]=correctCoordinates(zeros(2,length(trh(jj).x)),trh(jj),GNSS_height,ch_x(jj),ch_y(jj),1,1,3);
    end
    % determine traces that are valid in all channels:
    deltraces=~all(deltra==0,2); % valid=0, delete trace=1

    % set channel number
    for i=1:numchannels:sum(deltraces==0)*numchannels
        pos(i:i+numchannels-1,3)=1:numchannels; % fill in channel number
    end

    ind_all=find(deltraces==0);
    for jj=1:numchannels   % for each channel
        % only take those traces given by general deltraces
        ind_chan=find(deltra(:,jj)==0);
        ind_chan(:,2)=0;
        for ii=1:length(ind_chan)
            if any(ind_chan(ii,1)==ind_all) % if channel-index is also in general index list -> ok
                ind_chan(ii,2)=1;
            else
                ind_chan(ii,2)=0; % not in both lists -> delete this trace later
            end
        end
        % set corrected coordinates in pos-variable:
        pos(pos(:,3)==jj,4)=trh(jj).x(ind_chan(:,2)==1);
        pos(pos(:,3)==jj,5)=trh(jj).y(ind_chan(:,2)==1);
        pos(pos(:,3)==jj,6)=trh(jj).z(ind_chan(:,2)==1);
    end
    % delete empty traces in pos:
    pos=pos(pos(:,4)~=0,:);

    num_rows=length(pos(:,1));

    % pos: profile number, trace number, channel number, x, y, z
else
    traces=[];
    dt=[];
    ns=[];
    x=[];
    y=[];
    z=[];
    numchannels=[];
    return;
end
    

%%% Read radar data file
traces1=cell(numchannels,1); % initialize cell array (one cell for each channel)
for i=1:numchannels
    if profile_num<10
        if i<10
            fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_A0',int2str(i),'.iprb']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_A',int2str(i),'.iprb']),'r');
        end
    elseif profile_num>=10 && profile_num<100
        if i<10
            fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_A0',int2str(i),'.iprb']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_A',int2str(i),'.iprb']),'r');
        end
    else
        if i<10
            fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_A0',int2str(i),'.iprb']),'r');
        else
            fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_A',int2str(i),'.iprb']),'r');
        end
    end

    traces1{i}=fread(fid,[ns num_traces],'int16','l');    % read traces

    fclose(fid);

    % delete traces if invalid coordinates:
    traces1{i}=traces1{i}(:,~deltraces);
end

% adjust num_traces2:
num_traces2=size(traces1{1},2)*numchannels; % number of traces for all channels

% sort in one matrix
numtr2=num_traces2/numchannels; % number of traces per channel without filtered traces (end and beginning of profile when measured with PPS)
traces=zeros(ns,num_traces2);
x=zeros(1,num_traces2);
y=zeros(1,num_traces2);
z=zeros(1,num_traces2);
for i=1:numchannels
    traces(:,(i-1)*numtr2+i-(i-1):i*numtr2)=traces1{i};

    tempx=pos(pos(:,3)==i,4);
    tempy=pos(pos(:,3)==i,5);
    tempz=pos(pos(:,3)==i,6);
    nums=1:length(tempx);
    tempx(isnan(tempx))=interp1(nums(~isnan(tempx)),tempx(~isnan(tempx)),nums(isnan(tempx)),'linear','extrap'); % interpolate nan coords
    tempy(isnan(tempy))=interp1(nums(~isnan(tempy)),tempy(~isnan(tempy)),nums(isnan(tempy)),'linear','extrap'); % interpolate nan coords
    tempz(isnan(tempz))=interp1(nums(~isnan(tempz)),tempz(~isnan(tempz)),nums(isnan(tempz)),'linear','extrap'); % interpolate nan coords
    x(1,(i-1)*numtr2+i-(i-1):i*numtr2)=tempx;
    y(1,(i-1)*numtr2+i-(i-1):i*numtr2)=tempy;
    z(1,(i-1)*numtr2+i-(i-1):i*numtr2)=tempz;
end

dt=range/ns;