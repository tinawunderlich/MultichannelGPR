function [traces,dt,ns,x,y,z,numchannels,flag]=readmala4parfor(foldername,name,profile_num,changeDir,add_Yoffset,GNSS_height)

% Read Mala rd3-data (all channels in one file)
% [traces,dt,ns,x,y,z,numchannels]=readmala4parfor(foldername,name,profile_num,changeDir,add_Yoffset,GNSS_height)
%
% Dr. Tina Wunderlich, CAU Kiel 2020-2025, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% foldername: complete rSlicer-foldername and path
% name: Name of datafiles without '_...'
% profile_num: Numbers of profiles
% changeDir: flag if you want to change the +/- of the
% y-antenna-GPS-offset, yes=1, no=0
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
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.rad']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.rad']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.rad']),'r');
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
    
    if strfind(temp{i},'SAMPLES:')
        ns=str2double(temp{i}(length('SAMPLES:')+1:end));
    elseif strfind(temp{i},'TIMEWINDOW:')
        range=str2double(temp{i}(length('TIMEWINDOW:')+1:end));
    elseif strfind(temp{i},'LAST TRACE:')
        num_traces=str2double(temp{i}(length('LAST TRACE:')+1:end));
    elseif strfind(temp{i},'NUMBER_OF_CH:')
        numchannels=str2double(temp{i}(length('NUMBER_OF_CH:')+1:end));
    elseif strfind(temp{i},'CH_X_OFFSETS:')
        tmp=temp{i}(length('CH_X_OFFSETS:')+1:end);
        ind=[1 strfind(tmp,' ')];   % indices of separating spaces
        ch_x=zeros(numchannels,1);
        for j=1:length(ind)-1
            ch_x(j)=str2double(tmp(ind(j):ind(j+1)));  % channel offsets in x-direction
        end
    elseif strfind(temp{i},'CH_Y_OFFSETS:')
        tmp=temp{i}(length('CH_Y_OFFSETS:')+1:end);
        ind=[1 strfind(tmp,' ')];   % indices of separating spaces
        ch_y=zeros(numchannels,1);
        for j=1:length(ind)-1
            ch_y(j)=str2double(tmp(ind(j):ind(j+1))); % channel offsets in y-direction
        end
        if changeDir==1
            ch_y=ch_y.*(-1); % change +/- sign for offset
        end
        ch_y=ch_y+add_Yoffset; % add additional offset
    end
    i=i+1;
end
fclose(fid);


%%% Read positioning file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_G01.pos']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_G01.pos']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_G01.pos']),'r');
end
if fid~=-1
    temp=textscan(fid,'%f%f%f%f','Headerlines',1);
    fclose(fid);
    pos_orig=[temp{1} temp{2} temp{3} temp{4}];  % trace number, x, y, z
    
    if ~isempty(pos_orig) && num_traces/numchannels>=15  % if number of traces in file is too small -> omit file (=do not save in position-matrix)
        % interpolate between traces
        pos_new=[[pos_orig(1,1):pos_orig(end,1)]' interp1(pos_orig(:,1),pos_orig(:,3),[pos_orig(1,1):pos_orig(end,1)]') interp1(pos_orig(:,1),pos_orig(:,2),[pos_orig(1,1):pos_orig(end,1)]') interp1(pos_orig(:,1),pos_orig(:,4),[pos_orig(1,1):pos_orig(end,1)]')];  % trace number, x(now Rw!), y(now Hw!), z
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
end
    

%%% Read radar data file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.rd3']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.rd3']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.rd3']),'r');
end

traces1=cell(numchannels,1); % initialize cell array (one cell for each channel)
for i=1:numchannels
    traces1{i}=zeros(ns,num_traces/numchannels);
end
anz=1;
for ii=1:num_traces    % for all traces in this file
    temp1=fread(fid,ns,'int16','l');    % read trace
    currChan=mod(ii-1,numchannels)+1;  % current channel number
    traces1{currChan}(:,anz)=temp1;
    if currChan==numchannels
        anz=anz+1;
    end
end
fclose(fid);

% delete traces if invalid coordinates:
for jj=1:numchannels
    traces1{jj}=traces1{jj}(:,~deltraces);
end

% adjust num_traces2:
num_traces2=size(traces1{1},2)*numchannels;

% sort in one matrix
numtr2=num_traces2/numchannels; % number of traces per channel without filtered traces (end and beginning of profile when measured with PPS)
traces=zeros(ns,num_traces2);
x=zeros(1,num_traces2);
y=zeros(1,num_traces2);
z=zeros(1,num_traces2);
for i=1:numchannels
    traces(:,(i-1)*numtr2+i-(i-1):i*numtr2)=traces1{i};
    x(1,(i-1)*numtr2+i-(i-1):i*numtr2)=pos(pos(:,3)==i,4);
    y(1,(i-1)*numtr2+i-(i-1):i*numtr2)=pos(pos(:,3)==i,5);
    z(1,(i-1)*numtr2+i-(i-1):i*numtr2)=pos(pos(:,3)==i,6);
end

dt=range/ns;