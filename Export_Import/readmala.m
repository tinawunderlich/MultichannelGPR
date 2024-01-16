function [traces,dt,ns,x,y,z,numchannels]=readmala(foldername,name,profile_num,changeDir,add_Yoffset)

%
% Read Mala rd3-data (all channels in one file)
% [traces,dt,ns,x,y,z,numchannels]=readmala(foldername,name,profile_num,changeDir,add_Yoffset)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires helmert.m
%
% Input:
% foldername: complete rSlicer-foldername and path
% name: Name of datafiles without '_...'
% profile_num: Numbers of profiles
% changeDir: flag if you want to change the +/- of the
% y-antenna-GPS-offset, yes=1, no=0
% add_Yoffset: additional offset for y-antenna-GPS-offset
%
% Output:
% traces: matrix with traces of one profile (channels in cells)
% dt: sampling interval in ns
% ns: Number of samples
% x,y,z: matrices with coordinates in m (channels in cells)
% numchannels: number of channels

if nargin==3
    changeDir=0;
    add_Yoffset=0;
elseif nargin==4
    add_Yoffset=0;
end


%%% Read header file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.rad']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.rad']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.rad']),'r');
end

if fid==-1
    traces=[];
    dt=[];
    ns=[];
    x=[];
    y=[];
    z=[];
    numchannels=[];
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
        num_traces2=length(pos_new(:,1))*numchannels; % update number of traces, if less valid coordinates are present

        % create positioning file
        pos=zeros(num_traces2,8);   % profile number, trace number, channel number, x, y, z, ns, dt
        pos(:,1)=profile_num;   % profile number
        
        tr=pos_new(1,1);
        for i=1:numchannels:num_traces2
            pos(i:i+numchannels-1,2)=tr;   % fill in trace number
            pos(i:i+numchannels-1,6)=pos_new(tr-pos_new(1,1)+1,4);     % fill in z
            tr=tr+1;
        end
        
        % apply channel offsets
        for jj=1:numchannels   % for each channel
            anz=15; % number of points for direction determination
            for ii=1:num_traces2/numchannels-anz
                dist=sqrt((pos_new(ii,2)-pos_new(ii+anz,2))^2+(pos_new(ii,3)-pos_new(ii+anz,3))^2);
                pos(jj+(ii-1)*numchannels,4:5)=helmert([ch_x(jj) ch_y(jj)],[0 0; 0 dist],[pos_new(ii,2) pos_new(ii,3); pos_new(ii+anz,2) pos_new(ii+anz,3)]);
                pos(jj+(ii-1)*numchannels,3)=jj;  % channel number
            end
            anz1=anz;
            % calculation for the last few traces
            anz=anz-1;
            for ii=num_traces2/numchannels-anz1+1:num_traces2/numchannels-1
                dist=sqrt((pos_new(ii,2)-pos_new(ii+anz,2))^2+(pos_new(ii,3)-pos_new(ii+anz,3))^2);
                pos(jj+(ii-1)*numchannels,4:5)=helmert([ch_x(jj) ch_y(jj)],[0 0; 0 dist],[pos_new(ii,2) pos_new(ii,3); pos_new(ii+anz,2) pos_new(ii+anz,3)]);
                pos(jj+(ii-1)*numchannels,3)=jj;  % channel number
                anz=anz-1;
            end
            % extrapolate for last trace
            pos(jj+(ii-1)*numchannels+numchannels,4)=interp1([1 2],[pos(jj+(ii-1)*numchannels-numchannels,4) pos(jj+(ii-1)*numchannels,4)],3,'linear','extrap');
            pos(jj+(ii-1)*numchannels+numchannels,5)=interp1([1 2],[pos(jj+(ii-1)*numchannels-numchannels,5) pos(jj+(ii-1)*numchannels,5)],3,'linear','extrap');
            pos(jj+(ii-1)*numchannels+numchannels,3)=jj;  % channel number
        end
        
        num_rows=length(pos(:,1));
        
       % pos: profile number, trace number, channel number, x, y, z   
    end
end
for i=1:numchannels
    x{i}=pos(pos(:,3)==i,4);
    y{i}=pos(pos(:,3)==i,5);
    z{i}=pos(pos(:,3)==i,6);
end
    

%%% Read radar data file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.rd3']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.rd3']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.rd3']),'r');
end

traces=cell(numchannels,1); % initialize cell array (one cell for each channel)

for i=1:numchannels
    traces{i}=zeros(ns,num_traces/numchannels);
end
anz=1;
for ii=1:num_traces    % for all traces in this file (also without coords)
    temp1=fread(fid,ns,'int16','l');    % read trace
    currChan=mod(ii-1,numchannels)+1;  % current channel number
    traces{currChan}(:,anz)=temp1;
    if currChan==numchannels
        anz=anz+1;
    end
end
fclose(fid);

% delete traces without coords:
in=pos_new(:,1); % valid trace numbers with coordinates
for i=1:numchannels
    traces{i}=traces{i}(:,in);
end

dt=range/ns;