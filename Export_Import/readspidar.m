function [traces,dt,ns,x,y,z,numchannels]=readspidar(foldername,name,profile_num,antennaIDs,offsetGPS_X,offsetGPS_Y,gps_channel,zone,dataplot)

% Read Sensors&Software Spidar-data (HD,DT,GPS)
% [traces,dt,ns,x,y,z,numchannels]=readspidar(foldername,name,profile_num,antennaIDs,offsetGPS_X,offsetGPS_Y,gps_channel,zone,dataplot)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% foldername: complete foldername and path
% name: cell array with {1} name of datafile and {2} e.g. 'Line'
% profile_num: Number of profile
% antennaIDs: IDs of antennas/channels (same as in the file names) (vector of length(numchannels)
% offsetGPS_x/offset_GPS_y: GPS offsets in x and y direction for each channel (vectors)
% gps_channel: which channel is synchronized with the GPS? (-> same filename-number)
% zone: utm zone, e.g. 32
% dataplot: if =1: plotting for control, if =0 no plot
%
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

%%
% read gps-file
if gps_channel<10
    gpsname_part1=[name{1},'0',int2str(gps_channel)];
else
    gpsname_part1=[name{1},int2str(gps_channel)];
end
if profile_num<10
    gpsname_part2=['_',name{2},'00',int2str(profile_num),'.GPS'];
elseif profile_num<100
    gpsname_part2=['_',name{2},'0',int2str(profile_num),'.GPS'];
else
    gpsname_part2=['_',name{2},int2str(profile_num),'.GPS'];
end
temp=readlines(fullfile(foldername,[gpsname_part1,gpsname_part2]));
g=cellfun(@(x) startsWith(x,'$GPGGA'),temp); % find GGA rows
tr=cellfun(@(x) startsWith(x,'Trace'),temp); % find trace rows
temp=temp(g | tr);
n=1;
for i=1:2:length(temp)
    trnum(n)=str2num(extractBetween(temp(i),'#',' ')); % trace number
    xtr(n)=str2num(extractAfter(temp(i),'position ')); % x position of this trace
    test=strsplit(temp(i+1),','); % GGA string
    % extract the coordinate information
    if strcmp(test{6},'E')
        xx(n) = str2double(test{5});
    else
        xx(n) = -str2double(test{5}); % West of Greenwich -> negative
    end
    yy(n) = str2double(test{3});
    zz(n) = str2double(test{10});
    n=n+1;
end
% convert to utm
for i=1:length(xx)
    lontemp=num2str(xx(i),'%.8f');
    lattemp=num2str(yy(i),'%.8f');
    for j=1:length(lattemp(:,1))
        temp=strsplit(lattemp(j,:),'.');
        lat(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
        temp=strsplit(lontemp(j,:),'.');
        if strcmp(temp{1}(1),'-') % negative longitude
            lon(j)=-1*(str2num(temp{1}(2:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60);
        else
            lon(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
        end
    end
    [xneu(i),yneu(i)]=wgs2utm(lat,lon,zone,'N'); % utm coordinates
end


%%
for i=1:length(antennaIDs) % for each channel:
    % Read header file
    if antennaIDs(i)<10
        HDname_part1=[name{1},'0',int2str(antennaIDs(i))];
    else
        HDname_part1=[name{1},int2str(antennaIDs(i))];
    end
    if profile_num<10
        HDname_part2=['_',name{2},'00',int2str(profile_num),'.HD'];
    elseif profile_num<100
        HDname_part2=['_',name{2},'0',int2str(profile_num),'.HD'];
    else
        HDname_part2=['_',name{2},int2str(profile_num),'.HD'];
    end
    temp=readlines(fullfile(foldername,[HDname_part1,HDname_part2]));
    for line = temp'
        if contains(line, 'NUMBER OF TRACES')
            h.ntr(i) = str2num(extractAfter(line, '= '));
        elseif contains(line, 'TOTAL TIME WINDOW')
            h.tmax(i) = str2num(extractAfter(line, '= '));
        elseif contains(line, 'STARTING POSITION')
            h.start_pos(i) = str2num(extractAfter(line, '= '));
        elseif contains(line, 'FINAL POSITION')
            h.ending_pos(i) = str2num(extractAfter(line, '= '));
        elseif contains(line, 'STEP SIZE USED')
            h.dx(i) = str2num(extractAfter(line, '= '));
        end
    end
end

% interpolate coordinates for every trace number (and make constant trace
% number for all channels)
xn=interp1(trnum,xneu,1:min(h.ntr));
yn=interp1(trnum,yneu,1:min(h.ntr));
zzn=interp1(trnum,zz,1:min(h.ntr));

traces=[];
x=[];
y=[];
z=[];
for i=1:length(antennaIDs) % for each channel:
    % Read data file:
    if antennaIDs(i)<10
        spidarname_part1=[name{1},'0',int2str(antennaIDs(i))];
    else
        spidarname_part1=[name{1},int2str(antennaIDs(i))];
    end
    if profile_num<10
        spidarname_part2=['_',name{2},'00',int2str(profile_num),'.DT1'];
    elseif profile_num<100
        spidarname_part2=['_',name{2},'0',int2str(profile_num),'.DT1'];
    else
        spidarname_part2=['_',name{2},int2str(profile_num),'.DT1'];
    end
    fid=fopen(fullfile(foldername,[spidarname_part1,spidarname_part2]),'r');
    for j=1:min(h.ntr) % read every trace
        % header
        trh.tracenum(j)=fread(fid,1,'float');
        trh.pos(j)=fread(fid,1,'float');
        trh.ns(j)=fread(fid,1,'float');
        tmp=fread(fid,3,'float');
        trh.topo(j)=tmp(1);
        trh.twindow(j)=fread(fid,1,'float');
        trh.nstacks(j)=fread(fid,1,'float');
        trh.x(j)=fread(fid,1,'double');
        trh.y(j)=fread(fid,1,'double');
        trh.z(j)=fread(fid,1,'double');
        tmp=fread(fid,7,'float');
        tmp=fread(fid,11,'float');
        trh.timezero(j)=tmp(1);
        trh.zeroflag(j)=tmp(2); % 0=data ok, 1=no data
        % data
        if j==1
            data=zeros(trh.ns(j),min(h.ntr));
        end
        data(:,j)=fread(fid,trh.ns(j),'int16');
    end
    fclose(fid);

    trh.x=xn;
    trh.y=yn;
    trh.z=zzn;

    if dataplot==1
        figure
        imagesc(1:length(data(1,:)),linspace(0,h.tmax(1),trh.ns(1)),data)
        colormap(flipud(gray))
        xlabel('Trace number')
        ylabel('t [ns]')
        coldata=sort(data(~isnan(data)));
        cmin=coldata(round(length(coldata)/100*1));
        cmax=coldata(end-round(length(coldata)/100*1));
        set(gca,'CLim',[cmin cmax])
    end


    % correct GPS-antenna offset:
    if offsetGPS_X(i)~=0 || offsetGPS_Y(i)~=0
        anz2=round(0.5/mean(sqrt(diff(trh.x).^2+diff(trh.y).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
        if anz2/2==round(anz2/2)
            anz2=anz2+1; % make odd
        end
        for ii=1:length(trh.x)-anz2
            dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
            temp=helmert([offsetGPS_X(i) offsetGPS_Y(i)],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
            trh.x(ii)=temp(1);
            trh.y(ii)=temp(2);
        end
        anz1=anz2;
        % calculation for the last few traces
        anz2=anz2-1;
        for ii=length(trh.x)-anz1+1:length(trh.x)-1
            dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
            temp=helmert([offsetGPS_X(i) offsetGPS_Y(i)],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
            trh.x(ii)=temp(1);
            trh.y(ii)=temp(2);
            anz2=anz2-1;
        end
        % extrapolate for last trace
        trh.x(end)=interp1([1 2],[trh.x(end-2) trh.x(end-1)],3,'linear','extrap');
        trh.y(end)=interp1([1 2],[trh.y(end-2) trh.y(end-1)],3,'linear','extrap');
    end

    traces=[traces data];
    x=[x trh.x];
    y=[y trh.y];
    z=[z trh.z];
    ns=trh.ns(1); % number of samples
    clear trh;
end

% prepare output:
dt=h.tmax(1)/ns; % sampleinterval
numchannels=length(antennaIDs);