function [traces,dt,ns,x,y,z] = readGSSI(foldername,UTMzone,profile_num)

% Read GSSI DZT-data (up to two channels)
% this is just an interface to provide functionality like the readmala.m
% function by Dr. Tina Wunderlich
% Input:
% foldername: Data foldername and path
% name: Name of datafiles without '_...'
% profile_num: Numbers of profile to open
%
% Output:
% traces: matrix with traces of one profile (channels in cells)
% dt: sampling interval in ns
% ns: Number of samples
% x,y,z: matrices with coordinates in m (channels in cells)

% get filename of profile
list = dir(fullfile(foldername,'/*.DZT'));
dztfile = fullfile(foldername,list(profile_num).name);

% read data
[data,trh,h]=readdzt_4000(dztfile,0);
dt = h.dt;
ns = h.ns;

% UTM coversion
lontemp=num2str(trh.x','%.8f');
lattemp=num2str(trh.y','%.8f');
for j=1:length(lattemp(:,1))
    temp=strsplit(lattemp(j,:),'.');
    lat(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
    temp=strsplit(lontemp(j,:),'.');
    lon(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
end
[xneu,yneu]=wgs2utm(lat,lon,UTMzone,'N');

% in case of two channels, split the data
 if h.nchan == 2
     traces{1} = data(:,1:h.numtraces);
     traces{2} = data(:,h.numtraces+1:end);
     x{1} = xneu(1:h.numtraces)';
     x{2} = xneu(h.numtraces+1:end)';
     y{1} = yneu(1:h.numtraces)';
     y{2} = yneu(h.numtraces+1:end)';
     z{1} = trh.z(1:h.numtraces)';
     z{2} = trh.z(h.numtraces+1:end)';
else
    traces{1} = data;
    x{1} = xneu';
    y{1} = yneu';
    z{1} = trh.z';
end

end

