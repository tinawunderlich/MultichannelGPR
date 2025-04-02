function [traces,header,coords]=readProceq(filename)


%%% [traces,coords,header]=readProceq(filename)
%
% filename: Name of sgy file
%
% Output:
% traces: matrix with traces in columns
% header: Header for file
% coords: coordinates from csv-file
%
% Dr. Tina Wunderlich, CAU Kiel, April 2025
%%%


%% Read cvs file
[p,n,ext]=fileparts(filename);
temp=strsplit(n,'_');
csvname=[temp{1}];
for i=2:length(temp)-1
    csvname=[csvname ['_',temp{i}]];
end
csvname=[csvname,'.csv'];
lines=readlines(fullfile(p,csvname)); % all lines of file
% extract file name:
temp=char(lines(cellfun(@(x) startsWith(x,'File Name'),lines)));
test=~isstrprop(temp,'wspace');
fn=[];
i=length(test);
while test(i)==1
    fn=[temp(i) fn];
    i=i-1;
end
% extract Repetition Rate [scans/cm]:
temp=char(lines(cellfun(@(x) startsWith(x,'Repetition Rate [scans/cm]'),lines)));
test=~isstrprop(temp,'wspace');
trspacing=[];
i=length(test);
while test(i)==1
    trspacing=[temp(i) trspacing]; % trace spacing in cm
    i=i-1;
end
trspacing=str2num(trspacing);
% extract coordinates:
temp=char(lines(find(cellfun(@(x) startsWith(x,'SCAN DISTANCE [m]'),lines))+1));
test=~isstrprop(temp,'wspace');
starty=[];
startx=[];
len=[];
i=length(test);
while test(i)==1
    starty=[temp(i) starty];
    i=i-1;
end
i=i-1;
while test(i)==1
    startx=[temp(i) startx];
    i=i-1;
end
i=i-1;
while test(i)==1
    len=[temp(i) len];
    i=i-1;
end
startx=str2num(startx);
starty=str2num(starty);
len=str2num(len);




%% Read GPR data ---------------------------------------------------------
fid=fopen(filename,'r','b');

% get file size:
fseek(fid, 0, 'eof');
filesize = ftell(fid); % in byte

fseek(fid,0,'bof'); % go back to beginning

header.TextualFileHeader=fread(fid,3200,'uchar');         % 0-3200


%% segy header (400 byte)
header.Job=fread(fid,1,'int32');                           % after this line is: 3204
header.Line=fread(fid,1,'int32');                          % 3208
header.Reel=fread(fid,1,'int32');                          % 3212
header.DataTraces=fread(fid,1,'int16');        % 3214
header.AuxiliaryTraces=fread(fid,1,'uint16');   % 3216
header.dt=fread(fid,1,'uint16');    % in picoseconds              % 3218
header.dtOrig=fread(fid,1,'uint16');                      % 3220
header.ns=fread(fid,1,'uint16');                          % 3222
header.nsOrig=fread(fid,1,'uint16');                      % 3224
header.DataSampleFormat=fread(fid,1,'int16');            % 3226
% if =3: 16 bit fixed point (integer)
header.NumberOfTraces=fread(fid,1,'int16');       % 3228
fseek(fid,3252,'bof');
header.MeasurementSystem=fread(fid,1,'int16');               % 3254
% 1= meters, 2=feet


%% calculate additional things:
ntr=(filesize-3200-400)/(240+4*header.ns); % number of traces
t=0:header.dt*1e-3:(header.ns-1)*header.dt*1e-3; % time range in ns


%% trace header (each 240 byte) and data (4*ns byte)
traces=zeros(header.ns,ntr);
for i=1:ntr

    fseek(fid,3600+(i-1)*(240+4*header.ns),'bof'); % go to beginning of next trace header

    %%% trace header:
    temp = fread(fid,7,'int32'); % 4 bytes each
    trh(i).TraceSequenceLine = temp(1);    % before this line is: 0

    temp = fread(fid,4,'int16');
    trh(i).TraceIdCode = temp(1); % 28
    trh(i).NSummedTraces = temp(2); % 30
    trh(i).NStackedTraces = temp(3); % 32
    trh(i).DataUse  = temp(4); % 34

    temp = fread(fid,3,'float');
    trh(i).bla = temp(1);  %36
    trh(i).Altitude  = temp(2);  %40  % mean sea level
    trh(i).Height  = temp(3);  %44  % height of geoid above WGS84 ellipsoid
    h(i,1)=trh(i).Altitude;

    temp = fread(fid,1,'int32');
    trh(i).Direction = temp(1);  %48

    temp = fread(fid,1,'float');
    trh(i).DatumOffset = temp(1);  %52

    temp = fread(fid,8,'int16'); % 2 bytes each
    trh(i).SourceGroupScalar = temp(8);  %70

    temp = fread(fid,4,'float');
    trh(i).Longitude = temp(1); %72
    trh(i).Latitude = temp(2); %76
    trh(i).GroupX  = temp(3); %80
    trh(i).GroupY  = temp(4); %84

    temp = fread(fid,1,'int16');
    trh(i).CoordinateUnits = temp(1);  %88
    % 1: meters/feet, 2: arc seconds (DDMM.SSSS)

    temp = fread(fid,4,'int32');
    trh(i).GPSQuality  = temp(1);  %90

    temp = fread(fid,4,'int16');
    trh(i).DelayRecordingTime	 = temp(2);  %108   in picoseconds

    temp = fread(fid,2,'uint16');
    trh(i).ns = temp(1);  %114
    trh(i).dt = temp(2);  %116 in picoseconds

    temp = fread(fid,32,'int16');
    trh(i).HourOfDay					= temp(22);  %160
    trh(i).MinuteOfHour				= temp(23);  %162
    trh(i).SecondOfMinute				= temp(24);  %164
    trh(i).TimeBaseCode				= temp(25);  %166
    % 1:local, 2:GMT

    temp = fread(fid,2,'double');
    trh(i).LongitudeDouble=temp(1);     % 182
    trh(i).LatitudeDouble=temp(2);         % 190

    %%% data trace:
    fseek(fid,3600+i*240+(i-1)*(4*header.ns),'bof'); % go to beginning of next trace data

    traces(:,i) = fread(fid,header.ns,'single');

end

fclose(fid);

%% coordinates:
coords=zeros(ntr,2);
coords(:,1)=linspace(startx,len,ntr);