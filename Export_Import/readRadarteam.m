function [traces,coords,header]=readRadarteam(filename,filter_coords,m,nstd,show_filter_plot,coords_smooth,nsmooth,utmzone)


%%% [traces,coords,header]=readRadarteam(filename,filter_coords,m,nstd,show_filter_plot,coords_smooth,nsmooth,utmzone)
%
% filename: Name of sgy file
% filter_coords: if =1: coordinates are filtered and traces at
% approx. the same position will be removed (e.g. at the beginning or
% end of the profile): a local difference
% over m coordinates for x and y is calculated separately. If this
% absolute local difference is smaller than nstd*std(diff) it will be removed.
% m: number of samples for local difference calculation
% nstd: nstd*standard deviation as valid range
% show_filter_plot: if =1: show plot for quality control of filtering
% coords_smooth: if=1: smooth coordinates along profile
% nsmooth: number of traces for smoothing of coords
% utmzone: UTM zone for transforming into UTM coordinates
%
% Output:
% traces: matrix with traces in columns
% coords: matrix with trace number, Easting, Northing, Height
% header: Header for file
%
% Dr. Tina Wunderlich, CAU Kiel, September 2023
%%%


%%% Read GPR data ---------------------------------------------------------
fid=fopen(filename,'r');

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
ntr=(filesize-3200-400)/(240+2*header.ns); % number of traces
t=0:header.dt*1e-3:(header.ns-1)*header.dt*1e-3; % time range in ns


%% trace header (each 240 byte) and data (2*ns byte)
traces=zeros(header.ns,ntr);
for i=1:ntr

    fseek(fid,3600+(i-1)*(240+2*header.ns),'bof'); % go to beginning of next trace header

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
    fseek(fid,3600+i*240+(i-1)*(2*header.ns),'bof'); % go to beginning of next trace data

    traces(:,i) = fread(fid,512,'int16');

    %%% extract coordinates:
    % coordinates in degree as string:
    wgs_lat_string{i}=num2str(trh(i).LatitudeDouble,20);
    wgs_lon_string{i}=num2str(trh(i).LongitudeDouble,20);
end

fclose(fid);

% Convert Coordinates from string to double
temp=cellfun(@(x) strsplit(x,'.'),wgs_lat_string,'UniformOutput',false); % split at .
lat_deg=cellfun(@(x) str2num(x{1}(1:end-2))+str2num([x{1}(end-1:end),'.',x{2}])/60,temp);

temp=cellfun(@(x) strsplit(x,'.'),wgs_lon_string,'UniformOutput',false); % split at .
lon_deg=cellfun(@(x) str2num(x{1}(1:end-2))+str2num([x{1}(end-1:end),'.',x{2}])/60,temp);

% Coordinate transformation to utm
[Easting,Northing]=wgs2utm(lat_deg,lon_deg,utmzone,'N');
utm=[[1:ntr]' Easting' Northing' h];

% test if coordinates are very strange...
if max(abs(diff(utm(:,2))))<1e-6
    disp(['There is a coordinate problem with ',filename])
    disp(['  -> Skip this file.'])
    traces=[];
    coords=[];
    header=[];
    return
end

% filter coordinates
if filter_coords==1
    % local differences over m sample
    for j=(m-1)/2+1:length(utm(:,1))-(m-1)/2
        Ndiff(j)=diff(utm([j-(m-1)/2 j+(m-1)/2],3));
        Ediff(j)=diff(utm([j-(m-1)/2 j+(m-1)/2],2));
    end
    Ndiff(j+1:length(utm(:,1)))=0;
    Ediff(j+1:length(utm(:,1)))=0;

    % mean differences:
    meandiffE=mean(Ediff);
    meandiffN=mean(Ndiff);

    % standard deviation of local diffs
    Nstd=std(Ndiff);
    Estd=std(Ediff);

    % find good points
    inE=(Ediff>=meandiffE-Estd*nstd & Ediff<=meandiffE+Estd*nstd);
    inN=(Ndiff>=meandiffN-Nstd*nstd & Ndiff<=meandiffN+Nstd*nstd);

    utmfilt=utm(inE & inN,:); % filtered utm list

    if coords_smooth==1 && ~isempty(utmfilt)
        ind=[utmfilt(1,1):utmfilt(end,1)]';
        x=interp1(utmfilt(:,1),utmfilt(:,2),ind); % for all traces
        y=interp1(utmfilt(:,1),utmfilt(:,3),ind);
        z=interp1(utmfilt(:,1),utmfilt(:,4),ind);

        xs=smooth(x,nsmooth);
        ys=smooth(y,nsmooth);
        zs=smooth(z,nsmooth);

        utmfiltsmooth=zeros(size(utmfilt));
        for k=1:length(utmfilt(:,1))
            utmfiltsmooth(k,:)=[utmfilt(k,1) xs(utmfilt(k,1)==ind) ys(utmfilt(k,1)==ind) zs(utmfilt(k,1)==ind)];
        end
    end


    if show_filter_plot==1
        figure
        subplot(2,2,1)
        plot(Ndiff,'*')
        hold on
        yline(meandiffN+Nstd*nstd,'r','LineWidth',2)
        yline(meandiffN-Nstd*nstd,'r','LineWidth',2)
        xlabel('Sample number')
        ylabel('Northing difference [m]')
        set(gca,'Fontsize',20)
        grid on

        subplot(2,2,2)
        plot(Ediff,'*')
        hold on
        yline(meandiffE+Estd*nstd,'r','Linewidth',2)
        yline(meandiffE-Estd*nstd,'r','Linewidth',2)
        xlabel('Sample number')
        ylabel('Easting difference [m]')
        set(gca,'Fontsize',20)
        grid on

        subplot(2,2,3)
        plot(utm(:,2),utm(:,3),'b*','DisplayName','All points')
        hold on
        plot(utm(~inE,2),utm(~inE,3),'kd','DisplayName','filtered out via easting')
        plot(utm(~inN,2),utm(~inN,3),'r+','DisplayName','filtered out via northing')
        set(gca,'FontSize',20)
        xlabel('Easting [m]')
        ylabel('Northing [m]')
        grid on
        legend
        axis equal

        subplot(2,2,4)
        plot(utmfilt(:,2),utmfilt(:,3),'b*','DisplayName','All points filtered')
        hold on
        if exist('utmfiltsmooth','var')
            plot(utmfiltsmooth(:,2),utmfiltsmooth(:,3),'r+','DisplayName','smoothed')
        end
        set(gca,'FontSize',20)
        xlabel('Easting [m]')
        ylabel('Northing [m]')
        grid on
        legend
        axis equal
    end

    clear Ndiff;
    clear Ediff;
end

if exist('utmfiltsmooth','var')
    utm=utmfiltsmooth;
else
    utm=utmfilt;
end

% delete traces without coordinates in beginning and end:
if utm(end,1)<ntr
    traces(:,utm(end,1)+1:end)=[];
end
if utm(1,1)~=1
    traces(:,1:utm(1,1)-1)=[];
end

% interpolate coordinates in between (if some are missing):
coords=zeros(length(traces(1,:)),4);
coords(:,2)=interp1(utm(:,1),utm(:,2),[utm(1,1):utm(end,1)]);
coords(:,3)=interp1(utm(:,1),utm(:,3),[utm(1,1):utm(end,1)]);
coords(:,4)=interp1(utm(:,1),utm(:,4),[utm(1,1):utm(end,1)]);

% and renumber the traces:
coords(:,1)=[1:length(traces(1,:))]';

