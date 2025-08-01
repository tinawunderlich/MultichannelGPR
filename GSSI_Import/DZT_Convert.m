%------------------------ DZT_Convert -------------------------------------
%
% Converts DZT-files of all GSSI equipments to segy and mat (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, August 2020, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


clear all
close all
clc

app='Tablet';    % Equipment: SIR20 / SIR30 / SIR3000 / SIR4000 / Tablet / UtilityScan (UtilityScan with DF antenna only!)

dataplot=1; % plot radargram for controlling? 1=yes, 0=no

convert2utm=0; % convert WGS Lat/Long to UTM (=1 if measured with Stonex-GPS)
zone=33; % if convert2utm==1 -> give UTM-zone

offsetGPS_X=0; % [m] Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=0; % [m] Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)
h_GPS=0.0; % [m] height of GPS/prism above ground

% smoothing of GPS-coordinates before applying antenna-offsets:
smooth_coords=1; % yes=1, no=0
numsamp_smooth=5; % number of samples for moving median smoothing

% Options for calculating inline coordinates for each trace:
coords_opt=1;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
                % =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)

% Attention: still experimental! If you think that the output radargrams
% are wrong, set to 0!
removeOutliers=0; % do you want to remove coordinate outliers? 
    % 0= no, use raw coordinates
    % 1= in middle of profile
    % 2= at the end/beginning of profile

% remove start/end traces of profiles
removeStartEnd=0; % =1: yes, remove start and end traces of profiles (at same position)
                    % =0: No. Keep all traces.

% Export to other formats
export2mat=1; % export to Multichannel-GPR format for radargrams (mat-files)
export2segy=0; % export all radargrams as segy-files
constoff=0; % if=1: a constant coordinate offset will be subtracted and coordinates will be in mm accuracy in segy file (offsets will be saved in Inline3D (x) and Crossline3D (y))


%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
% 
% LOAD DATA

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));


test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end

% get file names
if ispc
    if exist('dzt.temp') % read last opened folder from temp.temp
        fid=fopen('dzt.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZT-file(s)');
        else
            pfad=uigetdir([],'Choose folder with DZT-file(s)');
        end
        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with DZT-file(s)'); % path to radargram-folder

        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.dzt.temp') % read last opened folder from temp.temp
        fid=fopen('.dzt.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZT-file(s)');
        else
            pfad=uigetdir([],'Choose folder with DZT-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with DZT-file(s)'); % path to radargram-folder
    end

    fid=fopen('.dzt.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*.DZT')); % GSSI format
% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end

if convert2utm==0 % no conversion
    zone=0;
end

radargrams=[];
global_coords=[];
x=[];
marker=[];
anz=1;
listrad.name=[];
listrad.chan=[];
for i=1:length(list)
    disp(['File ',int2str(i),' of ',int2str(length(list))])
    
    dztname=list(i).name;
        [~,name{i},ext]=fileparts(dztname);
        

        if strcmp(app,'SIR4000')
            [data,trh,h]=readdzt_4000(fullfile(pfad,[name{i},'.DZT']),dataplot,removeOutliers,zone);
        elseif strcmp(app,'Tablet')
            [data,trh,h]=readdzt_Tablet(fullfile(pfad,[name{i},'.DZT']),dataplot);
        elseif strcmp(app,'UtilityScan')
            [data,trh,h]=readdzt_UtilityScan(fullfile(pfad,[name{i},'.DZT']),dataplot);
        elseif strcmp(app,'SIR30')
            [data,trh,h]=readdzt_30(fullfile(pfad,[name{i},'.DZT']),dataplot);
        else
            [data,trh,h]=readdzt_all(fullfile(pfad,[name{i},'.DZT']),dataplot,app);
        end
        
        if ~strcmp(app,'SIR4000') && convert2utm==1 && exist([pfad,'/', name{i}, '.DZG'])
            lontemp=num2str(trh.x','%.8f');
            lattemp=num2str(trh.y','%.8f');
            for j=1:length(lattemp(:,1))
                temp=strsplit(lattemp(j,:),'.');
                lat(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
                temp=strsplit(lontemp(j,:),'.');
                lon(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
            end
            [xneu,yneu]=wgs2utm(lat,lon,zone,'N');
            trh.x=xneu;
            trh.y=yneu;
            clear lat lon;
        end
        
        % correct GPS-antenna offset:
        if offsetGPS_X~=0 || offsetGPS_Y~=0 || h_GPS~=0
            % if there are traces at the same position -> delete
            d1=diff(trh.x);
            d2=diff(trh.y);
            d1=[d1 d1(end)];
            d2=[d2 d2(end)];
            weg=(d1==0 & d2==0);
            if removeStartEnd==1
                diffe = sqrt(diff(trh.x).^2 + diff(trh.y).^2);
                standstill = diffe < mean(diffe);
                first = find(~standstill,1,'first')-1;
                last = length(trh.x) - find(fliplr(~standstill),1,'first')-1;
                weg(1:first) = 1;
                weg(last:end) = 1;
            end
            trh.x=trh.x(~weg);
            trh.y=trh.y(~weg);
            trh.z=trh.z(~weg);
            trh.time=trh.time(~weg);
            if isfield(trh,'quality')
                trh.quality=trh.quality(~weg);
            end
            trh.mark=trh.mark(~weg);
            trh.channum=trh.channum(~weg);
            trh.tracenum=1:length(trh.mark);
            data=data(:,~weg);
            % adjust header
            h.numtraces=size(data,2)/h.nchan;

            
            if smooth_coords==1
                for ii=unique(trh.channum)
                    trh.x(trh.channum==ii)=movmean(trh.x(trh.channum==ii),numsamp_smooth);
                    trh.y(trh.channum==ii)=movmean(trh.y(trh.channum==ii),numsamp_smooth);
                    trh.z(trh.channum==ii)=movmean(trh.z(trh.channum==ii),numsamp_smooth);
                end
            end


            % now correct for offsets:
            anz2=round(0.5/mean(sqrt(diff(trh.x).^2+diff(trh.y).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
            if anz2/2==round(anz2/2)
                anz2=anz2+1; % make odd
            end
            for ii=1:length(trh.x)-anz2
                dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2); % distance between two points anz2-traces away
                dist_xy=[trh.x(ii+anz2)-trh.x(ii) trh.y(ii+anz2)-trh.y(ii)]; % differences in x and y direction
                % correct for GPS height and associated xy-shift if
                % topography is present
                if any(trh.z~=0)
                    alpha=atand(abs(trh.z(ii+anz2)-trh.z(ii))/dist); % slope angle in degree
                    hnew=h_GPS/cosd(alpha); % new height of GPS above ground at measured point
                    b=hnew*sind(alpha); % shift along surface to correct point
                    d_x=b*cosd(alpha); % horizontal shift of point in profile direction
                    d_z=b*sind(alpha); % vertical shift of point in profile direction

                    % angle of profile direction towards north
                    if trh.x(ii)==trh.x(ii+anz2) % north/south
                        beta=0;
                    elseif trh.y(ii)==trh.y(ii+anz2) % west/east
                        beta=90;
                    else
                        beta=atand(abs(trh.x(ii+anz2)-trh.x(ii))/abs(trh.y(ii+anz2)-trh.y(ii)));
                    end
                    % split horizontal shift into x and y component:
                    d=[d_x*cosd(beta) d_x*sind(beta)];
                    % differentiate between different cases
                    updownflag=trh.z(ii)<trh.z(ii+anz2); % up=1, down=0
                    if all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards NNE
                        d_xx=min(d);
                        d_xy=max(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards NNE
                        d_xx=-min(d);
                        d_xy=-max(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards NEE
                        d_xx=max(d);
                        d_xy=min(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards NEE
                        d_xx=-max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards SEE
                        d_xx=max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards SEE
                        d_xx=-max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards SSE
                        d_xx=min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards SSE
                        d_xx=-min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards SSW
                        d_xx=-min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards SSW
                        d_xx=min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards SWW
                        d_xx=-max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards SWW
                        d_xx=max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards NWW
                        d_xx=-max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards NWW
                        d_xx=max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards NNW
                        d_xx=-min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards NNW
                        d_xx=min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==1 % up, towards N
                        d_xx=0;
                        d_xy=max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==0 % down, towards N
                        d_xx=0;
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==1 % up, towards S
                        d_xx=0;
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==0 % down, towards S
                        d_xx=0;
                        d_xy=max(d);
                    elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==1 % up, towards E
                        d_xx=max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==0 % down, towards E
                        d_xx=-max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==1 % up, towards W
                        d_xx=-max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==0 % down, towards W
                        d_xx=max(d);
                        d_xy=0;
                    else
                        d_xx=NaN;
                        d_xy=NaN;
                    end

                    % coordinate transformation
                    temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii)+d_xx trh.y(ii)+d_xy; trh.x(ii+anz2)+d_xx trh.y(ii+anz2)+d_xy]);

                    trh.xn(ii)=temp(1);
                    trh.yn(ii)=temp(2);
                    trh.zn(ii)=trh.z(ii)-hnew+d_z; % correct height
                else % no topography present
                    temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
                    trh.xn(ii)=temp(1);
                    trh.yn(ii)=temp(2);
                end
            end
            % calculation for the last (anz2) few traces -> flip
            % coordinates and start profile from back
            n=length(trh.x);
            tempx_flip=fliplr(trh.x);
            tempy_flip=fliplr(trh.y);
            tempz_flip=fliplr(trh.z);
            tempdx=-offsetGPS_X;
            tempdy=-offsetGPS_Y;
            for ii=1:anz2
                dist=sqrt((tempx_flip(ii)-tempx_flip(ii+anz2))^2+(tempy_flip(ii)-tempy_flip(ii+anz2))^2); % distance between two points anz2-traces away
                dist_xy=[tempx_flip(ii+anz2)-tempx_flip(ii) tempy_flip(ii+anz2)-tempy_flip(ii)]; % differences in x and y direction
                % correct for GPS height and associated xy-shift if
                % topography is present
                if any(tempz_flip~=0)
                    alpha=atand(abs(tempz_flip(ii+anz2)-tempz_flip(ii))/dist); % slope angle in degree
                    hnew=h_GPS/cosd(alpha); % new height of GPS above ground at measured point
                    b=hnew*sind(alpha); % shift along surface to correct point
                    d_x=b*cosd(alpha); % horizontal shift of point in profile direction
                    d_z=b*sind(alpha); % vertical shift of point in profile direction

                    % angle of profile direction towards north
                    if tempx_flip(ii)==tempx_flip(ii+anz2) % north/south
                        beta=0;
                    elseif tempy_flip(ii)==tempy_flip(ii+anz2) % west/east
                        beta=90;
                    else
                        beta=atand(abs(tempx_flip(ii+anz2)-tempx_flip(ii))/abs(tempy_flip(ii+anz2)-tempy_flip(ii)));
                    end
                    % split horizontal shift into x and y component:
                    d=[d_x*cosd(beta) d_x*sind(beta)];
                    % differentiate between different cases
                    updownflag=tempz_flip(ii)<tempz_flip(ii+anz2); % up=1, down=0
                    if all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards NNE
                        d_xx=min(d);
                        d_xy=max(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards NNE
                        d_xx=-min(d);
                        d_xy=-max(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards NEE
                        d_xx=max(d);
                        d_xy=min(d);
                    elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards NEE
                        d_xx=-max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards SEE
                        d_xx=max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards SEE
                        d_xx=-max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards SSE
                        d_xx=min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards SSE
                        d_xx=-min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards SSW
                        d_xx=-min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards SSW
                        d_xx=min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards SWW
                        d_xx=-max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards SWW
                        d_xx=max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards NWW
                        d_xx=-max(d);
                        d_xy=min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards NWW
                        d_xx=max(d);
                        d_xy=-min(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards NNW
                        d_xx=-min(d);
                        d_xy=max(d);
                    elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards NNW
                        d_xx=min(d);
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==1 % up, towards N
                        d_xx=0;
                        d_xy=max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==0 % down, towards N
                        d_xx=0;
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==1 % up, towards S
                        d_xx=0;
                        d_xy=-max(d);
                    elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==0 % down, towards S
                        d_xx=0;
                        d_xy=max(d);
                    elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==1 % up, towards E
                        d_xx=max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==0 % down, towards E
                        d_xx=-max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==1 % up, towards W
                        d_xx=-max(d);
                        d_xy=0;
                    elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==0 % down, towards W
                        d_xx=max(d);
                        d_xy=0;
                    else
                        d_xx=NaN;
                        d_xy=NaN;
                    end

                    % coordinate transformation
                    temp=helmert([tempdx tempdy],[0 0; 0 dist],[tempx_flip(ii)+d_xx tempy_flip(ii)+d_xy; tempx_flip(ii+anz2)+d_xx tempy_flip(ii+anz2)+d_xy]);

                    trh.xn(n-ii+1)=temp(1);
                    trh.yn(n-ii+1)=temp(2);
                    trh.zn(n-ii+1)=trh.z(n-ii+1)-hnew+d_z; % correct height
                else % no topography present
                    temp=helmert([tempdx tempdy],[0 0; 0 dist],[tempx_flip(ii) tempy_flip(ii); tempx_flip(ii+anz2) tempy_flip(ii+anz2)]);
                    trh.xn(n-ii+1)=temp(1);
                    trh.yn(n-ii+1)=temp(2);
                end
            end

            % set new/corrected coordinates in original fields in trace
            % header
            trh.x=trh.xn;
            trh.y=trh.yn;
            trh.z=trh.zn;
        end
        
        
        % if readdzt_all was used no field time is created
        if ~isfield(trh,'time')
            trh.time = zeros(size(trh.x));
            trh.quality = NaN(size(trh.x)); % NaN because zero has a definition: (0=Fix not valid)
        end

        if ~isfield(trh,'quality')
            trh.quality = NaN(size(trh.x)); % NaN because zero has a definition: (0=Fix not valid)
        end
        
        if h.nchan==2
            if i==1
                radargrams=[{data(:,trh.channum==1)}];
                radargrams=[radargrams; {data(:,trh.channum==2)}];
                global_coords=[{[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                global_coords=[global_coords; {[trh.x(trh.channum==2)' trh.y(trh.channum==2)' trh.z(trh.channum==2)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                xtemp=trh.x(trh.channum==2);
                ytemp=trh.y(trh.channum==2);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                time = [{trh.time(trh.channum==1)}];
                time = [time; {trh.time(trh.channum==2)}];
                quality = [{trh.quality(trh.channum==1)}];
                quality = [quality; {trh.quality(trh.channum==2)}];
                if ~all(trh.mark==0) % marker present
                    marker=[{trh.mark(trh.channum==1)}];
                    marker=[marker; {trh.mark(trh.channum==1)}]; % marker are only present in first channel
                end
                listrad.name=[{name{i}}; {name{i}}];
                listrad.chan=[1; 2];
            else
                radargrams=[radargrams; {data(:,trh.channum==1)}];
                radargrams=[radargrams; {data(:,trh.channum==2)}];
                global_coords=[global_coords; {[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                global_coords=[global_coords; {[trh.x(trh.channum==2)' trh.y(trh.channum==2)' trh.z(trh.channum==2)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];     
                end
                xtemp=trh.x(trh.channum==2);
                ytemp=trh.y(trh.channum==2);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];  
                end
                time = [time; {trh.time(trh.channum==1)}];
                time = [time; {trh.time(trh.channum==2)}];
                quality = [quality; {trh.quality(trh.channum==1)}];
                quality = [quality; {trh.quality(trh.channum==2)}];
                if ~all(trh.mark==0) % marker present
                    marker=[marker; {trh.mark(trh.channum==1)}];
                    marker=[marker; {trh.mark(trh.channum==1)}];
                end
                listrad.name=[listrad.name; {name{i}}; {name{i}}];
                listrad.chan=[listrad.chan; 1; 2];
            end
            t=double(0:h.dt:h.dt*(h.ns-1));
            if max(t)<1 % for DF-Antenna is dt in s 
                t=t.*1e9; % convert to ns
            end

            if isfield(h,'dt2')
                t2=double(0:h.dt2:h.dt2*(h.ns-1));
                if max(t2)<1
                    t2=t2.*1e9; % convert to ns for DF antenna
                end
            end
            
            chan(anz)=1;
            chan(anz+1)=2;
            anz=anz+2;
        else % 1 channel
            if i==1
                radargrams=[{data(:,trh.channum==1)}];
                global_coords=[{[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                
                time = [{trh.time(trh.channum==1)}];
                quality = [{trh.quality(trh.channum==1)}];

                if ~all(trh.mark==0) % marker present
                    marker=[{trh.mark(trh.channum==1)}];
                end
                
                listrad.name=[{name{i}}];
                listrad.chan=[1];
            else
                radargrams=[radargrams; {data(:,trh.channum==1)}];
                global_coords=[global_coords; {[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];    
                end
                
                time = [time; {trh.time(trh.channum==1)}];
                quality = [quality; {trh.quality(trh.channum==1)}];
                if ~all(trh.mark==0) % marker present
                    marker=[marker; {trh.mark(trh.channum==1)}];
                end
                
                listrad.name=[listrad.name; {name{i}}];
                listrad.chan=[listrad.chan; 1];
            end
            
            t=double(0:h.dt:h.dt*(h.ns-1));
            if max(t)<1 % for DF-Antenna is dt in s 
                t=t.*1e9; % convert to ns
            end
            chan(anz)=1;
            anz=anz+1;
        end
       % collect headers in cell
       headers{i} = h;

       clear trh;

end
%------------------------------ WRITE DATA --------------------------------

if export2mat==1
   
    % delete empty cells
    radargrams(~cellfun('isempty',radargrams));
    
    disp('Start saving of mat-files...')
    if matlab==1
        % save all profiles in one variable
        if isfield(h,'dt2') % split into different folders for different time ranges (DF antenna)
            % temporary variables for splitting:
            rtemp=radargrams;
            xtemp=x;
            ctemp=global_coords;
            mtemp=marker;
            ttemp=t;
            timetemp = time;
            qualitytemp = quality;
            clear global_coords;
            clear x;
            clear radargrams;
            clear marker;
            clear time;
            clear quality;
            
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch1'),'dir')
                mkdir(fullfile(pfad,'mat_ch1'));
            end
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch2'),'dir')
                mkdir(fullfile(pfad,'mat_ch2'));
            end
            
            % Channel1 
            anz=1;
            for i=1:length(chan)
                if chan(i)==1
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    quality{anz}=qualitytemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            save(fullfile(pfad,'mat_ch1','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat_ch1','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat_ch1','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat_ch1','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat_ch1','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat_ch1','time.mat'),'time','-v7.3');
            save(fullfile(pfad,'mat_ch1','quality.mat'),'quality','-v7.3');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch1','marker.mat'),'marker','-v7.3');
            end
            fid=fopen(fullfile(pfad,'mat_ch1','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            n=1;
            for i=1:length(listrad.name)
                if chan(i)==1
                    fprintf(fid,'%d\t',n);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                    n=n+1;
                end
            end
            fclose(fid);
            
            % Channel2
            anz=1;
            for i=1:length(chan)
                if chan(i)==2
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    quality{anz}=qualitytemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            t=t2;
            save(fullfile(pfad,'mat_ch2','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat_ch2','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat_ch2','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat_ch2','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat_ch2','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat_ch2','time.mat'),'time','-v7.3');
            save(fullfile(pfad,'mat_ch2','quality.mat'),'quality','-v7.3');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch2','marker.mat'),'marker','-v7.3');
            end
            fid=fopen(fullfile(pfad,'mat_ch2','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            n=1;
            for i=1:length(listrad.name)
                if chan(i)==2
                    fprintf(fid,'%d\t',n);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                    n=n+1;
                end
            end
            fclose(fid);
            
            % set variables back to original (for sgy export)
            global_coords=ctemp;
            t=ttemp;
            marker=mtemp;
            radargrams=rtemp;
            x=xtemp;
            
        else  % for one channel/not DF antenna
            % make new folder
            if ~exist(fullfile(pfad,'mat'),'dir')
                mkdir(fullfile(pfad,'mat'));
            end
            save(fullfile(pfad,'mat','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat','time.mat'),'time','-v7.3');
            save(fullfile(pfad,'mat','quality.mat'),'quality','-v7.3');
            if ~isempty(marker)
                save(fullfile(pfad,'mat','marker.mat'),'marker','-v7.3');
            end
            
            fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                fprintf(fid,'%d\t',i);
                fprintf(fid,'%s\t',listrad.name{i});
                fprintf(fid,'%d\n',listrad.chan(i));
            end
            fclose(fid);
        end
    else  % Octave
        % save all profiles in one variable
        if isfield(h,'dt2') % split into different folders for different time ranges (DF antenna)
            % temporary variables for splitting:
            rtemp=radargrams;
            xtemp=x;
            ctemp=global_coords;
            mtemp=marker;
            ttemp=t;
            timetemp = time;
            qualitytemp = quality;
            clear global_coords;
            clear x;
            clear radargrams;
            clear marker;
            clear time;
            clear quality;

            % make new folder
            if ~exist(fullfile(pfad,'mat_ch1'),'dir')
                mkdir(fullfile(pfad,'mat_ch1'));
            end
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch2'),'dir')
                mkdir(fullfile(pfad,'mat_ch2'));
            end
            
            % Channel1 
            anz=1;
            for i=1:length(chan)
                if chan(i)==1
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    quality{anz}=qualitytemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            save(fullfile(pfad,'mat_ch1','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat_ch1','t.mat'),'t');
            save(fullfile(pfad,'mat_ch1','x.mat'),'x');
            save(fullfile(pfad,'mat_ch1','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat_ch1','h.mat'),'headers');
            save(fullfile(pfad,'mat_ch1','time.mat'),'time');
            save(fullfile(pfad,'mat_ch1','quality.mat'),'quality');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch1','marker.mat'),'marker');
            end
            fid=fopen(fullfile(pfad,'mat_ch1','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==1
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % Channel2
            anz=1;
            for i=1:length(chan)
                if chan(i)==2
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    quality{anz}=qualitytemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            t=t2;
            save(fullfile(pfad,'mat_ch2','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat_ch2','t.mat'),'t');
            save(fullfile(pfad,'mat_ch2','x.mat'),'x');
            save(fullfile(pfad,'mat_ch2','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat_ch2','h.mat'),'headers');
            save(fullfile(pfad,'mat_ch2','time.mat'),'time');
            save(fullfile(pfad,'mat_ch2','quality.mat'),'quality');

            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch2','marker.mat'),'marker');
            end
            fid=fopen(fullfile(pfad,'mat_ch2','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==2
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % set variables back to original (for sgy export)
            global_coords=ctemp;
            t=ttemp;
            marker=mtemp;
            radargrams=rtemp;
            x=xtemp;
            
        else  % for one channel/not DF antenna
            % make new folder
            if ~exist(fullfile(pfad,'mat'),'dir')
                mkdir(fullfile(pfad,'mat'));
            end
            save(fullfile(pfad,'mat','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat','t.mat'),'t');
            save(fullfile(pfad,'mat','x.mat'),'x');
            save(fullfile(pfad,'mat','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat','h.mat'),'headers');
            save(fullfile(pfad,'mat','time.mat'),'time');
            save(fullfile(pfad,'mat','quality.mat'),'quality');
            if ~isempty(marker)
                save(fullfile(pfad,'mat','marker.mat'),'marker');
            end
            
            fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                fprintf(fid,'%d\t',i);
                fprintf(fid,'%s\t',listrad.name{i});
                fprintf(fid,'%d\n',listrad.chan(i));
            end
            fclose(fid);
        end
    end
   
    disp('   Finished!')
end


if export2segy==1
    disp('Start saving as sgy...')
    % make new folder
    if ~exist(fullfile(pfad,'SEGY'),'dir')
        mkdir(fullfile(pfad,'SEGY'));
    end

    if isfield(h,'dt2') % is DF antenna
        h.dt=h.dt*1e9; % convert to ns
        h.dt2=h.dt2*1e9;
        for i=1:length(headers)
            headers{i}.dt=headers{i}.dt*1e9;
            headers{i}.dt2=headers{i}.dt2*1e9;
        end
    end
    
    for i=1:length(radargrams)
        if length(global_coords{i}(1,:))<3
            global_coords{i}(:,3)=zeros(size(global_coords{i}(:,1))); % if no topography present, set to zero
        end
        if chan(i)==2 && isfield(h,'dt2')
            export2sgy2D(radargrams{i},headers{i/2}.dt2,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);
        else
            if isfield(h,'dt2')
                export2sgy2D(radargrams{i},headers{(i-1)/2+1}.dt,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);
            else
                export2sgy2D(radargrams{i},headers{i}.dt,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);
            end
        end
    end
    
    disp('   Finished!')
end

% Exporting settings:
disp('Exported settings_DZTconvert.txt!')

fid=fopen(fullfile(pfad,'mat','settings_DZTconvert.txt'),'wt');
if fid==-1
    fid=fopen(fullfile(pfad,'mat_ch1','settings_DZTconvert.txt'),'wt');
end
fprintf(fid,'Equipment:\n');
fprintf(fid,'  app=%s\n\n',app);    % Equipment: SIR20 / SIR30 / SIR3000 / SIR4000 / Tablet / UtilityScan (UtilityScan with DF antenna only!)

fprintf(fid,'Options for GNSS:\n');
fprintf(fid,'  convert2utm=%d\n',convert2utm); % convert WGS Lat/Long to UTM (=1 if measured with Stonex-GPS)
fprintf(fid,'  zone=%d\n',zone); % if convert2utm==1 -> give UTM-zone
fprintf(fid,'  offsetGPS_X=%f\n',offsetGPS_X); % [m] Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
fprintf(fid,'  offsetGPS_Y=%f\n',offsetGPS_Y); % [m] Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)
fprintf(fid,'  h_GPS=%f\n\n',h_GPS); % [m] height of GPS/prism above ground

% smoothing of GPS-coordinates before applying antenna-offsets:
fprintf(fid,'Smoothing of coordinates:\n');
fprintf(fid,'  smooth_coords=%d\n',smooth_coords); % yes=1, no=0
fprintf(fid,'  numsamp_smooth=%d\n\n',numsamp_smooth); % number of samples for moving median smoothing

% Options for calculating inline coordinates for each trace:
fprintf(fid,'Profile coordinate calculation:\n');
fprintf(fid,'  coords_opt=%d\n\n',coords_opt);   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
                % =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)

% Attention: still experimental! If you think that the output radargrams
% are wrong, set to 0!
fprintf(fid,'Remove outliers (still experimental):\n');
fprintf(fid,'  removeOutliers=%d\n\n',removeOutliers); % do you want to remove coordinate outliers? 
    % 0= no, use raw coordinates
    % 1= in middle of profile
    % 2= at the end/beginning of profile

% remove start/end traces of profiles
fprintf(fid,'Remove start/end of profiles:\n');
fprintf(fid,'  removeStartEnd=%d\n\n',removeStartEnd); % =1: yes, remove start and end traces of profiles (at same position)
                    % =0: No. Keep all traces.

% Export to other formats
fprintf(fid,'Export to other formats:\n');
fprintf(fid,'  export2mat=%d\n',export2mat); % export to Multichannel-GPR format for radargrams (mat-files)
fprintf(fid,'  export2segy=%d\n',export2segy); % export all radargrams as segy-files
fprintf(fid,'  constoff=%d\n\n',constoff); % if=1: a constant coordinate offset will be subtracted and coordinates will be in mm accuracy in segy file (offsets will be saved in Inline3D (x) and Crossline3D (y))
fclose(fid);


% set original path
path(oldpath);
