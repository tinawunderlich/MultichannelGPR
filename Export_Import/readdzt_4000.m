function [data,trh,h]=readdzt_4000(dztfile,plotten,removeOutliers,zone)


%%% [data,trh,h]=readdzt(dztfile,plotten)
%
% dztfile: Name der dzt-Datei
% plotten: wenn =1 -> Bild wird geplottet, sonst =0
% removeOutliers: (optional) flag if outliers should be deleted with different
% options: 0=no, 1=in middle of profile, 2=at end/beginning of profile
% zone: UTM zone (optional) for transforming into UTM coordinates
%
% Output:
% data: Datensamples in einer Matrix
% trh: Traceheader f�r jede Spur
% h: Header f�r die gesamte Datei
% falls zwei Kan�le: alle Daten zusammen in data, aber im Header gibt es
% Kanalnummer
%
% Dr. Tina Wunderlich, CAU, November 2017-August 2022
%%%

if nargin==2
    removeOutliers=0;
    zone=0;
elseif nargin==3
    zone=0;
end

test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end


% name zerlegen
[pathstr,fname,ext] = fileparts(dztfile);


% Koordinatendatei
if exist([pathstr,'/',fname,'.DZG'])    % falls koordinatendatei vorhanden

    if ~exist([pathstr,'/',fname,'_backup.DZG'])  % if backup already exists, some changes have been made before, so do not do the changes again.
        % PART 1: check DZG-file if corrupted (sometimes if measured with Stonex-GPS)
        copyfile([pathstr,'/',fname,'.DZG'],[pathstr,'/',fname,'_backup.DZG']); % Save backup file

        fileIDorig = fopen([pathstr,'/',fname,'_backup.DZG'],'r');
        fileIDnew = fopen([pathstr,'/',fname,'.DZG'],'wt');

        dataorig = textscan(fileIDorig,'%s','Delimiter','\n');
        dataorig = dataorig{1};
        ndata = length(dataorig);
        for j = 1:ndata

            tmpline = strsplit(dataorig{j},',','CollapseDelimiters',false);
            nentry = length(tmpline);
            if strcmp(tmpline{1},'$GSSIS') && nentry == 3 && ...
                    ~isempty(str2double(tmpline{2})) && ...
                    ~isempty(str2double(tmpline{3}))
                fprintf(fileIDnew,'%s\n',dataorig{j});
            elseif strcmp(tmpline{1},'$GPGGA') && nentry == 15
                fprintf(fileIDnew,'%s\n',dataorig{j});
            elseif isempty(tmpline{1})
                fprintf(fileIDnew,'%s','');
            end
        end
        fclose(fileIDorig);
        fclose(fileIDnew);
        s1=dir([pathstr,'/',fname,'_backup.DZG']);
        s2=dir([pathstr,'/',fname,'.DZG']);
        if s1.bytes==s2.bytes  % if same size, delete backup, because no changes
            delete([pathstr,'/',fname,'_backup.DZG']);
        end
        % file.DZG is now without GPSerrors
    end

    % PART 2: Read DZG-file
    fid=fopen([pathstr,'/',fname,'.DZG']);
    temp=textscan(fid,'%s');
    fclose(fid);

    if ~isempty(temp{1})
        count = 1;
        % loop over lines of file
        for j=1:length(temp{1})

            % the first line contains the trace number for the coordinate
            if strfind(temp{1}{j},'$GSSIS')
                if matlab
                    GSSILine=split(temp{1}{j},',');
                else
                    GSSILine=strsplit(temp{1}{j},',');
                end


                if j < length(temp{1})
                    % the coordinates are found in the second line
                    if strfind(temp{1}{j+1},'$GPGGA')
                        if matlab
                            GPGGALine=split(temp{1}{j+1},',');
                        else
                            GPGGALine=strsplit(temp{1}{j+1},',');
                        end
                        % extract the coordinate information
                        if strcmp(GPGGALine{6},'E')
                            x(count) = str2double(GPGGALine{5});
                        else
                            x(count) = -str2double(GPGGALine{5}); % West of Greenwich -> negative
                        end
                        y(count) = str2double(GPGGALine{3});
                        z(count) = str2double(GPGGALine{10});

                        % extract the time of each sample
                        if matlab
                            t(count) = datetime(GPGGALine{2},'InputFormat','HHmmss.SS');
                        end

                        % extract the gps quality
                        q(count) = str2double(GPGGALine{7});

                        % extract the trace number
                        trn(count)=str2double(GSSILine{2});
                        count = count + 1;
                        
                    end
                end
            end

        end
        [trnUnique, ia, ~] = unique(trn);

        % save the data into cell structures
        trnum = trnUnique;
        if matlab
            time = t(ia);
        else
            time = zeros(size(trnum));
        end

        quality = q(ia);

        xyz = [ trnum', x(ia)', y(ia)', z(ia)' ];
    else
        xyz=[];
    end


    % PART 3: Read Markerdatei (DZX)
    fid=fopen([pathstr,'/',fname,'.DZX']);
    temp=textscan(fid,'%s');
    DF=0;
    for i=1:length(temp{1})
        if strfind(temp{1}{i},'<scanRange>')
            scanrange(1)=str2num(temp{1}{i}(12));
            for j=14:19
                if ~isempty(str2num(temp{1}{i}(14:j)))
                    scanrange(2)=str2num(temp{1}{i}(14:j));
                end
            end
        end
        if strfind(temp{1}{i},'<systemMode>')
            DFtest=temp{1}{i}(13:14);
            if strcmp(DFtest,'DF')
                DF=1;
            end
        end
        if strfind(temp{1}{i},'<depthRange>')
            dr=str2num(temp{1}{i}(13:18));
        end
        if strfind(temp{1}{i},'<depthRangeCh2>')
            dr2=str2num(temp{1}{i}(16:21));
        end
    end
    fclose(fid);

    if ~isempty(xyz)
        mark=zeros(size(xyz(:,1)))';
    else
        mark=[];
    end
elseif exist(fullfile(pathstr,[fname,'.DZX']))
    % Markerdatei
    fid=fopen(fullfile(pathstr,[fname,'.DZX']));
    temp=textscan(fid,'%s');
    DF=0;
    anz=1;
    anzc=1;
    for i=1:length(temp{1})
        if strfind(temp{1}{i},'<scanRange>')
            scanrange(1)=str2num(temp{1}{i}(12));
            for j=14:19
                if ~isempty(str2num(temp{1}{i}(14:j)))
                    scanrange(2)=str2num(temp{1}{i}(14:j));
                end
            end
        end
        if strfind(temp{1}{i},'<systemMode>')
            DFtest=temp{1}{i}(13:14);
            if strcmp(DFtest,'DF')
                DF=1;
            end
        end
        if strfind(temp{1}{i},'<depthRange>')
            dr=str2num(temp{1}{i}(13:18));
        end
        if strfind(temp{1}{i},'<depthRangeCh2>')
            dr2=str2num(temp{1}{i}(16:21));
        end
        if strfind(temp{1}{i},'<scan>')
            for j=7:12
                if ~isempty(str2num(temp{1}{i}(7:j)))
                    mark(anz)=str2num(temp{1}{i}(7:j));
                end
            end
            anz=anz+1;
        end
        if strfind(temp{1}{i},'<localCoords>')
            for j=13:21
                if ~isempty(str2num(temp{1}{i}(13:j)))
                    coords(anzc)=str2num(temp{1}{i}(13:j));
                end
            end
            anzc=anzc+1;
        end
        if strfind(temp{1}{i},'<mark>User</mark>')
            marktype{anz-1}='User';
        elseif strfind(temp{1}{i},'<mark>Distance</mark>')
            marktype{anz-1}='Distance';
        else
            marktype{1}='NA';
        end
    end
    fclose(fid);
    mark(length(marktype)+1:end)=[];

    xyz=zeros(length(mark),4);
    xyz(:,1)=mark';
else
    xyz=[];
    DF=0;
    mark=[];
end
if ~exist('marktype','var')
    marktype='NA';
end

% PART 4: read radar data (DZT)
fid=fopen(dztfile,'r','l');

test=fread(fid,1,'int16','l');

%%% Radar Header auslesen:

% Apparatur-Tag
SIR10=0;
SIR20=0;
SIR30=0;
SIR3000=0;
SIR4000=0;
% dann die entsprechende Apparatur raussuchen
if test(1)==2047
    SIR30=1;
    SIR20=1;
    SIR4000=1;
elseif test(1)==255
    SIR10=1;
elseif test(1)==1792
    SIR3000=1;
end
rh_tag=test(1);

MINHEADSIZE=1024;

%---------------------------------------------------------
test=fread(fid,3,'uint16','l');   % weiter lesen...

rh_data=test(1);    % Offset to data from beginning of file
rh_nsamp=test(2);   % Samples per Scan
rh_bit=test(3);     % Bits per data word (16 = uint16 f�r sp�teres Dateneinlesen)

%---------------------------------------------------------
test=fread(fid,1,'short','l');  % weiter lesen...

rh_zero=test(1);    % Binary offset (0x80 bei 8 bit, 0x8000 bei 16 bit)

%---------------------------------------------------------
test=fread(fid,5,'single=>float','l');   % weiter lesen...

rh_sps=test(1); % Scans per Second
rh_mpm=test(2); % Scans per meter
rh_position=test(3);    % Position start in ns (Offset 1. Sample zum Nullpunkt)
rh_range=test(5);   % Range in ns

%---------------------------------------------------------
test=fread(fid,1,'short=>short','l');  % weiter lesen...

rh_npass=test(1);   % number of passes for 2D-Files (?)

%---------------------------------------------------------
test=fread(fid,2,'uint32=>uint32','l');  % weiter lesen...

rh_create=test(1);  % Datum created
rh_modif=test(2);   % Datum modified

% add date to datetime object: time
dt = readTime(rh_create);
time(:).Year = dt.Year;
time(:).Month = dt.Month;
time(:).Day = dt.Day;

%---------------------------------------------------------
test=fread(fid,7,'uint16','l');  % weiter lesen...

rh_rgain=test(1);   % Offset to range gain function
rh_nrgain=test(2);  % Size of range gain function
rh_text=test(3);    % Offset to text
rh_ntext=test(4);   % Size of text
rh_proc=test(5);    % Processing history
rh_nproc=test(6);
rh_nchan=test(7);   % Number of channels

%---------------------------------------------------------
test=fread(fid,3,'single','l');  % weiter lesen...

rh_epsr=test(1);    % Average permittivity
rh_top=test(2);     % Position in meters
rh_depth=test(3);   % Range in meters

%---------------------------------------------------------
test=fread(fid,1,'double','l');  % weiter lesen...

rh_coordX=test(1);  % x-Coord

%---------------------------------------------------------
test=fread(fid,1,'single','l');  % weiter lesen...

rhf_servo_level=test(1);    % Gain servo level

%---------------------------------------------------------
test=fread(fid,4,'uchar=>uchar','l');  % weiter lesen...

rh_accomp=test(4);  % ?

%---------------------------------------------------------
test=fread(fid,1,'int16','l');  % weiter lesen...

rh_sconfig=test(1); % setup configuration number

%---------------------------------------------------------
test=fread(fid,2,'uint16','l');  % weiter lesen...

rh_spp=test(1); % scans per pass
rh_linenum=test(2); % line number

%---------------------------------------------------------
test=fread(fid,1,'double','l');  % weiter lesen...

rh_coordY=test(1);  % Y-coordinate

%---------------------------------------------------------
test=fread(fid,2,'uchar=>uchar','l');  % weiter lesen...

rh_lineorder=test(1);   % Line order
rh_dtype=test(2);   % ?

%---------------------------------------------------------
test=fread(fid,14,'char=>char','l');  % weiter lesen...

rh_antname=test;    % Antennenname

%---------------------------------------------------------
test=fread(fid,2,'uchar=>uchar','l');  % weiter lesen...

rh_pass0TX=test(1); % Active transmit mask
rh_version=test(2); % 1: no GPS, 2: PPS

%---------------------------------------------------------
test=fread(fid,12,'char=>char','l');  % weiter lesen...

rh_name=test;   % File-Name

%---------------------------------------------------------
test=fread(fid,1,'int16','l');  % weiter lesen...

rh_chksum=test(1);  % checksum for header

%---------------------------------------------------------
% test=fread(fid,1024-128-2*(9),'uchar=>uchar','l');  % weiter lesen...
%
% variable=test;  % rh_rgain, rh_text und rh_proc

%---------------------------------------------------------
% Springen auf Anfang der Daten:
if rh_data<MINHEADSIZE
    offset=MINHEADSIZE*rh_data;
else
    offset=MINHEADSIZE*rh_nchan;
end

%fseek(fid,rh_data*rh_nchan+rh_nsamp*rh_bit*(rh_zero+1),'bof');
fseek(fid,offset,'bof');
if rh_bit==16
    test=fread(fid,'uint16','l'); % Rest einlesen zum gucken wie viele Spuren
elseif rh_bit==32
    test=fread(fid,'int32','l'); % Rest einlesen zum gucken wie viele Spuren
end

rh_ntraces=ceil(length(test)/rh_nsamp);  % Anzahl Spuren (richtig!)

%---------------------------------------------------------
%fseek(fid,rh_data*rh_nchan+rh_nsamp*rh_bit*2,'bof'); % auf Datenanfang zurück setzen
%fseek(fid,rh_data*rh_nchan+rh_nsamp*rh_bit*(rh_zero+1),'bof');
fseek(fid,offset,'bof');

data=zeros(rh_nsamp,rh_ntraces*rh_nchan); % Datenmatrix initialisieren

% Daten auslesen
datatemp=fread(fid,[rh_nsamp,rh_ntraces*rh_nchan],'int32','l');
if any(strcmp(marktype,'NA')) && isempty(mark)
    mark=find(datatemp(2,:)~=0); % markers are non-zero numbers in second row
    mark(1)=[]; % delete marker on first trace
elseif any(strcmp(marktype,'NA')) && length(mark)==length(marktype)
    if mark(strcmp(marktype,'NA'))==1 % if marker on first trace with 'NA' -> not real marker
        mark(1)=[];
        marktype(1)=[];
    else
        marktype(strcmp(marktype,'NA'))={'User'}; % valid marker -> set correct type
    end
end
% ersten zwei Zeilen Null setzen
datatemp(1:2,:)=0;
% /1000 und Mittelwert abziehen
datatemp=datatemp./1000;
datatemp=datatemp-mean(datatemp(:));


if rh_nchan==1
    data=datatemp;
else
    data=datatemp(:,1:2:end);
    data=[data datatemp(:,2:2:end)]; % h�nge Kanal 2 hinten dran
end


% falls am Anfang oder Ende Leerspuren sind -> l�schen
weg=find(data(150,:)==0);
data(:,weg)=[];


% delete line if coordinates nan
if isempty(xyz)
    xyz=zeros(length(data(1,:)),3);
end
xyz(isnan(xyz(:,2)) | isnan(xyz(:,3)),:)=[];

if zone~=0 && ~all(all(xyz(:,2:3)==0)) % convert to UTM
    lontemp=num2str(xyz(:,2),'%.8f');
    lattemp=num2str(xyz(:,3),'%.8f');
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
    [xneu,yneu]=wgs2utm(lat,lon,zone,'N');
    xyz(:,2)=xneu;
    xyz(:,3)=yneu;
end


% Koordinaten interpolieren und Stillstandspuren l�schen
if ~isempty(xyz) && ~all(xyz(:,2)==0)
    % check for large coordinate jumps in the middle of the profile and delete them:
    if removeOutliers==1

%         figure
%         subplot(2,2,[1 3])
%         plot(xyz(:,2),xyz(:,3),'*')
%         hold on
%         xlabel('x')
%         ylabel('y')
%         subplot(2,2,2)
%         plot(xyz(:,1),xyz(:,2),'*')
%         hold on
%         xlabel('Trace number')
%         ylabel('x')
%         subplot(2,2,4)
%         plot(xyz(:,1),xyz(:,3),'*')
%         hold on
%         xlabel('Trace number')
%         ylabel('y')
        
        % compare smoothed and input coordinates and find points with
        % larger deviations = outliers:
        diffx=abs(xyz(:,2)-smooth(xyz(:,2),5));
        diffy=abs(xyz(:,3)-smooth(xyz(:,3),5));
        % make histogram with 10 bins:
        [hx,distx]=hist(diffx,10);
        [hy,disty]=hist(diffy,10);
        % calculate parts of all differences
        prozx=cumsum(hx)./sum(hx);
        prozy=cumsum(hy)./sum(hy);
        % go through numbers and find the one when it stays constant for
        % two entries:
        tempx=prozx(1);
        tempy=prozy(1);
        i=2;
        while i<=length(prozx) && (tempx~=prozx(i) && tempy~=prozy(i))
            tempx=prozx(i);
            tempy=prozy(i);
            i=i+1;
        end
        % find all differences larger than distx(i) and disty(i) and delete
        % them:
        weg=(diffx>distx(min(i,10)) | diffy>disty(min(i,10)));
        xyz(weg,:)=[];


%         subplot(2,2,[1 3])
%         plot(xyz(:,2),xyz(:,3),'*')
%         legend('original','filtered')
%         subplot(2,2,2)
%         plot(xyz(:,1),xyz(:,2),'*')
%         subplot(2,2,4)
%         plot(xyz(:,1),xyz(:,3),'*')
        
    end
    
    % remove outliers at end/beginning of profile
    if removeOutliers==2
        
%         figure
%         subplot(2,2,[1 3])
%         plot(xyz(:,2),xyz(:,3),'*')
%         hold on
%         xlabel('x')
%         ylabel('y')
%         subplot(2,2,2)
%         plot(xyz(:,1),xyz(:,2),'*')
%         hold on
%         xlabel('Trace number')
%         ylabel('x')
%         subplot(2,2,4)
%         plot(xyz(:,1),xyz(:,3),'*')
%         hold on
%         xlabel('Trace number')
%         ylabel('y')
        
        % make linear fit through coordinate points for x and y separately
        p1=polyfit(xyz(:,1),xyz(:,2),1);
        p2=polyfit(xyz(:,1),xyz(:,3),1);
        
        % calculate differences between points and line
        diffx=abs(xyz(:,2)-polyval(p1,xyz(:,1)));
        diffy=abs(xyz(:,3)-polyval(p2,xyz(:,1)));
        
        % make histogram with 10 bins:
        [hx,distx]=hist(diffx,10);
        [hy,disty]=hist(diffy,10);
        % calculate parts of all differences
        prozx=cumsum(hx)./sum(hx);
        prozy=cumsum(hy)./sum(hy);
        % go through numbers and find the one when it stays constant for
        % two entries:
        tempx=prozx(1);
        tempy=prozy(1);
        i=2;
        while tempx~=prozx(i) && tempy~=prozy(i)
            tempx=prozx(i);
            tempy=prozy(i);
            i=i+1;
        end
        % find all differences larger than distx(i) and disty(i) and delete
        % them:
        weg=(diffx>distx(i) | diffy>disty(i));
        xyz(weg,:)=[];
        
%         subplot(2,2,[1 3])
%         plot(xyz(:,2),xyz(:,3),'*')
%         legend('original','filtered')
%         subplot(2,2,2)
%         plot(xyz(:,1),xyz(:,2),'*')
%         subplot(2,2,4)
%         plot(xyz(:,1),xyz(:,3),'*')
    end

    if size(quality,2) ~= size(xyz,1) && exist("weg",'var')
        quality(weg) = [];
    end

    trnum1=1:length(data(1,:))/rh_nchan;
    [a,b]=unique(xyz(:,1));
    xyz=xyz(b,:); % delete non-unique points
    quality=quality(b);


    % remove outliers 3
    if removeOutliers==3
        % window size of median filter: 2 percent of all traces
        wind = round(max(trnum)*0.02);
        [xyz(:,2),xyz(:,3),xyz(:,4)] = rmOutliers(xyz(:,2),xyz(:,3),xyz(:,4),wind);
    end
    
    % smooth coordinates over 15 sample
    sx=smooth(xyz(:,2),15);
    sy=smooth(xyz(:,3),15);
    
    % Calculate rate of change between coordinates:
    rate=[0; sqrt(diff(sx).^2+diff(sy).^2)./diff(xyz(:,1))];
    klein = (rate<mean(rate)/3); % Rate too small -> delete points
    xyz_temp=xyz;
    xyz(klein,:)=[];

    if length(xyz(:,1))<2 % if too many points were deleted -> recover original xyz (for interpolation)
        xyz=xyz_temp;
        klein_flag=1;
    else 
        klein_flag=0;
    end
    quality(klein) = [];
    q = quality;
    quality = NaN(size(trnum1));
    quality(xyz(:,1)) = q;
    quality=interp1(xyz(:,1),q,trnum1,'nearest',NaN);
   
    x=interp1(xyz(:,1),xyz(:,2),trnum1,'linear',NaN);
    y=interp1(xyz(:,1),xyz(:,3),trnum1,'linear',NaN);
    z=interp1(xyz(:,1),xyz(:,4),trnum1,'linear',NaN);

    % convert time to numerical value
    if matlab
        posTime = convertTo(time,'posixtime');
        if klein_flag==0
            t=interp1(xyz(:,1),posTime(~klein),trnum1,'linear',NaN);
        else
            t=interp1(xyz(:,1),posTime,trnum1,'linear',NaN);
        end
    else
        t = zeros(size(trnum1));
    end

else
    x=zeros(1,length(data(1,:))/rh_nchan);
    y=zeros(1,length(data(1,:))/rh_nchan);
    z=zeros(1,length(data(1,:))/rh_nchan);
    t=zeros(1,length(data(1,:))/rh_nchan);
    quality=zeros(1,length(data(1,:))/rh_nchan);
end
% Spuren ohne Koordinaten löschen
weg=find(isnan(x) | isnan(y));
x(weg)=[];
y(weg)=[];
z(weg)=[];
t(weg)=[];
quality(weg) = [];

if rh_nchan==2
    data(:,rh_ntraces/2+weg)=[];
end
data(:,weg)=[];
rh_ntraces=length(x);
if matlab
    time = datetime(t,'ConvertFrom','posixtime');
else
    time = t;
end

% File schlie�en:
fclose(fid);



%---------------------------------------------------------


%%% Header erstellen
% Profilname -> h.profname
h.profname=fname;

% Apparaturname -> h.appname
h.appname='SIR4000';

% Spuranzahl -> h.numtraces
if rh_nchan==2
    h.numtraces=length(data(1,:))/rh_nchan;
else
    h.numtraces=length(data(1,:));
end

% Anzahl Sample/Spur -> h.ns
h.ns=rh_nsamp;


if DF==1 % dualfrequ-Antenna
    h.range=2*dr*sqrt(8.854*1e-12*rh_epsr*4*pi*1e-7);
    h.range2=2*dr2*sqrt(8.854*1e-12*rh_epsr*4*pi*1e-7);
    h.dt=h.range/(rh_nsamp-1);
    h.dt2=h.range2/(rh_nsamp-1);
else
    % Aufzeichnungsdauer in ns -> h.range
    h.range=rh_range;

    % Sampleintervall in ns -> h.dt
    h.dt=rh_range/(rh_nsamp-1);
end

% Anzahl Marker -> h.anzmark
if ~all(mark==0)
    h.anzmark=length(mark);
else
    h.anzmark=0;
end

% Anzahl Kan�le
h.nchan=rh_nchan;



%%% Traceheader erstellen
% Spurnummer -> trh.tracenum
% Kanalnummer -> trh.channum
if rh_nchan==1
    trh.tracenum=1:h.numtraces;
    trh.channum=ones(size(trh.tracenum));
else
    trh.tracenum=1:h.numtraces;
    trh.channum(1:h.numtraces)=1;

    trh.tracenum(h.numtraces+1:h.numtraces*2)=[h.numtraces+1:h.numtraces*2]-h.numtraces;
    trh.channum(h.numtraces+1:h.numtraces*2)=ones(1,h.numtraces).*2;
end


% Koordinaten  setzen
% X-Koordinate -> trh.x
% Y-Koordinate -> trh.y
% Z-Koordinate -> trh.z
if ~isempty(xyz) && length(xyz(:,1))==1  % CMP
    if rh_nchan==1
        trh.x=zeros(size(trh.tracenum));
        trh.y=zeros(size(trh.tracenum));
        trh.z=zeros(size(trh.tracenum));
        if all(mark==0)
            trh.mark=zeros(size(trh.tracenum));
        else
            trh.mark=mark;
        end
        trh.time = NaT(size(trh.tracenum));
        trh.quality = NaN(size(trh.tracenum));
    else
        trh.x=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.y=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.z=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        if all(mark==0)
            trh.mark=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        else
            trh.mark=[mark mark];
        end
        trh.time = [NaT(size(trh.tracenum)) NaT(size(trh.tracenum))];
        trh.quality = [NaN(size(trh.tracenum)) NaN(size(trh.tracenum))];
    end
else
    if rh_nchan==2
        trh.x=[x x];
        trh.y=[y y];
        trh.z=[z z];
        trh.time = [time time];
        trh.quality = [quality quality];

    else
        trh.x=x;
        trh.y=y;
        trh.z=z;
        trh.time=time;
        trh.quality=quality;
    end

end


% Marker -> trh.mark
b=1;
for i=1:h.numtraces*rh_nchan
    if b<=length(mark) && mark(b)==i
        trh.mark(i)=1;
        b=b+1;
    else
        trh.mark(i)=0;
    end
end


if plotten==1
    t=0:h.dt:h.dt*(h.ns-1);
    for i=1:rh_nchan
        if i==1
            if DF==1
                t=t.*1e9;
            end
            figure
            imagesc([1:length(data(1,:))/rh_nchan],t,data(:,1:rh_ntraces))
            hold on
            colorbar
            if ~all(mark==0) && length(mark)<rh_ntraces/50
                for ii=1:length(mark)
                    plot([mark(ii) mark(ii)],[0 10],'k')
                end
            end
            colormap(flipud(gray))
            xlabel('Trace number')
            ylabel('t [ns]')
        else
            if DF==0
                figure
                imagesc([1:length(data(1,:))/rh_nchan],t,data(:,rh_ntraces+1:end))
                hold on
                colorbar
                if ~all(mark==0) && length(mark)<rh_ntraces/50
                    for ii=1:length(mark)
                        plot([mark(ii) mark(ii)],[0 10],'k')
                    end
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            else
                t2=(0:h.dt2:h.dt2*(h.ns-1)).*1e9;
                figure
                imagesc([1:length(data(1,:))/rh_nchan],t2,data(:,rh_ntraces+1:end))
                hold on
                colorbar
                if ~all(mark==0) && length(mark)<rh_ntraces/50
                    for ii=1:length(mark)
                        plot([mark(ii) mark(ii)],[0 10],'k')
                    end
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end
        end
    end
end
end

function dtime = readTime(binTime)
% convert the binary coded time data
% into a matlab datetime object
%
% Source:
% GSSI-SIR 3000 Manual S. 55 (Addendum A)
% MN72-433 Rev L 02.20.2015



% convert to bit-string
sT= dec2bin(binTime);

% split bit-string into parts
yr = bin2dec(sT(1:6)) + 1980;
mo = bin2dec(sT(7:10));
day = bin2dec(sT(11:15));
hr = bin2dec(sT(16:20));
min = bin2dec(sT(21:26));
sec = bin2dec(sT(27:31)) *2;

% return datetime object
dtime = datetime(yr, mo, day, hr, min, sec);
end

function [x,y,z] = rmOutliers(x,y,z,wind)
	% RMOUTLIERS: removing outliers from coordinates 
	% Outliers over a threshold are removed and
	% replaced by interpolated values.
	% Morten Harms, CAU Kiel 2022, morten.harms@ifg.uni-kiel.de

	% median filter data
	xfilt = medfilt1(x,wind);
	yfilt = medfilt1(y,wind);
	zfilt = medfilt1(z,wind);
	xfilt(1) = x(1);
	yfilt(1) = y(1);
	zfilt(1) = z(1);

	% detect outliers by comparing to filtered data
	residuals = sqrt((x-xfilt) .^ 2 + (y - yfilt) .^ 2);

	% find all outliers greater than the threshold
	threshold = mean(residuals)*3; 

	valid = residuals' < threshold;

	% find outliers in z axis
	zresiduals = abs(z-zfilt)';
	zthreshold = mean(zresiduals)*3;
	zValid = zresiduals < zthreshold;

	% interpolate the data at the outliers
	summed = cumsum(valid - diff([1,valid])/2);
	zsummed = cumsum(zValid - diff([1,zValid])/2);
	xnew = interp1(1:nnz(valid),x(valid),summed,"linear","extrap")';
	ynew = interp1(1:nnz(valid),y(valid),summed,"linear","extrap")';
	znew = interp1(1:nnz(zValid),z(zValid),zsummed,"linear","extrap")';

	x = xnew;
	y = ynew;
	z = znew;
end
