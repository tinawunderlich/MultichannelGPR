function [data,trh,h]=readdzt_UtilityScan(dztfile,plotten)



%%% [data,trh,h]=readdzt_UtilityScan(dztfile,plotflag)
% for reading dzt-file from UtilityScan with DF antenna
%
% dztfile: Name der dzt-Datei
% plotten: wenn =1 -> Bild wird geplottet, sonst =0
%
% Output:
% data: Datensamples in einer Matrix
% trh: Traceheader f�r jede Spur
% h: Header f�r die gesamte Datei
% falls zwei Kan�le: alle Daten zusammen in data, aber im Header gibt es
% Kanalnummer
%
% Dr. Tina Wunderlich, CAU, 2024
%%%



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
                        % ectract the coordinate information
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

    xyz = [ trnum', x(ia)', y(ia)', z(ia)' ];


    %xyz=[trnum(1:length(xyz(:,1)))' xyz];


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
    
    mark=zeros(size(xyz(:,1)))';
else
    % Markerdatei
    fid=fopen([pathstr,'/',fname,'.DZX']);
    temp=textscan(fid,'%s');
    DF=1; % script for DF antenna explicitly!
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
        if strfind(temp{1}{i},'<unitsPerScan>')
            dx=str2num(temp{1}{i}(15:22));
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
if rh_nchan==2
    fseek(fid,1026,'bof');
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
    rh_range2=test(5);   % Range in ns

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
    rh_depth2=test(3);   % Range in meters

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
end

%---------------------------------------------------------
% fseek(fid,rh_data+rh_nsamp*rh_bit*8,'bof');   % Springen auf Anfang der Daten
% if rh_bit==16
%     test=fread(fid,'uint16','l'); % Rest einlesen zum gucken wie viele Spuren
% elseif rh_bit==32
%     test=fread(fid,'int32','l'); % Rest einlesen zum gucken wie viele Spuren
% end

fseek(fid,rh_data*rh_nchan+rh_nsamp*rh_bit*8,'bof');
if rh_bit==16
    test=fread(fid,'uint16','l'); % Rest einlesen zum gucken wie viele Spuren
elseif rh_bit==32
    test=fread(fid,'int32','l'); % Rest einlesen zum gucken wie viele Spuren
end

rh_ntraces=ceil(length(test)/rh_nsamp);  % Anzahl Spuren (richtig!)

%---------------------------------------------------------

fseek(fid,rh_data*rh_nchan+(rh_nsamp-1)*rh_bit*8,'bof'); % auf Datenanfang zurück setzen

data=zeros(rh_nsamp,rh_ntraces*rh_nchan); % Datenmatrix initialisieren

% Daten auslesen
datatemp=fread(fid,[rh_nsamp,rh_ntraces*rh_nchan],'int32','l');
% ersten zwei Zeilen Null setzen
datatemp(1:2,:)=0;
% /1000 und Mittelwert abziehen
datatemp=datatemp./1000;
datatemp=datatemp-mean(datatemp(:));
% ersten zwei Zeilen Null setzen
datatemp(1:2,:)=0;

if rh_nchan==1
    data=datatemp;
else
    data=datatemp(:,1:2:end);
    data=[data datatemp(:,2:2:end)]; % h�nge Kanal 2 hinten dran
end


% falls am Anfang oder Ende Leerspuren sind -> l�schen
weg=find(data(20,:)==0);
data(:,weg)=[];


% Koordinaten interpolieren und Stillstandspuren l�schen
if ~all(xyz(:,2)==0)
    %xyz(xyz(:,1)==0,:)=[]; % delete traces in beginning with tracenumber 0
    trnum1=1:length(data(1,:))/rh_nchan;
    [a,b]=unique(xyz(:,1));

  

    x=interp1(xyz(b,1),xyz(b,2),trnum1,'linear',NaN);
    y=interp1(xyz(b,1),xyz(b,3),trnum1,'linear',NaN);
    z=interp1(xyz(b,1),xyz(b,4),trnum1,'linear',NaN);

      % convert time to numerical value
    if matlab
        posTime = convertTo(time,'posixtime');
        t=interp1(xyz(b,1),posTime(b),trnum1,'linear',NaN);
    else
        t = zeros(size(trnum1));
    end
   
else
    x=zeros(1,length(data(1,:))/rh_nchan);
    y=zeros(1,length(data(1,:))/rh_nchan);
    z=zeros(1,length(data(1,:))/rh_nchan);
    t=zeros(1,length(data(1,:))/rh_nchan);
end
% Spuren ohne Koordinaten löschen
weg=find(isnan(x));
x(weg)=[];
y(weg)=[];
z(weg)=[];
t(weg)=[];
if rh_nchan==2
    data(:,weg*2)=[];
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
DF=1; % Tablet wird nur für DF-Antenne benutzt!

%%% Header erstellen
% Profilname -> h.profname
h.profname=fname;

% Apparaturname -> h.appname
h.appname='UtilityScan';

% Spuranzahl -> h.numtraces
if rh_nchan==2
    h.numtraces=length(data(1,:))/rh_nchan;
else
    h.numtraces=length(data(1,:));
end

% Anzahl Sample/Spur -> h.ns
h.ns=rh_nsamp;


if DF==1 % dualfrequ-Antenna
    h.range=rh_range;
    h.range2=rh_range2;
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
if length(xyz(:,1))==1  && all([all(x==0) all(y==0) all(z==0)])% CMP
    if rh_nchan==1
        trh.x=zeros(size(trh.tracenum));
        trh.y=zeros(size(trh.tracenum));
        trh.z=zeros(size(trh.tracenum));
        trh.mark=zeros(size(trh.tracenum));
        trh.time = NaT(size(trh.tracenum));
    else
        trh.x=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.y=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.z=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.mark=[zeros(size(trh.tracenum)) zeros(size(trh.tracenum))];
        trh.time = [NaT(size(trh.tracenum)) NaT(size(trh.tracenum))];
    end
else
    if all([all(x==0) all(y==0) all(z==0)]) && exist('dx','var')
        x=0:dx:(rh_ntraces-1)*dx;    %create coordinate vector from trace spacing
    end

    if rh_nchan==2
        trh.x=[x x];
        trh.y=[y y];
        trh.z=[z z];
        trh.time = [time time];

    else
        trh.x=x;
        trh.y=y;
        trh.z=z;
        trh.time=time;
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
            figure
            imagesc([1:length(data(1,:))/rh_nchan],t,data(:,1:rh_ntraces))
            hold on
            colorbar
            if ~all(mark==0)
                for ii=1:length(mark)
                    plot([mark(ii) mark(ii)],[0 10],'k')
                end
            end
        else
            if DF==0
                figure
                imagesc([1:length(data(1,:))/rh_nchan],t,data(:,rh_ntraces+1:end))
                hold on
                colorbar
                if ~all(mark==0)
                    for ii=1:length(mark)
                        plot([mark(ii) mark(ii)],[0 10],'k')
                    end
                end
            else
                t2=0:h.dt2:h.dt2*(h.ns-1);
                figure
                imagesc([1:length(data(1,:))/rh_nchan],t2,data(:,rh_ntraces+1:end))
                hold on
                colorbar
                if ~all(mark==0)
                    for ii=1:length(mark)
                        plot([mark(ii) mark(ii)],[0 10],'k')
                    end
                end
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
    % Morten Harms, CAU, Februar 2022
    

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