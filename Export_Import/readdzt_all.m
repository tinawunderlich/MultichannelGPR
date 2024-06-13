function [data,trh,h]=readdzt_all(dztfile,plotten,app) 

%%% [data,trh,h]=readdzt(dztfile,plotten,app) 
%
% dztfile: Name der dzt-Datei
% plotten: wenn =1 -> Bild wird geplottet, sonst =0
% app: Apparatur ('SIR3000' or 'SIR20')
%
% Output: 
% data: Datensamples in einer Matrix
% trh: Traceheader f?r jede Spur
% h: Header f?r die gesamte Datei
%
% Dr. Tina Wunderlich, CAU, November 2017
%%%

% File ?ffnen:
fid=fopen(dztfile,'r','l');

test=fread(fid,1,'int16','l');

% name zerlegen
[pathstr,fname,ext] = fileparts(dztfile);

%%% Radar Header auslesen:

% Apparatur-Tag
SIR10=0;
SIR20=0;
SIR30=0;
SIR3000=0;
SIR4000=0;
if strcmp(app,'SIR3000')
    SIR3000=1;
elseif strcmp(app,'SIR20')
    SIR20=1;
end
% % dann die entsprechende Apparatur raussuchen
% if test(1)==2047
%     SIR30=1;
%     SIR20=1;
%     SIR4000=1;
% elseif test(1)==255
%     SIR10=1;
% elseif test(1)==1792
%     SIR3000=1;
% end
rh_tag=test(1);

%---------------------------------------------------------
test=fread(fid,3,'uint16','l');   % weiter lesen...

rh_data=test(1);    % Offset to data from beginning of file
rh_nsamp=test(2);   % Samples per Scan
rh_bit=test(3);     % Bits per data word (16 = uint16 f?r sp?teres Dateneinlesen)

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

rh_antname=test;    % Transmitrate (?)

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
test=fread(fid,1024-128-2*(9),'uchar=>uchar','l');  % weiter lesen...

variable=test;  % rh_rgain, rh_text und rh_proc

%---------------------------------------------------------
if rh_nchan==1   % f?r einen Kanal
    fseek(fid,rh_data+rh_nsamp*64,'bof');   % Springen auf Anfang der Daten
    if rh_bit==16
        test=fread(fid,'uint16','l'); % Rest einlesen zum gucken wie viele Spuren
    elseif rh_bit==32
        test=fread(fid,'uint32','l'); % Rest einlesen zum gucken wie viele Spuren
    end

    rh_ntraces=length(test)/rh_nsamp;  % Anzahl Spuren
    
    %---------------------------------------------------------
    fseek(fid,rh_data+rh_nsamp*64,'bof');   % Zur?cksetzen auf Anfang der Daten

    % Daten auslesen
    if SIR30==1 && rh_bit==32
        fseek(fid,rh_data+rh_nsamp*64-64*2,'bof'); % keine Ahnung warum... ausprobiert!
        data=fread(fid,[rh_nsamp,rh_ntraces],'int32','l');
    else
        data=fread(fid,[rh_nsamp,rh_ntraces],'uint16','l')-32768;
    end

    if SIR30==1 %|| SIR3000==1
        data(1:2,:)=zeros(size(data(1:2,:)));   % die ersten zwei Zeilen 0 setzen
    end
    
elseif rh_nchan==2   % f?r 2 Kan?le
    fseek(fid,rh_data*2,'bof');   % Springen auf Anfang der Daten
    if rh_bit==16
        test=fread(fid,'uint16','l'); % Rest einlesen zum gucken wie viele Spuren
    elseif rh_bit==32
        test=fread(fid,'uint32','l'); % Rest einlesen zum gucken wie viele Spuren
    end

    rh_ntraces=length(test)/rh_nsamp;  % Anzahl Spuren (f?r beide Kan?le zusammen)
    
    %---------------------------------------------------------
    fseek(fid,rh_data*2,'bof');   % Zur?cksetzen auf Anfang der Daten
    
    % Daten auslesen
    temp=fread(fid,[rh_nsamp,rh_ntraces],'uint16','l')-32768;   % immer abwechselnd 1 Spur Kanal 1, dann 1 Spur Kanal 2
    % auf zwei Kan?le aufteilen
    data=temp(:,1:2:end);
    data=[data temp(:,2:2:end)]; % h�nge Kanal 2 hinten dran

end

weg=find(data(150,:)==0);
data(:,weg)=[];
%---------------------------------------------------------


% Zeitvektor berechnen
dt=rh_range/rh_nsamp;
t=[0:dt:rh_range-dt]';

rh_ntraces=length(data(1,:))/rh_nchan; % number of traces for one channel

% Marker suchen (nur in Kanal 1)
if SIR20==1

    % finde alle stellen mit 1
    test=data(1,1:rh_ntraces);
    eins=find(test==1);

    % Spr?nge suchen
    m=1;
    anfang=1;
    for i=1:length(eins)-1
        if eins(i+1)-eins(i)>2
            ende(m)=i;
            anfang(m+1)=i+1;
            m=m+1;
        end
    end
    ende(m)=length(eins);


    % Bereiche durchgehen
    if ~isempty(eins)
        for i=1:length(ende)
            if mean(test(eins(anfang(i)):eins(ende(i))))==1   % nur Einsen
                mark(i)=round((eins(ende(i))-eins(anfang(i)))/2)+eins(anfang(i));
            else
                temp=find(test(eins(anfang(i)):eins(ende(i)))~=1);
                mark(i)=eins(anfang(i))+temp(1)-1;
            end
        end
    end
elseif SIR3000==1
    val=unique(data(2,:));
    test=[sum(data(2,:)==val(1)) sum(data(2,:)==val(2))]; % get number of occurences
    mark=find(data(2,:)==val(min(test)==test));

%     fseek(fid,rh_data*2,'bof');   % Zur?cksetzen auf Anfang der Daten
%     for i=1:rh_ntraces
%         test=fread(fid,4,'char*1=>uint8','l');
%         marker(i)=test(4);
%         fseek(fid,rh_data+rh_nsamp*i*2,'bof');    % zum n?chsten Spuranfang setzen
%     end
% 
%     mark=find(marker<mean(marker));
end
    
if ~ exist('mark','var')
    mark=zeros(1,length(data(1,:)));
end

if plotten==1
    for i=1:rh_nchan
        if i==1
            figure
            imagesc([1:length(data(1,:))/rh_nchan],t,data(:,1:rh_ntraces))
            hold on
            colorbar
            if ~all(mark==0)
                hold on
                for i=1:length(mark)
                    plot([mark(i) mark(i)],[0 10],'k')
                end
            end
        else
            figure
            imagesc([1:length(data(1,:))/rh_nchan],t,data(:,rh_ntraces+1:end))
            hold on
            colorbar
            if ~all(mark==0)
                hold on
                for i=1:length(mark)
                    plot([mark(i) mark(i)],[0 10],'k')
                end
            end
        end
    end
end



%---------------------------------------------------------
% File schlie?en:
fclose(fid);



%%% Header erstellen
% Profilname -> h.profname
h.profname=fname;

% Apparaturname -> h.appname
if SIR20==1
    h.appname='SIR20';
elseif SIR30==1
    h.appname='SIR30';
elseif SIR3000==1
    h.appname='SIR3000';
end

% Spuranzahl -> h.numtraces (pro channel)
h.numtraces=rh_ntraces;

% Anzahl Sample/Spur -> h.ns
h.ns=rh_nsamp;

% Aufzeichnungsdauer in ns -> h.range
h.range=rh_range;

% Sampleintervall in ns -> h.dt
h.dt=rh_range/rh_nsamp;

% Anzahl Marker -> h.anzmark (pro channel)
h.anzmark=length(mark);

if rh_nchan>1
    marker=[];
    for i=1:rh_nchan
        marker=[marker mark];
    end
    mark=marker;    % overwrite mark with complete markers for both channels
end

% Anzahl Kan�le
h.nchan=rh_nchan;


%%% Traceheader erstellen
% Spurnummer -> trh.tracenum
if rh_nchan==1
    trh.tracenum=1:h.numtraces;
    trh.channum=ones(size(trh.tracenum));
else
    trh.tracenum=1:h.numtraces;
    trh.channum(1:h.numtraces)=1;
    
    trh.tracenum(h.numtraces+1:h.numtraces*2)=[h.numtraces+1:h.numtraces*2]-h.numtraces;
    trh.channum(h.numtraces+1:h.numtraces*2)=ones(1,h.numtraces).*2;
end

% Koordinaten erstmal 0 setzen
% X-Koordinate -> trh.x
% Y-Koordinate -> trh.y
trh.x(1:h.numtraces*rh_nchan)=0;
trh.y(1:h.numtraces*rh_nchan)=0;
trh.z(1:h.numtraces*rh_nchan)=0;

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