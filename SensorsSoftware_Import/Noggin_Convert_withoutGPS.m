%------------------------ EKKO_Convert -------------------------------------
%
% Converts Noggin HD- & DT1-files of Sensors&Software EKKO equipments to segy and mat (compatible with Multichannel-GPR)
% measured with survey wheel and without GPS
%
% Dr. Tina Wunderlich, November 2023, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


clear all
close all
clc

name='liney'; % name of data files

line_spacing=1; % spacing between lines [m]

dataplot=0; % plot radargram for controlling? 1=yes, 0=no

% if you find out that the survey wheel was not calibrated, use following
% information to re-calibrate it (i.e. calculate new trace spacing for all
% profiles)
calibration=1; % do you want to recalibrate? 1=yes, 0=no
profile1length=52.9; % length of first profile [m]

% Export to other formats
export2mat=1; % export to Multichannel-GPR format for radargrams (mat-files)
export2segy=0; % export all radargrams as segy-files
constoff=0; % for sgy: if=1: a constant coordinate offset will be subtracted and coordinates will be in mm accuracy in segy file (offsets will be saved in Inline3D (x) and Crossline3D (y))


%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
% 
% LOAD DATA

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));


% get file names
if ispc
    if exist('ekko.temp') % read last opened folder from temp.temp
        fid=fopen('ekko.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with EKKO-file(s)');
        else
            pfad=uigetdir([],'Choose folder with EKKO-file(s)');
        end
        fid=fopen('ekko.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with EKKO-file(s)'); % path to radargram-folder

        fid=fopen('ekko.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.ekko.temp') % read last opened folder from temp.temp
        fid=fopen('.ekko.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with EKKO-file(s)');
        else
            pfad=uigetdir([],'Choose folder with EKKO-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with EKKO-file(s)'); % path to radargram-folder
    end

    fid=fopen('.ekko.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% get list of HD-files in this folder
list_HD=dir(fullfile(pfad,'*.HD')); % header format
% check for names starting with .
ii=1;
while ii<length(list_HD)
    if strcmp(list_HD(ii).name(1),'.')
        list_HD(ii,:)=[];
    else
        ii=ii+1;
    end
end
% get profile number
for i=1:length(list_HD)
    a=isstrprop(list_HD(i).name,'digit'); % find numbers in string
    list_HD(i).number=str2num(list_HD(i).name(a));
end

% get list of DT1-files in this folder
list_DT=dir(fullfile(pfad,'*.DT1')); % data format
% check for names starting with .
ii=1;
while ii<length(list_DT)
    if strcmp(list_DT(ii).name(1),'.')
        list_DT(ii,:)=[];
    else
        ii=ii+1;
    end
end
% get profile number
for i=1:length(list_DT)
    temp2=extractBefore(list_DT(i).name,'.DT1');
    a=isstrprop(temp2,'digit'); % find numbers in string
    list_DT(i).number=str2num(list_DT(i).name(a));
end
% sort list
DT=struct2table(list_DT);
DT=sortrows(DT,'number');
list_DT=table2struct(DT); % sorted list after profile number


radargrams=[];
global_coords=[];
x=[];
marker=[];
anz=1;
listrad.name=[];
listrad.chan=[];
for i=1:length(list_DT)

    disp(['File ',int2str(i),' of ',int2str(length(list_DT))])
    
    % Read header file
    for j=1:length(list_HD)
        if list_DT(i).number==list_HD(j).number
            HDname=list_HD(j).name; % get correct file name for header
        end
    end
    temp=readlines(fullfile(pfad,HDname));
    % number of traces
    a=cellfun(@(x) contains(x,'NUMBER OF TRACES'),temp);
    h.ntr=str2num(extractAfter(temp{a},'= '));
    % time window
    a=cellfun(@(x) contains(x,'TOTAL TIME WINDOW'),temp);
    h.tmax=str2num(extractAfter(temp{a},'= '));
    % starting pos
    a=cellfun(@(x) contains(x,'STARTING POSITION'),temp);
    h.start_pos=str2num(extractAfter(temp{a},'= '));
    % ending pos
    a=cellfun(@(x) contains(x,'FINAL POSITION'),temp);
    h.ending_pos=str2num(extractAfter(temp{a},'= '));
    % step size
    a=cellfun(@(x) contains(x,'STEP SIZE USED'),temp);
    h.dx=str2num(extractAfter(temp{a},'= '));

    % re-calibrate?
    if calibration==1 && i==1
        h.dx=profile1length/h.ntr; % determine dx for all profiles
        dx_all=h.dx;
        % adjust ending pos:
        h.ending_pos=h.dx*(h.ntr-1);
    end
    if calibration==1 && i>1 % set this dx for all other profiles
        if h.dx<0
            % reverse line
            h.dx=-dx_all; 
            % also adjust starting_pos:
            h.start_pos=dx_all*(h.ntr-1);
        else
            h.dx=dx_all;
            % also adjust ending pos:
            h.ending_pos=h.dx*(h.ntr-1);
        end
    end

    % -----------------------------------------------------
    % Read data file:
    ekkoname=list_DT(i).name;

    fid=fopen(fullfile(pfad,ekkoname),'r');
    for j=1:h.ntr % read every trace
        % header
        trh(j).tracenum=fread(fid,1,'float');
        trh(j).pos=fread(fid,1,'float');
        trh(j).ns=fread(fid,1,'float');
        tmp=fread(fid,3,'float');
        trh(j).topo=tmp(1);
        trh(j).twindow=fread(fid,1,'float');
        trh(j).nstacks=fread(fid,1,'float');
        trh(j).x=fread(fid,1,'double');
        trh(j).y=fread(fid,1,'double');
        trh(j).z=fread(fid,1,'double');
        tmp=fread(fid,7,'float');
        trh(j).xr=tmp(1); % Receiver
        trh(j).yr=tmp(2);
        trh(j).zr=tmp(3);
        trh(j).xt=tmp(4); % transmitter
        trh(j).yt=tmp(5);
        trh(j).zt=tmp(6);
        tmp=fread(fid,11,'float');
        trh(j).timezero=tmp(1);
        trh(j).zeroflag=tmp(2); % 0=data ok, 1=no data
        % data
        if j==1
            data=zeros(trh(j).ns,h.ntr);
        end
        data(:,j)=fread(fid,trh(j).ns,'int16');
    end
    trh=struct2table(trh);
    fclose(fid);

    
    % current line coordinate (x or y, will be determined automatically)
    line_coord=line_spacing*(i-1);
    % calculate trace coords:
    if all(trh.xr==0)
        trh.y=linspace(h.start_pos,h.ending_pos,h.ntr)';
        trh.x=zeros(size(trh.x))+line_coord;
    else
        trh.x=linspace(h.start_pos,h.ending_pos,h.ntr)';
        trh.y=zeros(size(trh.y))+line_coord;
    end

    if dataplot==1
        figure
        imagesc(table2array(trh(:,2)),linspace(0,h.tmax,table2array(trh(1,3))),data)
        colormap(flipud(gray))
        xlabel('x [m]')
        ylabel('t [ns]')
    end

    if i==1
        radargrams=[{data}];
        global_coords=[{[trh.x trh.y trh.z]}];
        xtemp=trh.x;
        ytemp=trh.y;
        x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2))]}];

        listrad.name=[{ekkoname}];
        listrad.chan=[1];
    else
        radargrams=[radargrams; {data}];
        global_coords=[global_coords; {[trh.x trh.y trh.z]}];
        xtemp=trh.x;
        ytemp=trh.y;
        x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2))]}];


        listrad.name=[listrad.name; {ekkoname}];
        listrad.chan=[listrad.chan; 1];
    end

    t=linspace(0,h.tmax,table2array(trh(1,3)));
    if max(t)<1 % for DF-Antenna is dt in s
        t=t.*1e9; % convert to ns
    end
    anz=anz+1;

    % collect headers in cell
    headers{i} = h;

    clear trh;

end

%------------------------------ WRITE DATA --------------------------------

if export2mat==1

    % delete empty cells
    radargrams(~cellfun('isempty',radargrams));

    disp('Start saving of mat-files...')

    % save all profiles in one variable
    % for one channel/not DF antenna
    % make new folder
    if ~exist(fullfile(pfad,'mat'),'dir')
        mkdir(fullfile(pfad,'mat'));
    end
    save(fullfile(pfad,'mat','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pfad,'mat','t.mat'),'t','-v7.3');
    save(fullfile(pfad,'mat','x.mat'),'x','-v7.3');
    save(fullfile(pfad,'mat','global_coords.mat'),'global_coords','-v7.3');
    save(fullfile(pfad,'mat','h.mat'),'headers','-v7.3');

    fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
    fprintf(fid,'Nr.\tName\tChannel\n');
    for i=1:length(listrad.name)
        fprintf(fid,'%d\t',i);
        fprintf(fid,'%s\t',listrad.name{i});
        fprintf(fid,'%d\n',listrad.chan(i));
    end
    fclose(fid);

    disp('   Finished!')
end


if export2segy==1
    disp('Start saving as sgy...')
    % make new folder
    if ~exist(fullfile(pfad,'SEGY'),'dir')
        mkdir(fullfile(pfad,'SEGY'));
    end
    
    for i=1:length(radargrams)
        if length(global_coords{i}(1,:))<3
            global_coords{i}(:,3)=zeros(size(global_coords{i}(:,1))); % if no topography present, set to zero
        end

        export2sgy2D(radargrams{i},h.dt,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);

    end
    
    disp('   Finished!')
end

% set original path
path(oldpath);
