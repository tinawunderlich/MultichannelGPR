%------------------------ EKKO_Convert -------------------------------------
%
% Converts Spidar HD- & DT1- & GPS-files of Sensors&Software equipment to segy and mat (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, CAU Kiel, 2024, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------
% STILL UNDER CONSTRUCTION! PLEASE DO NOT USE!

clear all
close all
clc

name='NIC'; % name of data files

% number of channels
antenna_IDs=[1:4 6]; 
gps_channel=1; % which channel is synchronized with the GPS? (-> same filename-number)

dataplot=1; % plot radargram for controlling? 1=yes, 0=no

offsetGPS_X=[-0.5,-0.25,0,0.25,0.5]; % Offset between GPS and antenna midpoint crossline (for each antenna midpoint, same order as antenna_IDs) (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=[0,0,0,0,0]; % Offset between GPS and antenna midpoint in profile direction (for each antenna midpoint, same order as antenna_IDs) (if GPS behind antenna midpoint -> positive)

% Options for calculating inline coordinates for each trace:
coords_opt=3;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
                % =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)
                % =3: use distance in GPS-file

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
    if exist('spidar.temp') % read last opened folder from temp.temp
        fid=fopen('spidar.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with spidar-file(s)');
        else
            pfad=uigetdir([],'Choose folder with spidar-file(s)');
        end
        fileattrib('spidar.temp','-h');
        fid=fopen('spidar.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('spidar.temp','+h');
    else
        pfad=uigetdir([],'Choose folder with spidar-file(s)'); % path to radargram-folder

        fid=fopen('spidar.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('spidar.temp','+h');
    end
else
    if exist('.spidar.temp') % read last opened folder from temp.temp
        fid=fopen('.spidar.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with spidar-file(s)');
        else
            pfad=uigetdir([],'Choose folder with spidar-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with spidar-file(s)'); % path to radargram-folder
    end

    fid=fopen('.spidar.temp','wt');
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
    temp=extractAfter(list_HD(i).name,name);
    list_HD(i).channum=str2num(temp(1:2)); % antenna id/channel number
    temp1=extractBetween(temp,'Line','.HD');
    list_HD(i).pnumber=str2num(temp1{1}); % Profile number
end
HD=struct2table(list_HD);

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
    temp=extractAfter(list_DT(i).name,name);
    list_DT(i).channum=str2num(temp(1:2)); % antenna id/channel number
    temp1=extractBetween(temp,'Line','.DT1');
    list_DT(i).pnumber=str2num(temp1{1}); % Profile number
end
% sort list
DT=struct2table(list_DT);


% get GPS coordinates:
list_gps=dir(fullfile(pfad,'*.GPS')); % data format
% check for names starting with .
ii=1;
while ii<length(list_gps)
    if strcmp(list_gps(ii).name(1),'.')
        list_gps(ii,:)=[];
    else
        ii=ii+1;
    end
end
% get profile number
for i=1:length(list_gps)
    temp1=extractBetween(list_gps(i).name,'Line','.GPS');
    list_gps(i).pnumber=str2num(temp1{1}); % Profile number
end


radargrams=[];
global_coords=[];
x=[];
marker=[];
anz=1;
listrad.name=[];
listrad.chan=[];
for i=1:length(list_gps) % go through all profiles
    disp(['Profile ',int2str(list_gps(i).pnumber),' of ',int2str(length(list_gps))])
    
    for j=1:length(antenna_IDs) % for each channel

        % Read header file
        HDname=table2array(HD(HD.(8)==list_gps(i).pnumber & antenna_IDs(j)==HD.(7),1)); % get correct file name for header
        temp=readlines(fullfile(pfad,HDname{1}));
        for line = temp'
            if contains(line, 'NUMBER OF TRACES')
                h.ntr = str2num(extractAfter(line, '= '));
            elseif contains(line, 'TOTAL TIME WINDOW')
                h.tmax = str2num(extractAfter(line, '= '));
            elseif contains(line, 'STARTING POSITION')
                h.start_pos = str2num(extractAfter(line, '= '));
            elseif contains(line, 'FINAL POSITION')
                h.ending_pos = str2num(extractAfter(line, '= '));
            elseif contains(line, 'STEP SIZE USED')
                h.dx = str2num(extractAfter(line, '= '));
            end
        end

        % Read data file:
        spidarname=table2array(DT(DT.(8)==list_gps(i).pnumber & antenna_IDs(j)==DT.(7),1));
        fid=fopen(fullfile(pfad,spidarname{1}),'r');
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

        if dataplot==1
            figure
            imagesc(table2array(trh(:,2)),linspace(0,h.tmax,table2array(trh(1,3))),data)
            colormap(flipud(gray))
            xlabel('x [m]')
            ylabel('t [ns]')
            coldata=sort(data(~isnan(data)));
            cmin=coldata(round(length(coldata)/100*1));
            cmax=coldata(end-round(length(coldata)/100*1));
            set(gca,'CLim',[cmin cmax])
        end
        
        % correct GPS-antenna offset:
        if offsetGPS_X~=0 || offsetGPS_Y~=0
            anz2=round(0.5/mean(sqrt(diff(trh.x).^2+diff(trh.y).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
            if anz2/2==round(anz2/2)
                anz2=anz2+1; % make odd
            end
            for ii=1:length(trh.x)-anz2
                dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
            end
            anz1=anz2;
            % calculation for the last few traces
            anz2=anz2-1;
            for ii=length(trh.x)-anz1+1:length(trh.x)-1
                dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
                anz2=anz2-1;
            end
            % extrapolate for last trace
            trh.x(end)=interp1([1 2],[trh.x(end-2) trh.x(end-1)],3,'linear','extrap');
            trh.y(end)=interp1([1 2],[trh.y(end-2) trh.y(end-1)],3,'linear','extrap');
        end

    trh=struct2table(trh);

    if i==1
        radargrams=[{data}];
        global_coords=[{[trh.x trh.y trh.z]}];
        xtemp=trh.x;
        ytemp=trh.y;
        if coords_opt==1 % difference to beginning
            x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
        elseif coords_opt==2 % cumulative sum
            x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
        end


        listrad.name=[{ekkoname}];
        listrad.chan=[1];
    else
        radargrams=[radargrams; {data}];
        global_coords=[global_coords; {[trh.x trh.y trh.z]}];
        xtemp=trh.x;
        ytemp=trh.y;
        if coords_opt==1 % difference to beginning
            x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
        elseif coords_opt==2 % cumulative sum
            x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
        end


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
