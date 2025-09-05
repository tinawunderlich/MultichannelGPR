%------------------------ MalaSingle_Convert -------------------------------------
%
% Converts single antenna Mala data to segy and mat (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, August 2025, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


clear all
close all
clc

dataplot=1; % plot radargram for controlling? 1=yes, 0=no

use_GPS=1; % if =1: GPS was used, if =0: local coordinates
convert2utm=1; % convert WGS84 Lat/Long to UTM
zone=36; % if convert2utm==1 -> give UTM-zone
offsetGPS_X=0; % Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=0; % Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)

% Options for calculating inline coordinates for each trace:
coords_opt=1;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
                % =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)

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
    if exist('mala.temp') % read last opened folder from temp.temp
        fid=fopen('mala.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with Mala-file(s)');
        else
            pfad=uigetdir([],'Choose folder with Mala-file(s)');
        end
        fid=fopen('mala.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with Mala-file(s)'); % path to radargram-folder

        fid=fopen('mala.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.mala.temp') % read last opened folder from temp.temp
        fid=fopen('.mala.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with Mala-file(s)');
        else
            pfad=uigetdir([],'Choose folder with Mala-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with Mala-file(s)'); % path to radargram-folder
    end

    fid=fopen('.mala.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% get list of rad-files in this folder
if use_GPS==0
    list=dir(fullfile(pfad,'*.RAD'));
else
    list=dir(fullfile(pfad,'*.rad'));
end

% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end
% name
name=extractBefore(list(1).name,'_');
% get profile number
temp11=extractBefore(list(1).name,'.');
temp12=strsplit(temp11,'_'); % parts of first file name
temp21=extractBefore(list(2).name,'.');
temp22=strsplit(temp21,'_'); % parts of second file name
for i=1:length(temp22)
    if ~strcmp(temp12{i},temp22{i})
        % different -> file number
        w=i;
        break;
    end
end
for i=1:length(list)
    % profile number
    temp=extractBefore(list(i).name,'.');
    temp2=strsplit(temp,'_'); % parts of file name
    list(i).number=str2num(temp2{w});
    number(i)=list(i).number;
    % name without ending
    list(i).name=temp;
end


radargrams=[];
global_coords=[];
x=[];
marker=[];
anz=1;
listrad.name=[];
for i=1:length(list)
    disp(['File ',int2str(i),' of ',int2str(length(list))])
    
    [traces,dt1,ns,xx,yy,zz]=readmala_single(pfad,list(i).name,use_GPS);

    if convert2utm==1
        [xx,yy]=wgs2utm(yy,xx,zone,'N');
        clear lat lon;
    end

    if dt1~=0
        dt=dt1; % correct range is only written in first file -> take dt from this file
    end

    if dataplot==1
        figure
        dist=sqrt((xx-xx(1)).^2+(yy-yy(1)).^2);
        imagesc(dist,0:dt:(ns-1)*dt,traces)
        colormap(flipud(gray))
        xlabel('x [m]')
        ylabel('t [ns]')
    end

    % correct GPS-antenna offset:
    if use_GPS==1 && (offsetGPS_X~=0 || offsetGPS_Y~=0)
        anz2=round(0.5/mean(sqrt(diff(xx).^2+diff(yy).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
        if anz2/2==round(anz2/2)
            anz2=anz2+1; % make odd
        end
        for ii=1:length(xx)-anz2
            dist=sqrt((xx(ii)-xx(ii+anz2))^2+(yy(ii)-yy(ii+anz2))^2);
            temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[xx(ii) yy(ii); xx(ii+anz2) yy(ii+anz2)]);
            xx(ii)=temp(1);
            yy(ii)=temp(2);
        end
        anz1=anz2;
        % calculation for the last few traces
        anz2=anz2-1;
        for ii=length(xx)-anz1+1:length(xx)-1
            dist=sqrt((xx(ii)-xx(ii+anz2))^2+(yy(ii)-yy(ii+anz2))^2);
            temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[xx(ii) yy(ii); xx(ii+anz2) yy(ii+anz2)]);
            xx(ii)=temp(1);
            yy(ii)=temp(2);
            anz2=anz2-1;
        end
        % extrapolate for last trace
        xx(end)=interp1([1 2],[xx(end-2) xx(end-1)],3,'linear','extrap');
        yy(end)=interp1([1 2],[yy(end-2) yy(end-1)],3,'linear','extrap');
    end

    if i==1
        radargrams=[{traces}];
        global_coords=[{[xx yy zz]}];
        xtemp=xx;
        ytemp=yy;
        if coords_opt==1 % difference to beginning
            x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
        elseif coords_opt==2 % cumulative sum
            x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
        end

        listrad.name=[{list(i).name}];
    else
        radargrams=[radargrams; {traces}];
        global_coords=[global_coords; {[xx yy zz]}];
        xtemp=xx;
        ytemp=yy;
        if coords_opt==1 % difference to beginning
            x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
        elseif coords_opt==2 % cumulative sum
            x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
        end

        listrad.name=[listrad.name; {list(i).name}];
    end

    t=0:dt:dt*(ns-1);
    anz=anz+1;

end

%------------------------------ WRITE DATA --------------------------------

if export2mat==1

    % delete empty cells
    radargrams(~cellfun('isempty',radargrams));

    disp('Start saving of mat-files...')

    % save all profiles in one variable
    % for one channel
    % make new folder
    if ~exist(fullfile(pfad,'mat'),'dir')
        mkdir(fullfile(pfad,'mat'));
    end
    save(fullfile(pfad,'mat','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pfad,'mat','t.mat'),'t','-v7.3');
    save(fullfile(pfad,'mat','x.mat'),'x','-v7.3');
    save(fullfile(pfad,'mat','global_coords.mat'),'global_coords','-v7.3');

    fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
    fprintf(fid,'Nr.\tName\n');
    for i=1:length(listrad.name)
        fprintf(fid,'%d\t',i);
        fprintf(fid,'%s\n',listrad.name{i});
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

        export2sgy2D(radargrams{i},dt,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'.sgy']),global_coords{i}(:,3),constoff);

    end
    
    disp('   Finished!')
end

% set original path
path(oldpath);
