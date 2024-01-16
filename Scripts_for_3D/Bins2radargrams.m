clear all
close all
clc


% Use 3D-binned data to extract radargrams along bin intervals in both
% directions (these radargrams can be visualized e.g. using Tsl_Radargrams.m)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Requires binned 3D data in folders 3D_Grid_R*, if asked give rSlicer
% folder (folder with 3D_Grid_R*-folders inside)


rectangles=1:12; % give number of rectangles

% If you want to get only some profiles in North-South or East-West
% direction, specify Easting and/or Northing. If you want to have all
% radargrams along all bins, set both to []. Only possible if area has not
% been rotated!
Easting=[465570, 465600];
Northing=[6060229];
coordtrans=1; % if bin coordinates are local and Easting/Northing are given in global coords, set coordtrans=1

proc=0;  % if you want to use processed 3D_blocks in 3D_Grid_R*/processed -> =1, else =0

% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name 
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose rSlicer-folder');
        else
            pfad_rad=uigetdir([],'Choose rSlicer-folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad_rad=uigetdir([],'Choose rSlicer-folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose rSlicer-folder');
        else
            pfad_rad=uigetdir([],'Choose rSlicer-folder');
        end
    else
        pfad_rad=uigetdir([],'Choose rSlicer-folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));


% read time vector
if proc==0
    ttmp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
else
    ttmp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'processed','t.mat'));
end
t=ttmp.t';

% read coordinates
for i=1:length(rectangles)
    if proc==0
        temp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(i))],'x.mat'));
    else
        temp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(i))],'processed','x.mat'));
    end
    xx{i}=temp.x;
    if proc==0
        temp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(i))],'y.mat'));
    else
        temp=load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(i))],'processed','y.mat'));
    end
    yy{i}=temp.y;

    
    if i==1
        dx=xx{i}(1,2)-xx{i}(1,1);   % dx of 3D-block
    end
    
    xmin(i)=min(xx{i}(1,:));
    xmax(i)=max(xx{i}(1,:));
    ymin(i)=min(yy{i}(:,1));
    ymax(i)=max(yy{i}(:,1));
end

% make large grid for all rectangles
[xgrid,ygrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)]);

if exist(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'coordtrans.mat'),'file')
    load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'coordtrans.mat'));
    flag_co=1;
else
    flag_co=0;
end

% x and y coordinates for profiles
if isempty(Easting)
    xlines=unique(xgrid); % all bins
else
    if flag_co==0
        xlines=Easting; % just some profiles in North-South direction
    else
        for i=1:length(Easting)
            temp=helmert([Easting(i) 0],coordtrans(:,3:4),coordtrans(:,1:2));
            xlines(i)=temp(1);
        end
    end
end
if isempty(Northing)
    ylines=unique(ygrid); % all bins
else
    if flag_co==0
        ylines=Northing; % just some profiles in West-East direction
    else
        for i=1:length(Northing)
            temp=helmert([0 Northing(i)],coordtrans(:,3:4),coordtrans(:,1:2));
            ylines(i)=temp(2);
        end
    end
end

% open rectangle data:
for j=1:length(rectangles)
    if proc==1
        matFileObj{j}=matfile(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'processed',['data.mat']));
    else
        matFileObj{j}=matfile(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'data.mat'));
    end
end


% profiles for each x-coordinate (along y, in NS direction)
disp('Reading data... Please wait.');
h=waitbar(0,'Creating profiles along y.');
for i=1:length(xlines)
    radargrams{i}=NaN(length(t),length(xgrid(:,1)));
    global_coords{i}=[xlines(i) ygrid(1,1); xlines(i) ygrid(end,1)];
    x{i}=ygrid(:,1)-ygrid(1,1);
   
    % set data in radargram:
    for j=1:length(rectangles)
        if xmin(j)<=xlines(i) && xmax(j)>=xlines(i)
            ind=find(round(xlines(i)*100)==round(xx{j}(1,:).*100));
            tmp=matFileObj{j}.data(:,ind,:);
            radargrams{i}(:,round(ygrid(:,1).*100)>=round(ymin(j).*100) & round(ygrid(:,1).*100)<=round(ymax(j).*100))=permute(tmp,[3 1 2]);
        end
        waitbar(((i-1)*length(rectangles)+j)/(length(xlines)*length(rectangles)),h);
    end
end
close(h);
    

% profiles for each y-coordinate (along x, in WE direction)
disp('Reading data... Please wait.');
h=waitbar(0,'Creating profiles along x.');
n=length(radargrams);
for i=1+n:n+length(ylines)
    radargrams{i}=NaN(length(t),length(xgrid(1,:)));
    global_coords{i}=[xgrid(1,1) ylines(i-n); xgrid(1,end) ylines(i-n)];
    x{i}=xgrid(1,:)-xgrid(1,1);
   
    % set data in radargram:
    for j=1:length(rectangles)
        if ymin(j)<=ylines(i-n) && ymax(j)>=ylines(i-n)
            ind=find(round(ylines(i-n)*100)==round(yy{j}(:,1).*100));
            tmp=matFileObj{j}.data(ind,:,:);
            indstart=find(round(xgrid(1,:).*100)>=round(xmin(j).*100),1,'first'); % startindex
            radargrams{i}(:,indstart:indstart+length(tmp(1,:,1))-1)=permute(tmp,[3 1 2]);
            %radargrams{i}(:,round(xgrid(1,:).*100)>=round(xmin(j).*100) & round(xgrid(1,:).*100)<=round(xmax(j).*100))=permute(tmp,[3 1 2]);
        end
        waitbar(((i-1-n)*length(rectangles)+j)/(length(xlines)*length(rectangles)),h);
    end
end
close(h);
    



% Save radargrams
if proc==0
    if ~exist(fullfile(pfad_rad,'Radargrams_bins'))
        mkdir(fullfile(pfad_rad,'Radargrams_bins'));
    end
    save(fullfile(pfad_rad,'Radargrams_bins','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins','global_coords.mat'),'global_coords','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins','t.mat'),'t','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins','x.mat'),'x','-v7.3');
else
    if ~exist(fullfile(pfad_rad,'Radargrams_bins_proc'))
        mkdir(fullfile(pfad_rad,'Radargrams_bins_proc'));
    end
    save(fullfile(pfad_rad,'Radargrams_bins_proc','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins_proc','global_coords.mat'),'global_coords','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins_proc','t.mat'),'t','-v7.3');
    save(fullfile(pfad_rad,'Radargrams_bins_proc','x.mat'),'x','-v7.3');
end

% set original path
path(oldpath);