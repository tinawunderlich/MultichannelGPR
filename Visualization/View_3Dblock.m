clear all
close all
clc


% Display Mala-GPR-data as 3D block, preprocessed by Mala3D
%
% Dr. Tina Wunderlich, CAU Kiel 2019, tina.wunderlich@ifg.uni-kiel.de
%
% requires binned data in rectangles 3D_Grid_R*
%
% requires folder Subfunctions and Plotting



% Number of rectangles with binned data
rectangles=1; % e.g. 1:3

% depth or time domain?
dsl=1;  % =1: depth, =0: time

% use processed data (if =1, then use data in "3D_Grid_R*/processed" folder)
proc=1;

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
            pfad=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pfad=uigetdir([],'Choose rSlicer folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pfad=uigetdir([],'Choose rSlicer folder');
        end
    else
        pfad=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'),fullfile(curFold,'Subfunctions'));


% read time vector (or depth vector for dsl=1
if proc==0
    load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
else
    load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'processed','t.mat'));
end


% read coordinates
for i=1:length(rectangles)
    if proc==0
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'x.mat'));
        x{i}=temp.x;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'y.mat'));
        y{i}=temp.y;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'z.mat'));
        z{i}=temp.z;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'mask.mat'));
        mask{i}=temp.mask;
    else
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'processed','x.mat'));
        x{i}=temp.x;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'processed','y.mat'));
        y{i}=temp.y;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'processed','z.mat'));
        z{i}=temp.z;
        temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(i))],'processed','mask.mat'));
        mask{i}=temp.mask;
    end
    
    if i==1
        dx=x{i}(1,2)-x{i}(1,1);   % dx of 3D-block
    end
    
    xmin(i)=min(x{i}(1,:));
    xmax(i)=max(x{i}(1,:));
    ymin(i)=min(y{i}(:,1));
    ymax(i)=max(y{i}(:,1));
end

% make grids for complete block
[xgrid,ygrid,tgrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)],t);

data=zeros(size(xgrid));
% set data in big block
for i=1:length(rectangles)
    if proc==0
        matFileObj=matfile(fullfile(pfad,['3D_Grid_R',int2str(i)],'data.mat'));
    else
        matFileObj=matfile(fullfile(pfad,['3D_Grid_R',int2str(i)],'processed','data.mat'));
    end
    
    indx=find(xmin(i)<=xgrid(1,:,1) & xmax(i)>=xgrid(1,:,1));
    indy=find(ymin(i)<=ygrid(:,1,1) & ymax(i)>=ygrid(:,1,1));
    data(indy,indx,:)=matFileObj.data(1:length(indy),1:length(indx),:);
end

% plot in GUI
if dsl==0
    plot_3Dblock(data,xgrid,ygrid,tgrid,pfad);
else
    plot_3Dblock_depth(data,xgrid,ygrid,tgrid,pfad);
end

waitfor(gcf);

% set original path
path(oldpath);