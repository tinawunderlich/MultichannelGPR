clear all
close all
clc


% Use 3D-binned data and extract a smaller 3D-block from this. This block
% will be saved in folder Small_3Dblock/3D_Grid_R1 (for compatibility with the other scripts also in a 3D_Grid_R1 folder).
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Requires binned 3D data in folders 3D_Grid_R*, if asked give rSlicer
% folder (folder with 3D_Grid_R*-folders inside)


rectangles=1:15; % give number of rectangles

% give coordinate of lower left corner of small area, and width/height
corner=[585254 5971678]; % Easting/Northing
width=5; % width of small area (Easting/x) in m
height=15; % height of small area (Northing/y) in m
local_global=1; % =1: width/height in local(rotated) coordinate system; =2: width/height in global coordinate system

proc=0;  % if you want to use processed 3D_blocks in 3D_Grid_R*/processed -> =1, else =0

% Interpolation
interp=1; % =1: do linear interpolation for every slice, =0: no interpolation
radius=3; % radius for interpolation in bins

% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name 
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose rSlicer-folder');
        else
            pfad_rad=uigetdir([],'Choose rSlicer-folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
    else
        pfad_rad=uigetdir([],'Choose rSlicer-folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
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

if exist(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'coordtrans.mat'),'file')
    load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(1))],'coordtrans.mat'));
    flag_co=1;
else
    flag_co=0;
end

% define 4 corner-points of small area
if local_global==1 % local system
    x_corner=[corner(1) corner(1) corner(1)+width corner(1)+width]; % x-coords of corners
    y_corner=[corner(2) corner(2)+height corner(2)+height corner(2)]; % y-coords of corners
    % coordinate grids:
    [x,y]=meshgrid(x_corner(1):dx:x_corner(3),y_corner(1):dx:y_corner(2));
else % global system (Easting/Northing)
    if flag_co==0 % area not rotated -> global system
        x_corner=[corner(1) corner(1) corner(1)+width corner(1)+width]; % Easting-coords of corners
        y_corner=[corner(2) corner(2)+height corner(2)+height corner(2)]; % Northing-coords of corners
        % coordinate grids:
        [x,y]=meshgrid(x_corner(1):dx:x_corner(3),y_corner(1):dx:y_corner(2));
    else % area rotated and width/heigth in global system
        x_corner=[corner(1) corner(1) corner(1)+width corner(1)+width]; % Easting-coords of corners
        y_corner=[corner(2) corner(2)+height corner(2)+height corner(2)]; % Northing-coords of corners
        % coordinate grids:
        [x,y]=meshgrid(min(x_corner):dx:max(x_corner),min(y_corner):dx:max(y_corner)); % Easting/northing!
    end
end
xrg=x(1,1)-dx/2:dx:x(1,end)+dx/2;
yrg=y(1,1)-dx/2:dx:y(end,1)+dx/2;

if flag_co==1
    % convert xmin/xmax/ymin/ymax of rectangles to global system
    temp=helmert([xmin' ymin'; xmax' ymax'],coordtrans(:,1:2),coordtrans(:,3:4)); % convert local2global
    xmin=temp(1:length(xmin),1)';
    xmax=temp(length(xmin)+1:end,1);
    ymin=temp(1:length(xmin),2)';
    ymax=temp(length(xmin)+1:end,2);
end


% open rectangle data:
for j=1:length(rectangles)
    if proc==1
        matFileObj{j}=matfile(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'processed','data.mat'));
    else
        matFileObj{j}=matfile(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'data.mat'));
    end
end

disp('Reading data... Please wait.');
% read data and put in small area
for j=1:length(rectangles)
    disp(['Rectangle ',int2str(j),'/',int2str(length(rectangles))])
    
    if any(x(1,:)>=xmin(j)) && any(x(1,:)<=xmax(j)) && any(y(:,1)>=ymin(j)) && any(y(:,1)<=ymax(j))
        disp('   Found data in small area...') 
        tmp_data=reshape(permute(matFileObj{j}.data,[3 1 2]),length(t),[]);
        [r,c,~]=size(tmp_data);
        c=c/length(t);
        ns=length(t);
        ok=~isnan(tmp_data(1,:));
        tmp_data=tmp_data(:,ok); % only valid traces
        tmp_x=reshape(xx{j},r*c,1);
        tmp_y=reshape(yy{j},r*c,1);
        tmp_x=tmp_x(ok);
        tmp_y=tmp_y(ok);
        if flag_co==1
            co=helmert([tmp_x tmp_y],coordtrans(:,1:2),coordtrans(:,3:4));
            tmp_x=co(:,1);
            tmp_y=co(:,2);
        end
        temp=bindata3(tmp_data,tmp_x,tmp_y,xrg,yrg);
        if ~exist('data','var')
            data=zeros(size(temp),'single');
            anz=zeros(size(temp),'single');
            topo=zeros(size(x));
        end
        anz(~isnan(temp))=anz(~isnan(temp))+1;
        temp(isnan(temp))=0;
        data=data+temp;

        % topography
        if proc==1
            load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'processed','z.mat'));
            load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'mask','z.mat'));
        else
            load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'z.mat'));
            load(fullfile(pfad_rad,['3D_Grid_R',int2str(rectangles(j))],'mask.mat'));
        end
        tmp_x=xx{j}(mask==1 & ~isnan(z));
        tmp_y=yy{j}(mask==1 & ~isnan(z));
        if flag_co==1
            co=helmert([tmp_x tmp_y],coordtrans(:,1:2),coordtrans(:,3:4));
            tmp_x=co(:,1);
            tmp_y=co(:,2);
        end
        zt=bindata2(z(mask==1 & ~isnan(z)),tmp_x,tmp_y,xrg,yrg);
        zt(isnan(zt))=0;
        topo=topo+zt;
    end
end
data=data./anz;
topo=topo./anz(:,:,1);
z=topo;

if interp==1
    disp(' ')
    disp('Interpolating...')
    xt=x(anz(:,:,1)~=0);
    yt=y(anz(:,:,1)~=0);
    ok=anz(:,:,1)~=0;
    for i=1:ns
        if ~mod(i,20)
            disp([int2str(i),'/',int2str(ns)])
        end
        temp=data(:,:,i);
        F=scatteredInterpolant(xt,yt,double(temp(ok)));
        data(:,:,i)=F(x,y);
    end
    F=scatteredInterpolant(xt,yt,double(topo(ok)));
    z=F(x,y);
end

% mask
mask=zeros(size(x));
mask(anz(:,:,1)>0)=1;
if interp==1
    temp=ones(size(mask));
    temp(mask==1)=0;
    eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
    eucmap_interp=griddata(x(:),y(:),eucmap(:),x,y);
    eucmap_interp(isnan(eucmap_interp))=2*radius;
    mask_interp=ones(size(eucmap_interp)); % initialize new grid
    mask_interp(eucmap_interp>radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
    mask=mask_interp;
end


% Save small area
if proc==0
    if ~exist(fullfile(pfad_rad,'Small_3Dblock'))
        mkdir(fullfile(pfad_rad,'Small_3Dblock'));
        if ~exist(fullfile(pfad_rad,'Small_3Dblock','3D_Grid_R1'))
            mkdir(fullfile(pfad_rad,'Small_3Dblock','3D_Grid_R1'))
        end
    end
    save(fullfile(pfad_rad,'Small_3Dblock','data.mat'),'data','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock','mask.mat'),'mask','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock','t.mat'),'t','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock','x.mat'),'x','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock','y.mat'),'y','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock','z.mat'),'z','-v7.3');
else
    if ~exist(fullfile(pfad_rad,'Small_3Dblock_proc'))
        mkdir(fullfile(pfad_rad,'Small_3Dblock_proc'));
        if ~exist(fullfile(pfad_rad,'Small_3Dblock_proc','3D_Grid_R1'))
            mkdir(fullfile(pfad_rad,'Small_3Dblock_proc','3D_Grid_R1'))
        end
    end
    save(fullfile(pfad_rad,'Small_3Dblock_proc','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock_proc','mask.mat'),'mask','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock_proc','t.mat'),'t','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock_proc','x.mat'),'x','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock_proc','y.mat'),'y','-v7.3');
    save(fullfile(pfad_rad,'Small_3Dblock_proc','z.mat'),'z','-v7.3');
end

% set original path
path(oldpath);