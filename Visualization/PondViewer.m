clear all
close all
clc

% Use the binned 3D data blocks and visualize them with PondView (Grasmueck & Viggiano 2018)
% M. Grasmueck and D. Viggiano, "PondView: Intuitive and Efficient Visualization of 3D GPR Data," 2018 17th International Conference on Ground Penetrating Radar (GPR), Rapperswil, Switzerland, 2018, pp. 1-6, doi: 10.1109/ICGPR.2018.8441634.
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires PondView_GUI in folder GUIs and Subfunctions


% time or depth domain?
tz_flag=2; % =1: time, =2: depth

rectangles=[1:2]; % number of rectangles

radius=0.3; % radius in m for valid interpolation

use_slices=1;   % if =0: make new slices and save them,
% if =1: use already created slices in folder
% PondViewSlices (but nevertheless select folder containing
% the 3D_Grid... and PondViewSlices folders)

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning off

% path for slices/3D_Grids
if ~ispc; menu('Choose folder containing 3D_Grid_R*-folders','OK'); end
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Select folder containing 3D_Grid_R*-folders');
        else
            pfad=uigetdir([],'Select folder containing 3D_Grid_R*-folders');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Select folder containing 3D_Grid_R*-folders');

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
            pfad=uigetdir(fn{1}{1},'Select folder containing 3D_Grid_R*-folders');
        else
            pfad=uigetdir([],'Select folder containing 3D_Grid_R*-folders');
        end
    else
        pfad=uigetdir([],'Select folder containing 3D_Grid_R*-folders');
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% path for radargrams
if ~ispc; menu('Choose folder with radargrams','OK'); end
if ispc
    if exist('rtemp.temp') % read last opened folder from temp.temp
        fid=fopen('rtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Select folder containing radargrams.mat');
        else
            pfad_rad=uigetdir([],'Select folder containing radargrams.mat');
        end
        fileattrib('rtemp.temp','-h');
        fid=fopen('rtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('rtemp.temp','+h');
    else
        pfad_rad=uigetdir([],'Select folder containing radargrams.mat');

        fid=fopen('rtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('rtemp.temp','+h');
    end
else
    if exist('.rtemp.temp') % read last opened folder from temp.temp
        fid=fopen('.rtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Select folder containing radargrams.mat');
        else
            pfad_rad=uigetdir([],'Select folder containing radargrams.mat');
        end
    else
        pfad_rad=uigetdir([],'Select folder containing radargrams.mat');
    end

    fid=fopen('.rtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'),fullfile(curFold,'Subfunctions'));

%% Radargrams
temp=load(fullfile(pfad_rad,'global_coords.mat'));
coords=temp.global_coords;
temp=load(fullfile(pfad_rad,'radargrams.mat'));
radar=temp.radargrams;
temp=load(fullfile(pfad_rad,'t.mat'));
tr=temp.t;



%% SLices
if use_slices==0 % load new data and create slices
    %% load data
    if tz_flag==2 % depth
        if exist(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'depth.mat'),'file')
            temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'depth.mat'));
            t = temp.depth'; % depth is starting from 0, positive down
            maxElevation = temp.maxElevation;
        else % depth.mat does not exist, but depth is in t.mat (absolute depths!)
            t=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
            t=t.t';
            maxElevation=t(1); % absolute height at top
            t=abs(t-t(1)); % depth is starting from 0, positive down
        end
        fprintf("Depth vector boundaries in [m] are [%7.2f %7.2f] \n",min(t), max(t));
    else
        load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
        maxElevation=[];
        fprintf("Time vector boundaries in [ns] are [%7.2f %7.2f] \n",min(t), max(t));
    end

    % read coordinates
    for i=1:length(rectangles)

        p=['3D_Grid_R',int2str(rectangles(i))];

        temp=load(fullfile(pfad,p,'x.mat'));
        x{i}=temp.x;
        temp=load(fullfile(pfad,p,'y.mat'));
        y{i}=temp.y;
        temp=load(fullfile(pfad,p,'z.mat'));
        z{i}=temp.z;
        temp=load(fullfile(pfad,p,'mask.mat'));
        mask{i}=temp.mask;

        if i==1
            dx=x{i}(1,2)-x{i}(1,1);   % dx of 3D-block
        end

        xmin(i)=min(x{i}(1,:));
        xmax(i)=max(x{i}(1,:));
        ymin(i)=min(y{i}(:,1));
        ymax(i)=max(y{i}(:,1));
    end
    clear temp;

    % make grid for all bins
    [xgrid,ygrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)]);

    disp('Reading data... Please wait.');

    n=length(t);
    slice=NaN(size(xgrid));    % initialize tsl in first run
    for i=1:n   % for each sample -> make sampleslice
        if ~mod(i,50)
            disp(['Slice ',int2str(i),'/',int2str(n)])
        end

        for j=1:length(rectangles)  % in each rectangle...
            % open matfile
            p=['3D_Grid_R',int2str(rectangles(j)),'/data.mat'];
            matFileObj=matfile(fullfile(pfad,p));

            if i==1 % for first slice only
                % determine location in big area (find lower left corner)
                row_ind(j)=find(abs(ygrid(:,1)-ymin(j))==min(abs(ymin(j)-ygrid(:,1))));
                col_ind(j)=find(abs(xgrid(1,:)-xmin(j))==min(abs(xmin(j)-xgrid(1,:))));
            end

            datatemp=abs(matFileObj.data(:,:,i)); % abs(data) in this sample slice

            % put tsl_temp together in big tsl
            slice(row_ind(j):row_ind(j)+length(datatemp(:,1))-1,col_ind(j):col_ind(j)+length(datatemp(1,:))-1)=datatemp;
            slice=slice(1:length(xgrid(:,1)),1:length(xgrid(1,:)));   % resize if necessary
        end

        % Make mask
        if i==1 && tz_flag==1 % TIMEslice
            mask=zeros(size(slice));
            mask(~isnan(slice))=1;

            temp=ones(size(mask));
            temp(mask==1)=0;
            eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
            mask_interp=ones(size(eucmap)); % initialize new grid
            mask_interp(eucmap.*dx>radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
        elseif tz_flag==2 % DEPTHslices -> combine masks into maximum mask
            if i==1
                mask=zeros(size(slice)); % initialize mask
            end
            mask(~isnan(slice))=1;

            temp=ones(size(mask));
            temp(mask>=1)=0;
            eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
            mask_interp=ones(size(eucmap)); % initialize new grid
            mask_interp(eucmap.*dx>radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
        end


        % Interpolation
        if length(slice(~isnan(slice(:))))>= 3 % at least 3 values are necessary for triangulation
            F=scatteredInterpolant(xgrid(mask>0),ygrid(mask>0),slice(mask>0));
            slice=reshape(F(xgrid(:),ygrid(:)),size(xgrid));
            if isempty(slice) % problem if only some points and are in line (e.g. on one x value)
                slice=NaN(size(xgrid));
            else
                slice(mask_interp==0)=NaN;
            end
        else
            slice=NaN(size(xgrid));
        end

        % save slice
        if ~exist(fullfile(pfad,'PondViewSlices'),'dir')
            mkdir(fullfile(pfad,'PondViewSlices'))
        end
        save(fullfile(pfad,'PondViewSlices',['slice_',int2str(i),'.mat']),'slice','-v7.3');
    end
    if tz_flag==2
        mask(mask>=1)=1;
    end
    mask=logical(mask);

    save(fullfile(pfad,'PondViewSlices','mask_interp.mat'),'mask_interp','-v7.3');
    save(fullfile(pfad,'PondViewSlices','t.mat'),'t','-v7.3');
    save(fullfile(pfad,'PondViewSlices','maxElevation.mat'),'maxElevation','-v7.3');
    save(fullfile(pfad,'PondViewSlices','xgrid.mat'),'xgrid','-v7.3');
    save(fullfile(pfad,'PondViewSlices','ygrid.mat'),'ygrid','-v7.3');

    % write info file:
    fid=fopen(fullfile(pfad,'PondViewSlices','PondView_info.txt'),'wt');
    fprintf(fid,'Number of samples/slices: %d\n',length(t));
    fprintf(fid,'tz_flag: %d\n',tz_flag);
    fprintf(fid,'Number of rectangles: %d\n',length(rectangles));
    fprintf(fid,'Radius for valid interpolation: %.1f m',radius);
    fclose(fid);

else % if slices can be read
    if ~exist(fullfile(pfad,'PondViewSlices','slice_1.mat'),'file')
        disp('No slices found. Please check path and settings.')
        return
    else
        disp('Slices found.')
    end
end

% Read info on slices:
list=dir(fullfile(pfad,'PondViewSlices','slice_*.mat'));
load(fullfile(pfad,'PondViewSlices','t.mat'));
load(fullfile(pfad,'PondViewSlices','mask_interp.mat'));
load(fullfile(pfad,'PondViewSlices','maxElevation.mat'));
load(fullfile(pfad,'PondViewSlices','xgrid.mat'));
load(fullfile(pfad,'PondViewSlices','ygrid.mat'));


PondView_GUI(xgrid,ygrid,list,mask_interp,t,tz_flag,maxElevation,tr,radar,coords,fullfile(pfad,'PondViewSlices'));

waitfor(gcf);

% set original path
path(oldpath);