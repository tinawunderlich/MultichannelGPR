clear all
close all
clc


% If timeslices have been created with make_Timeslices and they were not 
% properly interpolated, you can re-interpolate them here without the need 
% to run make_Timeslices again
% 
% Select the Timeslices folder!
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%


% depth slices instead of time slices (input is in m)
dsl = 0; % =1: depth, =0: time

% masking options
nn_radius=3;    % radius to nearest neighbor (in bins) should be less than nn_radius to be valid

% Interpolation
griding=1;  % 1: Griddata (linear interpolation)
% 2: Inverse distance weighting =IDW (use also the following
% parameters: radius and power)
radius=6;   % Radius for IDW and masking of interpolated timeslice (in bins)
power=10;    % Power of IDW



%% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose Timeslices folder');
        else
            pfad=uigetdir([],'Choose Timeslices folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Choose Timeslices folder'); % path to radargram-folder

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
            pfad=uigetdir(fn{1}{1},'Choose Timeslices folder');
        else
            pfad=uigetdir([],'Choose Timeslices folder');
        end
    else
        pfad=uigetdir([],'Choose Timeslices folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));


% read data
disp('Reading data...')
load(fullfile(pfad,'mask.mat'));
load(fullfile(pfad,'topo.mat'));
load(fullfile(pfad,'tsl.mat'));
load(fullfile(pfad,'xgrid.mat'));
load(fullfile(pfad,'ygrid.mat'));
load(fullfile(pfad,'interpolated','xgrid_interp.mat'))
load(fullfile(pfad,'interpolated','ygrid_interp.mat'))

dx_tsl=abs(xgrid_interp(1,2)-xgrid_interp(1,1)); % dx of interpolated tsl
dx=abs(xgrid(1,2)-xgrid(1,1)); % dx of tsl

% Make mask for interpolated area
if dsl==0 % TIMEslice -> just one mask for all slices
    temp=ones(size(tsl{1}));
    temp(mask==1)=0;
    eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
    eucmap_interp=griddata(double(xgrid(:)),double(ygrid(:)),eucmap(:),double(xgrid_interp),double(ygrid_interp));
    eucmap_interp(isnan(eucmap_interp))=2*nn_radius;
    mask_interp=ones(size(eucmap_interp)); % initialize new grid
    mask_interp(eucmap_interp>nn_radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
    mask_topointerp=mask_interp;
else % DEPTHslices -> one mask for each slice
    for i=1:length(tsl)
        temp=ones(size(tsl{i}));
        temp(mask{i}==1)=0;
        eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
        eucmap_interp=griddata(double(xgrid(:)),double(ygrid(:)),eucmap(:),double(xgrid_interp),double(ygrid_interp));
        eucmap_interp(isnan(eucmap_interp))=2*nn_radius;
        mask_interp{i}=ones(size(eucmap_interp)); % initialize new grid
        mask_interp{i}(eucmap_interp>nn_radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
    end
    % topo:
    temp=ones(size(tsl{1}));
    temp(mask_topo==1)=0;
    eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
    eucmap_interp=griddata(double(xgrid(:)),double(ygrid(:)),eucmap(:),double(xgrid_interp),double(ygrid_interp));
    eucmap_interp(isnan(eucmap_interp))=2*nn_radius;
    mask_topointerp=ones(size(eucmap_interp)); % initialize new grid
    mask_topointerp(eucmap_interp>nn_radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius
end


% Interpolation
disp('Interpolation... Please wait.');
if griding==1
    h=waitbar(0,'Linear interpolation is running...');
    for i=1:length(tsl) % for each timeslice
        if dsl==0
            maske=mask_interp;
        else
            maske=mask_interp{i};
        end
        
        if i==1; c1=clock; end
        if length(tsl{i}(~isnan(tsl{i}(:))))>= 3 % at least 3 values are necessary for triangulation   
            tsl_interp{i}=griddata(double(xgrid(~isnan(tsl{i}(:)))),double(ygrid(~isnan(tsl{i}(:)))),tsl{i}(~isnan(tsl{i}(:))),double(xgrid_interp),double(ygrid_interp));
            if isempty(tsl_interp{i}) % problem if only some points and are in line (e.g. on one x value)
                tsl_interp{i}=NaN(size(xgrid_interp),'single');
            else
                tsl_interp{i}(~maske)=NaN;   % apply mask_interp
            end
        else
            tsl_interp{i}=NaN(size(xgrid_interp),'single');
        end
        if i==1
            c2=clock;
            diff=minutes(datetime(c2)-datetime(c1));    % time for one run in minutes
        end
        waitbar(i/length(tsl),h,['Approx. ',int2str(diff*length(tsl)-diff*i),' minutes remaining']);
    end
    topo_interp=griddata(double(xgrid(~isnan(topo(:)))),double(ygrid(~isnan(topo(:)))),topo(~isnan(topo(:))),double(xgrid_interp),double(ygrid_interp));
    topo_interp(~mask_topointerp)=NaN;   % apply mask_interp
    close(h);
elseif griding==2
    if dx~=dx_tsl
        % enlarge grid before idw
        topo=bindata2(topo(~isnan(topo)),xgrid(~isnan(topo)),ygrid(~isnan(topo)),linspace(xgrid_interp(1,1)-dx_tsl/2,xgrid_interp(1,end)+dx_tsl/2,length(xgrid_interp(1,:))+1),linspace(ygrid(1,1)-dx_tsl/2,ygrid(end,1)+dx_tsl/2,length(xgrid_interp(:,1))+1));
        for i=1:length(tsl)
            tsl{i}=bindata2(tsl{i}(~isnan(tsl{i})),xgrid(~isnan(tsl{i})),ygrid(~isnan(tsl{i})),linspace(xgrid_interp(1,1)-dx_tsl/2,xgrid_interp(1,end)+dx_tsl/2,length(xgrid_interp(1,:))+1),linspace(ygrid(1,1)-dx_tsl/2,ygrid(end,1)+dx_tsl/2,length(xgrid_interp(:,1))+1));
        end
        % apply idw
        tsl_interp=idw2d_tsl(tsl,xgrid_interp,ygrid_interp,eucmap_interp,radius,power);  % for all timeslices together in function, mask is applied automatically
        topo_interp=idw2d_tsl(topo,xgrid_interp,ygrid_interp,eucmap_interp,radius,power);
    else
        tsl_interp=idw2d_tsl(tsl,xgrid,ygrid,eucmap,radius,power);  % for all timeslices together in function, mask is applied automatically
        topo_interp=idw2d_tsl(topo,xgrid,ygrid,eucmap,radius,power);
    end
    temp=topo_interp{1};
    delete topo_interp;
    topo_interp=temp;
end


%% save data
disp('Saving data...')

% save interpolated timeslices
save(fullfile(pfad,'interpolated','tsl_interp.mat'),'tsl_interp','-v7.3');
save(fullfile(pfad,'interpolated','topo_interp.mat'),'topo_interp','-v7.3');
save(fullfile(pfad,'interpolated','mask_interp.mat'),'mask_interp','-v7.3');

disp('Finished!')

% set original path
path(oldpath);
