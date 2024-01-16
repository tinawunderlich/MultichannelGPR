% Script for joining Timeslices from e.g. measurements of several days and
% adjust them before joining
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files:
% ../day1/Timeslices and ../day2/Timeslices folders

clear all
close all
clc


% only use one of the following (or set both to 0!):
sub_mean=1;     % if ==1: subtract mean for each rectangle for each time sample (to reduce offsets)
sub_median=0;   % if ==1: subtract median for each rectangle for each time sample (to reduce offsets)

% masking options
nn_radius=2;    % radius to nearest neighbor (in bins) should be less than nn_radius to be valid


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

if ispc
    if exist('temp1.temp') % read last opened folder from temp.temp
        fid=fopen('temp1.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname1=uigetdir(fn{1}{1},'Choose folder with timeslices of day 1');
        else
            pathname1=uigetdir([],'Choose folder with timeslices of day 1');
        end
        fileattrib('temp1.temp','-h');
        fid=fopen('temp1.temp','wt');
        fprintf(fid,'%s',pathname1);
        fclose(fid);
        fileattrib('temp1.temp','+h');
    else
        pathname1=uigetdir([],'Choose folder with timeslices of day 1'); % path to radargram-folder

        fid=fopen('temp1.temp','wt');
        fprintf(fid,'%s',pathname1);
        fclose(fid);
        fileattrib('temp1.temp','+h');
    end

    if exist('temp2.temp') % read last opened folder from temp.temp
        fid=fopen('temp2.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname2=uigetdir(fn{1}{1},'Choose folder with timeslices of day 2');
        else
            pathname2=uigetdir([],'Choose folder with timeslices of day 2');
        end
        fileattrib('temp2.temp','-h');
        fid=fopen('temp2.temp','wt');
        fprintf(fid,'%s',pathname2);
        fclose(fid);
        fileattrib('temp2.temp','+h');
    else
        pathname2=uigetdir([],'Choose folder with timeslices of day 2'); % path to radargram-folder

        fid=fopen('temp2.temp','wt');
        fprintf(fid,'%s',pathname2);
        fclose(fid);
        fileattrib('temp2.temp','+h');
    end
else
    if exist('.temp1.temp') % read last opened folder from temp.temp
        fid=fopen('.temp1.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname1=uigetdir(fn{1}{1},'Choose folder with timeslices of day 1');
        else
            pathname1=uigetdir([],'Choose folder with timeslices of day 1');
        end
    else
        pathname1=uigetdir([],'Choose folder with timeslices of day 1'); % path to radargram-folder
    end

    fid=fopen('.temp1.temp','wt');
    fprintf(fid,'%s',pathname1);
    fclose(fid);

     if exist('.temp2.temp') % read last opened folder from temp.temp
        fid=fopen('.temp2.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname2=uigetdir(fn{1}{1},'Choose folder with timeslices of day 2');
        else
            pathname2=uigetdir([],'Choose folder with timeslices of day 2');
        end
    else
        pathname2=uigetdir([],'Choose folder with timeslices of day 2'); % path to radargram-folder
    end

    fid=fopen('.temp2.temp','wt');
    fprintf(fid,'%s',pathname2);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));


%% load Tsl day1:
load(fullfile(pathname1,'t_startende.mat'));
load(fullfile(pathname1,'coordtrans.mat'));
load(fullfile(pathname1,'mask.mat'));
load(fullfile(pathname1,'tsl.mat'));
load(fullfile(pathname1,'xgrid.mat'));
load(fullfile(pathname1,'ygrid.mat'));
load(fullfile(pathname1,'topo.mat'));

% transform coordinates:
xy1=helmert([xgrid(:) ygrid(:)],coordtrans(:,1:2),coordtrans(:,3:4));
maske1=mask;
tsl1=tsl;
topo1=topo;
startende1=t_startende;
dx1=xgrid(1,2)-xgrid(1,1);

%% load Tsl day2:
load(fullfile(pathname2,'t_startende.mat'));
load(fullfile(pathname2,'coordtrans.mat'));
load(fullfile(pathname2,'mask.mat'));
load(fullfile(pathname2,'tsl.mat'));
load(fullfile(pathname2,'xgrid.mat'));
load(fullfile(pathname2,'ygrid.mat'));
load(fullfile(pathname2,'topo.mat'));

% transform coordinates:
xy2=helmert([xgrid(:) ygrid(:)],coordtrans(:,1:2),coordtrans(:,3:4));
maske2=mask;
tsl2=tsl;
topo2=topo;
startende2=t_startende;
dx2=xgrid(1,2)-xgrid(1,1);

%% make new coordinate grid:
dx=min([dx1, dx2]);
[xgrid,ygrid]=meshgrid([min([xy1(:,1); xy2(:,1)]):dx:max([xy1(:,1); xy2(:,1)])],[min([xy1(:,2); xy2(:,2)]):dx:max([xy1(:,2); xy2(:,2)])]);
xgrid_interp=xgrid;
ygrid_interp=ygrid;

%% check if both Tsl have same times:
if ~all(startende1==startende2)
    disp('Timeslices do not have same time steps. Please correct and start again.');
    return;
end

%% adjust and join for each tsl
xy=[xy1; xy2];
for i=1:length(tsl1)
    % determine mean and median
    me1(i)=mean(tsl1{i}(~isnan(tsl1{i})));
    me2(i)=mean(tsl2{i}(~isnan(tsl2{i})));
    md1(i)=median(tsl1{i}(~isnan(tsl1{i})));
    md2(i)=median(tsl2{i}(~isnan(tsl2{i})));
    
    disp(['Tsl ',int2str(i),' of ',int2str(length(tsl1)),'...'])
    
    if sub_mean==1
        data=[tsl1{i}(:)-me1(i); tsl2{i}(:)-me2(i)];
    elseif sub_median==1
        data=[tsl1{i}(:)-md1(i); tsl2{i}(:)-md2(i)];
    else
        data=[tsl1{i}(:); tsl2{i}(:)]; % no adjustment
    end
    
    % delete nan before binning
    xytemp=xy(~isnan(data),:);
    data(isnan(data))=[];

    
    % binning:
    tsl{i}=bindata2(data,xytemp(:,1),xytemp(:,2),xgrid(1,1)-dx/2:dx:xgrid(1,end)+dx/2,ygrid(1,1)-dx/2:dx:ygrid(end,1)+dx/2);
    
    % interpolation:
    tsl_interp{i}=griddata(xytemp(:,1),xytemp(:,2),data,xgrid,ygrid);
end

% join mask:
disp('Join masks...')
datam=[maske1(:); maske2(:)];
% binning:
mask=bindata2(datam,xy(:,1),xy(:,2),xgrid(1,1)-dx/2:dx:xgrid(1,end)+dx/2,ygrid(1,1)-dx/2:dx:ygrid(end,1)+dx/2);
mask(mask>0)=1;
% Make mask for interpolated area
temp=ones(size(tsl{1}));
temp(~isnan(tsl{1}))=0;
eucmap=chamfer_DT(temp);  % approximated euclidian distance map (Distance to next neighbor in bins)
eucmap_interp=griddata(xgrid(:),ygrid(:),eucmap(:),xgrid,ygrid);
eucmap_interp(isnan(eucmap_interp))=2*nn_radius;
mask_interp=ones(size(eucmap_interp)); % initialize new grid
mask_interp(eucmap_interp>nn_radius)=0;   % set 0 for pixels with distance to nearest neighbor > radius

% apply mask to interpolated data
for i=1:length(tsl)
    tsl_interp{i}(mask_interp==0)=NaN;
end

%% coordtrans is not used, but stored as unity
coordtrans=ones(2,4);

%% topo
disp('Join Topography...')
datat=[topo1(:); topo2(:)];
% delete nan before binning
xytemp=xy(~isnan(datat),:);
datat(isnan(datat))=[];
% binning:
topo=bindata2(datat,xytemp(:,1),xytemp(:,2),xgrid(1,1)-dx/2:dx:xgrid(1,end)+dx/2,ygrid(1,1)-dx/2:dx:ygrid(end,1)+dx/2);
% interpolation:
topo_interp=griddata(xytemp(:,1),xytemp(:,2),datat,xgrid,ygrid);
topo_interp(mask_interp==0)=NaN;
    
%% write info file about processing:
if ~exist(fullfile(pathname1,'Tsl_Day1_Day2'),'dir')
    mkdir(fullfile(pathname1,'Tsl_Day1_Day2'));
end
fid=fopen(fullfile(pathname1,'Tsl_Day1_Day2','info.txt'),'wt');
fprintf(fid,['sub_mean=',int2str(sub_mean),'\nsub_median=',int2str(sub_median),'\n']);
fprintf(fid,'No.\tmean_day1\tmean_day2\tmedian_day1\tmedian_day2\n');
fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',[[1:length(me1)]' me1' me2' md1' md2']');
fclose(fid);

disp('Saving data...')
%% write data (uninterpolated)
save(fullfile(pathname1,'Tsl_Day1_Day2','xgrid.mat'),'xgrid','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','ygrid.mat'),'ygrid','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','t_startende.mat'),'t_startende','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','mask.mat'),'mask','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','tsl.mat'),'tsl','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','coordtrans.mat'),'coordtrans','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','topo.mat'),'topo','-v7.3');

%% write data (interpolated)
if ~exist(fullfile(pathname1,'Tsl_Day1_Day2','interpolated'),'dir')
    mkdir(fullfile(pathname1,'Tsl_Day1_Day2','interpolated'));
end
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','xgrid_interp.mat'),'xgrid_interp','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','ygrid_interp.mat'),'ygrid_interp','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','t_startende.mat'),'t_startende','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','mask_interp.mat'),'mask_interp','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','tsl_interp.mat'),'tsl_interp','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','coordtrans.mat'),'coordtrans','-v7.3');
save(fullfile(pathname1,'Tsl_Day1_Day2','interpolated','topo_interp.mat'),'topo_interp','-v7.3');

% set original path
path(oldpath);