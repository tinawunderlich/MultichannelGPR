clear all
close all
clc


% Prepare topo.mat and vgrid.mat for Migrations in Processing.m
%
% Dr. Tina Wunderlich, CAU Kiel 2021-2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires folder with radargrams.mat/global_coords.mat/x.mat/t.mat
% -> choose this folder when asked!


% Choose options for velocity model:
vopt=1;
% vopt=1: Constant velocity for all profiles in radargrams.mat
vconst=0.1; % v in m/ns
% vopt=2: Constant velocity (but different) for each profile in radargrams.mat
vall=[0.1 0.08]; % v for each profile in m/ns
% vopt=3: 1D velocity model for all profiles
v1d=[0.16 0.1]; % v at different times in m/ns
t1d=[0 18]; % corresponding times in ns
% vopt=4: 1D velocity model for all profiles from fitted function in
    % vrms.mat and tv.mat (in same folder as your data) (from Make_1Dv.m)


% Choose options for topography:
topoopt=1;
% topoopt=0: no topography file required
% topoopt=1: topography is already set in global_coords(:,3)
% topoopt=2: topography has to be set with file containing profile number (1. column),
% profile coordinate (2. Column) and topography (3. column), all in m, This
% file should be in the same folder as radargrams.mat
% topoopt=3: the topography is known at several points on the measured
% area, but not directly on the profiles. The topofile contains 3 columns:
% Easting, Northing and height
topofile='topofile.txt'; % file for topoopt==2 or topoopt==3
smooth_topo=0; % n>0: smooth topography over n samples, if no smoothing=0

removeOutliers=0; % if =1: remove outliers in topo data (is done before smoothing)
num=101; % odd(!) number of points for median calculation if removeOutliers==1
thresh=0.1; % threshold in m for difference between median and raw topography if removeOutliers==1

plottopo=0; % if =1: topo is plotted (raw and smoothed), if =0: no figures

%% -------------------------------------------------------------------------
% Do not change the following part!

% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with radargrams','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
        
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
    end
else
    if exist('.radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('.radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end
    
    fid=fopen('.radtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end

load(fullfile(pfad_rad,'global_coords.mat'));
load(fullfile(pfad_rad,'x.mat'));
load(fullfile(pfad_rad,'t.mat'));
if length(t(1,:))>1
    t=t'; % make column vector
end


%% prepare velocity file
if vopt==1
    for i=1:length(global_coords)
        v{i}=vconst;
    end
    maxdepth=max(t)/2*vconst; % maximum depth of radargram
elseif vopt==2
    for i=1:length(global_coords)
        v{i}=vall(i);
    end
    maxdepth=max(max(t)./2.*vall); % maximum depth of radargram
elseif vopt==3
    for i=1:length(global_coords)
        v{i}=repmat(interp1(t1d,v1d,t,'linear',v1d(end)),[1 length(global_coords{i}(:,1))]);
        % get depth for each profile
        d(i)=max(t)/2*v{i}(end);
    end
    maxdepth=max(d); % maximum depth of radargram
elseif vopt==4
    load(fullfile(pfad_rad,'vrms.mat'));
    if vrms>1 % if in cm/ns
        vrms=vrms./100;
    end
    load(fullfile(pfad_rad,'tv.mat'));
    for i=1:length(global_coords)
        v{i}=repmat(interp1(tv,vrms,t,'linear','extrap'),[1 length(global_coords{i}(:,1))]);
        % get depth for each profile
        d(i)=max(t)/2*v{i}(end);
    end
    maxdepth=max(d); % maximum depth of radargram
end
save(fullfile(pfad_rad,'vgrid.mat'),'v');


%% prepare topography file
if topoopt>0
    if topoopt==1
        if length(global_coords{i}(1,:))==2 % only x and y
            disp('global_coords does not contain topography. Please use topoopt=2 instead.')
            return;
        else
            if smooth_topo>0
                for i=1:length(global_coords)
                    gc=global_coords{i}(:,3);
                    if plottopo==1
                        figure
                        plot(gc,'Marker','*')
                        hold on
                        xlabel('Trace number')
                        ylabel('Topography [m]')
                        led=[{'raw data'}];
                        legend(led)
                        grid on
                        title(['Profile ',int2str(i)])
                    end
                    
                    if removeOutliers==1
                        gc=remout(gc,num,thresh);
                        if plottopo==1
                            plot(gc,'Marker','*')
                            led=[led; {'after removeOutliers'}];
                            legend(led)
                            grid on
                        end
                    end
                    topo{i}=smooth(gc,smooth_topo);
                    if plottopo==1
                        plot(topo{i},'Marker','*')
                        led=[led; {'after smoothing'}];
                        legend(led)
                        grid on
                    end
                end
            else
                for i=1:length(global_coords)
                    gc=global_coords{i}(:,3);
                    if plottopo==1
                        figure
                        plot(gc,'Marker','*')
                        hold on
                        xlabel('Trace number')
                        ylabel('Topography [m]')
                        led=[{'raw data'}];
                        legend(led)
                        grid on
                        title(['Profile ',int2str(i)])
                    end
                    if removeOutliers==1
                        gc=remout(gc,num,thresh);
                        if plottopo==1
                            plot(gc,'Marker','*')
                            led=[led; {'after removeOutliers'}];
                            legend(led)
                            grid on
                        end
                    end
                    topo{i}=gc;
                end
            end
        end
    elseif topoopt==2
        top=load(fullfile(pfad_rad,topofile));
        if smooth_topo>0
            for i=1:length(global_coords)
                topo{i}=interp1(top(top(:,1)==i,2),top(top(:,1)==i,3),x{i});
                if plottopo==1
                    figure
                    plot(topo{i},'Marker','*')
                    hold on
                    xlabel('Trace number')
                    ylabel('Topography [m]')
                    led=[{'raw data'}];
                    legend(led)
                    grid on
                    title(['Profile ',int2str(i)])
                end
                topo{i}=smooth(topo{i},smooth_topo);
                if plottopo==1
                    plot(topo{i},'Marker','*')
                    led=[led; {'after smoothing'}];
                    legend(led)
                    grid on
                end
            end
        else
            for i=1:length(global_coords)
                topo{i}=interp1(top(top(:,1)==i,2),top(top(:,1)==i,3),x{i});
                if plottopo==1
                    figure
                    plot(topo{i},'Marker','*')
                    hold on
                    xlabel('Trace number')
                    ylabel('Topography [m]')
                    led=[{'raw data'}];
                    legend(led)
                    grid on
                    title(['Profile ',int2str(i)])
                end
            end
        end
        % set topo in global_coords:
        for i=1:length(global_coords)
            global_coords{i}(:,3)=topo{i};
        end
        save(fullfile(pfad_rad,'global_coords.mat'),'global_coords');
   elseif topoopt==3
        top=load(fullfile(pfad_rad,topofile));
        % interpolate topography from area to radargrams:
        F=scatteredInterpolant(top(:,1),top(:,2),top(:,3),'linear','linear');
        if smooth_topo>0
            for i=1:length(global_coords)
                topo{i}=F(global_coords{i}(:,1),global_coords{i}(:,2));
                if plottopo==1
                    figure
                    plot(topo{i},'Marker','*')
                    hold on
                    xlabel('Trace number')
                    ylabel('Topography [m]')
                    led=[{'raw data'}];
                    legend(led)
                    grid on
                    title(['Profile ',int2str(i)])
                end
                topo{i}=smooth(topo{i},smooth_topo);
                if plottopo==1
                    plot(topo{i},'Marker','*')
                    led=[led; {'after smoothing'}];
                    legend(led)
                    grid on
                end
            end
        else
            for i=1:length(global_coords)
                topo{i}=F(global_coords{i}(:,1),global_coords{i}(:,2));
                if plottopo==1
                    figure
                    plot(topo{i},'Marker','*')
                    hold on
                    xlabel('Trace number')
                    ylabel('Topography [m]')
                    led=[{'raw data'}];
                    legend(led)
                    grid on
                    title(['Profile ',int2str(i)])
                end
            end
        end
        % set topo in global_coords:
        for i=1:length(global_coords)
            global_coords{i}(:,3)=topo{i};
        end
        save(fullfile(pfad_rad,'global_coords.mat'),'global_coords');
    end
    % get min/max of topo
    for i=1:length(topo)
        mintop(i)=min(topo{i});
        maxtop(i)=max(topo{i});
    end
    mintopo=min(mintop);
    maxtopo=max(maxtop);

    % Displaying min/max depth values for migration
    disp(['Maximum z of topography is ',num2str(maxtopo,5),' m.'])
    disp(['Minimum z of topography is ',num2str(mintopo,5),' m.'])
    disp(['Maximum depth range of all radargrams is ',num2str(maxdepth,5),' m.'])
    disp(['  -> Recommended values for zmin and zmax in settings.txt are: zmax = ',int2str(ceil(maxtopo)),' m and zmin = ',int2str(floor(mintopo-maxdepth)),' m.'])


    save(fullfile(pfad_rad,'topo.mat'),'topo');
end


%% function for removing outliers in topography
function c=remout(c,anz,grenze)
% c: topography-vector
% anz: number of values for median calculation
% grenze: height difference over anz values in m

% moving median calculation
for i=(anz-1)/2+1:length(c)-(anz-1)/2
    med(i)=median(c(i-(anz-1)/2:i+(anz-1)/2));
end
med(1:(anz-1)/2)=med((anz-1)/2+1);
med(length(c)-(anz-1)/2+1:length(c))=med(length(c)-(anz-1)/2);

% if difference of median and original values > grenze -> delete and interpolate
weg=abs(med-c')>grenze;
x=1:length(c);
c=interp1(x(~weg),med(~weg),x);
end