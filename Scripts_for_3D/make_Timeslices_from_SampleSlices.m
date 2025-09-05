clear all
close all
clc


% Use Sampleslices created by PondView or SampleSlices_MalaProc to create
% Time/Depthslices
%
% Dr. Tina Wunderlich, CAU Kiel 2025, tina.wunderlich@ifg.uni-kiel.de
%


% depth slices instead of time slices (input is in m)
dsl = 0; % =1: depth, =0: time

% if depth slice, cut horizontally (=0) or follow topography (=1)?
followTopo=0;

% starting time of first timeslice
t_start=0;  % in ns (or m if depth slice starting from top of data (=0m))

% thickness of timeslices
thick=5; % in ns (or m if dsl=1)

% overlap of timeslices
overlap = 0; % in ns (or m if dsl=1)

% ending time of timeslices
t_end=50; % in ns (or m if dsl=1, meters below top of data, positive!)

% normalize Timeslice to [0 1]
normTsl = 0; % if =1: yes

% method for creation of timeslices
method=3;   % 1: sum absolute amplitudes [=sum(abs(A))]
            % 2: rms of absolute amplitudes [=sqrt(sum(A^2))]
            % 3: sqrt of sum of abs amplitudes [=sqrt(sum(abs(A)))]


%% -------------------------------------------------------------------------
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
            pfad=uigetdir(fn{1}{1},'Choose folder containing sampleslices');
        else
            pfad=uigetdir([],'Choose folder containing sampleslices');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder containing sampleslices'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder containing sampleslices');
        else
            pfad=uigetdir([],'Choose folder containing sampleslices');
        end
    else
        pfad=uigetdir([],'Choose folder containing sampleslices'); % path to radargram-folder
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


% read time vector
if dsl
    if exist(fullfile(pfad,'maxElevation.mat'),'file')
        load(fullfile(pfad,'maxElevation.mat'));
        load(fullfile(pfad,'t.mat'));
    else % depth.mat does not exist, but depth is in t.mat (absolute depths!)
        load(fullfile(pfad,'t.mat'));
        maxElevation=t(1); % absolute height at top
        t=abs(t-t(1)); % depth is starting from 0, positive down
    end
    fprintf("Depth vector boundaries in [m] are [%7.4f %7.4f] \n",min(t), max(t));
else
    load(fullfile(pfad,'t.mat'));
    maxElevation=[];
    fprintf("Time vector boundaries in [ns] are [%7.4f %7.4f] \n",min(t), max(t));
end

% read coordinates
load(fullfile(pfad,'coordtrans.mat'));
load(fullfile(pfad,'xgrid.mat'));
load(fullfile(pfad,'ygrid.mat'));

dx=abs(xgrid(1,1)-xgrid(1,2));

% read mask
load(fullfile(pfad,'mask_interp.mat'));

% read topo
load(fullfile(pfad,'topo_interp.mat'));

% make grid for interpolated Timeslices
xgrid_interp=xgrid;
ygrid_interp=ygrid;

% calculate time indices for timeslices
start=t_start;  % start time of first tsl
ende=t_start+thick; % ending time of first tsl
dt=abs(t(2)-t(1)); % sample interval

% t_end needs to be smaller
if t_end > t(end)
    t_end = t(end);
end
if overlap * 2 > thick
    fprintf("overlap = %5.2f is greater than the half of thick (= %5.2f)\noverlap is changed to %5.2f\n",overlap,thick,thick/2);
    overlap = thick/2;
end

nTsl = abs(ceil((t_end-t_start)/(thick - overlap)));
t_startende=zeros(nTsl,2);

for i=1:nTsl
    t_ind{i}=find(t>=start & t<ende); % find indices of time interval
    t_startende(i,:)=[start ende];

    start=ende - overlap;
    ende=start+thick; % Tsl

    if ende > t_end
        break
    end
end

% remove empty entries from t_startende
if i < nTsl
    t_startende(i+1:end,:) = [];
end

tsl_interp=cell(length(t_ind),1);  % initialize cells for timeslices

disp('Reading data for timeslices...');

numbers=1:length(t);
for i=1:length(tsl_interp) % for each timeslice
    disp(['  Tsl ',int2str(i),' / ',int2str(length(tsl_interp))])

    % initialize tsl
    tsl_interp{i}=zeros(size(xgrid));

    if followTopo==0
        % load sampleslices
        nums=numbers(t>=t_startende(i,1) & t<t_startende(i,2)); % numbers of sample slices
        for j=1:length(nums)
            load(fullfile(pfad,['slice_',int2str(nums(j))]));
            if method==1  % sum absolute amplitudes
                tsl_interp{i}=tsl_interp{i}+abs(slice);
            elseif method==2 % rms of absolute amplitudes
                tsl_interp{i}=tsl_interp{i}+abs(slice).^2;
            elseif method==3 % sqrt(sum(abs(A)))
                tsl_interp{i}=tsl_interp{i}+abs(slice);
            end
        end
    else
        
    end

    if method==2
        tsl_interp{i}=sqrt(tsl_interp{i}/length(nums));
    elseif method==3
        tsl_interp{i}=sqrt(tsl_interp{i});
    end

    if normTsl==1
        tsl_interp{i}=(tsl_interp{i}-mean(tsl_interp{i}(:)))/std(tsl_interp{i}(:));
    end
end


%%% Create folder and save data
disp('Saving data and info files...')

mkdir(fullfile(pfad,'Timeslices_interpolated'));
% save original timeslices
save(fullfile(pfad,'Timeslices_interpolated','tsl_interp.mat'),'tsl_interp','-v7.3');
save(fullfile(pfad,'Timeslices_interpolated','xgrid_interp.mat'),'xgrid_interp','-v7.3');
save(fullfile(pfad,'Timeslices_interpolated','ygrid_interp.mat'),'ygrid_interp','-v7.3');
save(fullfile(pfad,'Timeslices_interpolated','topo_interp.mat'),'topo_interp','-v7.3');
save(fullfile(pfad,'Timeslices_interpolated','mask_interp.mat'),'mask_interp','-v7.3');
save(fullfile(pfad,'Timeslices_interpolated','t_startende.mat'),'t_startende','-v7.3');
if dsl
    depth=t;
    save(fullfile(pfad,'Timeslices_interpolated','depth.mat'),'depth','maxElevation','followTopo','-v7.3'); % aldo save depth vector/maxElevation and followTopo-flag
end

% write info-files
fid=fopen(fullfile(pfad,'Timeslices_interpolated','tslinfo.txt'),'wt');
fprintf(fid,'Created timeslices from sample slices.\n');
if dsl
    fprintf(fid,'Thickness of depthslices: %4.2f m\n',thick);
else
    fprintf(fid,'Thickness of timeslices: %4.2f ns\n',thick);
end
if method==1
    fprintf(fid,'Method 1: sum of absolute amplitudes\n');
elseif method==2
    fprintf(fid,'Method 2: rms of absolute amplitudes\n');
elseif method==3
    fprintf(fid,'Method 3: sqrt of sum of absolute amplitudes\n');
end
fprintf(fid,'Grid increment dx: %4.2f m\n',dx);
if dsl
    fprintf(fid,'Maximum Elevation (=0 m depth) is: %4.2f m\n',maxElevation);
    fprintf(fid,'followTopo = %d\n',followTopo);
end
fclose(fid);


% Apply mask_interp on interpolated tsl for plotting
for i=1:length(tsl_interp)
    if dsl==0 % time
        maske=mask_interp;
    else
        maske=mask_interp{i};
    end
    tsl_interp{i}=tsl_interp{i}.*maske;
end

% Plot timeslices
disp('Plot timeslices');
% plot timeslices with limited GUI options
if exist(fullfile(pfad,'coordtrans.mat'),'file')
    load(fullfile(pfad,'coordtrans.mat'));
    % save in Tsl folder:
    copyfile(fullfile(pfad,'coordtrans.mat'),fullfile(pfad,'Timeslices_interpolated','coordtrans.mat'));
    % plot:
    Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices_interpolated'),dsl,maxElevation,followTopo,coordtrans);
else
    coordtrans=[1 1 1 1; 2 2 2 2];
    Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices_interpolated'),dsl,maxElevation,followTopo,coordtrans);
end


waitfor(gcf);

% set original path
path(oldpath);
