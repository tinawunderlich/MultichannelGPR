clear all
close all
clc


% Make Timeslices from Mala-GPR-data, preprocessed by Mala3D (Select
% rSlicer folder!)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires binned data in rectangles 3D_Grid_R* or processed data in
% 3D_Grid_R*/processed/
% Number of rectangles with binned data
rectangles=1; % e.g. 1:3

% depth slices instead of time slices (input is in m)
dsl = 0; % =1: depth, =0: time

% if depth slice, cut horizontally (=0) or follow topography (=1)?
followTopo=0;

% starting time of first timeslice
t_start=0;  % in ns (or m if depth slice starting from top of data (=0m))

% thickness of timeslices
thick=2; % in ns (or m if dsl=1)

% overlap of timeslices
overlap = 0; % in ns (or m if dsl=1)

% ending time of timeslices
t_end=30; % in ns (or m if dsl=1, meters below top of data, positive!)

% dx of timeslices (<=dx of bins)
dx_tsl=0.05;    % in m

% use 3D processed data (if =1, then use data in /processed folder in 3Dbins)
proc=0;

% normalize Timeslice to [0 1]
normTsl = 0;

% method for creation of timeslices
method=1;   % 1: sum absolute amplitudes,
            % 2: rms of absolute amplitudes
            % 3: no summing, just amplitudes at certain times
            % 4: no summing, just absolute amplitudes at certain times

% masking options
nn_radius=6;    % radius to nearest neighbor (in bins) should be less than nn_radius to be valid

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
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pfad=uigetdir([],'Choose rSlicer folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

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


% read time vector
if dsl
    if proc==0
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
    else
        if exist(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'processed','depth.mat'),'file')
            temp=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'processed','depth.mat'));
            t = temp.depth';
            maxElevation = temp.maxElevation;
        else
            t=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'processed','t.mat'));
            t=t.t';
            maxElevation=t(1); % absolute height at top
            t=abs(t-t(1)); % depth is starting from 0, positive down
        end    
    end
    fprintf("Depth vector boundaries in [m] are [%7.4f %7.4f] \n",min(t), max(t));
else
    if proc==0
        t=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
    else
        t=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'processed','t.mat'));
    end
    t=t.t';
    maxElevation=[];
    fprintf("Time vector boundaries in [ns] are [%7.4f %7.4f] \n",min(t), max(t));
end

% read coordinates
for i=1:length(rectangles)
    if proc==0
        p=['3D_Grid_R',int2str(rectangles(i))];
    else
        p=['3D_Grid_R',int2str(rectangles(i)),'/processed'];
    end
    temp=load(fullfile(pfad,p,'x.mat'));
    x{i}=temp.x;
    temp=load(fullfile(pfad,p,'y.mat'));
    y{i}=temp.y;
    temp=load(fullfile(pfad,p,'z.mat'));
    z{i}=temp.z;
    temp=load(fullfile(pfad,p,'mask.mat'));
    mask{i}=temp.mask;
    if isempty(mask{i})
        if proc==0
            p=['3D_Grid_R',int2str(rectangles(i)),'/data.mat'];
        else
            p=['3D_Grid_R',int2str(rectangles(i)),'/processed/data.mat'];
        end
        matFileObj=matfile(fullfile(pfad,p));
        temp=matFileObj.data;
        mask{i}=zeros(size(x{i}));
        mask{i}(any(~isnan(temp),3))=1;
    end
    
    if i==1
        dx=x{i}(1,2)-x{i}(1,1);   % dx of 3D-block
    end
    
    xmin(i)=min(x{i}(1,:));
    xmax(i)=max(x{i}(1,:));
    ymin(i)=min(y{i}(:,1));
    ymax(i)=max(y{i}(:,1));
end

% make grid for Timeslices
[xgrid,ygrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)]);

% make grid for interpolated Timeslices
[xgrid_interp,ygrid_interp]=meshgrid([min(xmin):dx_tsl:max(xmax)],[min(ymin):dx_tsl:max(ymax)]);


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
    if method==1 || method==2
        t_ind{i}=find(t>=start & t<ende); % find indices of time interval
        
        t_startende(i,:)=[start ende];
    elseif method==3 || method==4
        tempt=find(start-dt/2<=t & start+dt/2>=t);  % find index of starting time
        t_ind{i}=tempt(1);
        
        t_startende(i,:)=[start start];
    end
    
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

tsl=cell(length(t_ind),1);  % initialize cells for timeslices

disp('Reading data for timeslices... Please wait.');

for j=1:length(rectangles)  % in each rectangle...
    disp(['   Rectangle ',int2str(j),'/',int2str(length(rectangles))])
    % open matfile
    if proc==0
        p=['3D_Grid_R',int2str(rectangles(j)),'/data.mat'];
    else
        p=['3D_Grid_R',int2str(rectangles(j)),'/processed/data.mat'];
    end
    matFileObj=matfile(fullfile(pfad,p));
    
    % determine location in big area (find lower left corner)
    row_ind=find(abs(ygrid(:,1)-ymin(j))==min(abs(ymin(j)-ygrid(:,1))));
    col_ind=find(abs(xgrid(1,:)-xmin(j))==min(abs(xmin(j)-xgrid(1,:))));
    
    for i=1:length(t_ind)   % number of Tsl
        disp(['      - Tsl ',int2str(i),'/',int2str(length(t_ind))])
        if isempty(tsl{i})
            tsl{i}=NaN(size(xgrid));    % initialize tsl in first run
        end
        
        if followTopo==0 % horizontal slices
            datatemp=matFileObj.data(:,:,t_ind{i});
        else % make Tsl following topography (for topocorr/mig data)
            z_absolut=maxElevation-t; % absolute z vector starting from 0 at top of radargram images
            datatemp=NaN(length(z{j}(:,1)),length(z{j}(1,:)),length(t_ind{i})); % initialize datatemp
            data=matFileObj.data; % load 3D cube
            % find rows in area with any ~isnan-value:
            kk=find(any(~isnan(z{j}),2)); % indices of rows with values
            % get data-cube-slice along columns:
            for k=1:length(kk)
                temp=permute(data(kk(k),:,:),[3 2 1]); % vertical data-slice
                ztemp=repmat(z{j}(kk(k),:),[length(t),1]); % topography along this slice (as matrix for each t sample)
                b=zeros(size(ztemp));
                for s=1:length(t_ind{i}) % go through depthsamples for this slice
                     b=b+double(min(abs(ztemp-(t_ind{i}(s)-1)*dt-z_absolut))==abs(ztemp-(t_ind{i}(s)-1)*dt-z_absolut)); % matrix where all samples in this slice are 1 and rest is 0
                end
                g=reshape(temp(logical(b)),length(t_ind{i}),[]); % all samples in this slices in one matrix (without nan-traces in between)
                % sometimes there are still complete Nan-traces -> delete them
                ind=find(~isnan(ztemp(1,:))); % trace-indices in data-cube slice that are not nan
                ind(find(all(isnan(g),1)))=[]; % sometimes there are still complete Nan-traces -> delete them in indices 
                g(:,all(isnan(g),1))=[]; % sometimes there are still complete Nan-traces -> delete them in trace matrix
                % sometimes there are some nan-values on top -> extrapolate
                w=find(isnan(g(1,:))); % columns with nan as first value
                for l=1:length(w)
                    firstval=g(find(~isnan(g(:,w(l))),1,'first'),w(l)); % first not-nan value of this trace
                    if length(find(~isnan(g(:,w(l)))))>1
                        g(:,w(l))=interp1(find(~isnan(g(:,w(l)))),g(~isnan(g(:,w(l))),w(l)),1:length(t_ind{i}),'linear',firstval);
                    else
                        g(:,w(l))=NaN(size(g(:,w(l))));
                    end
                end
                % sometimes there are some nan-values at the bottom -> extrapolate
                w=find(isnan(g(end,:))); % columns with nan as last value
                for l=1:length(w)
                    lastval=g(find(~isnan(g(:,w(l))),1,'last'),w(l)); % last non-nan value of this trace
                    if length(find(~isnan(g(:,w(l)))))>1
                        g(:,w(l))=interp1(find(~isnan(g(:,w(l)))),g(~isnan(g(:,w(l))),w(l)),1:length(t_ind{i}),'linear',lastval);
                    else
                        g(:,w(l))=NaN(size(g(:,w(l))));
                    end
                end
                % prepare new slice
                new=NaN(length(t_ind{i}),size(temp,2)); 
                new(:,ind)=g; % new dataslice in this slice (with missing traces in between)
                % set this into datatemp matrix:
                datatemp(kk(k),:,:)=permute(new,[3 2 1]);
            end
            len_ind=length(t_ind{i});
        end
        if method==1  % sum absolute amplitudes
            tsl_temp=sum(abs(datatemp),3);
        elseif method==2 % rms of absolute amplitudes
            tsl_temp=sqrt(sum(abs(datatemp).^2,3)./len_ind);
        elseif method==3
            tsl_temp=datatemp;
        elseif method==4
            tsl_temp=abs(datatemp);
        end
        
        if i==1
            if j==1
                topo=NaN(size(xgrid));
            end
            topo(row_ind:row_ind+length(z{j}(:,1))-1,col_ind:col_ind+length(z{j}(1,:))-1)=z{j};
            topo=topo(1:length(xgrid(:,1)),1:length(xgrid(1,:)));   % resize if necessary
        end
        
        % put tsl_temp together in big tsl
        tsl{i}(row_ind:row_ind+length(tsl_temp(:,1))-1,col_ind:col_ind+length(tsl_temp(1,:))-1)=tsl_temp;
        tsl{i}=tsl{i}(1:length(xgrid(:,1)),1:length(xgrid(1,:)));   % resize if necessary
    end
    
    clear tsl_temp;
end


% Make mask
if dsl==0 % TIMEslice -> just one mask for all slices
    mask=zeros(size(tsl{1}));
    for i=1:length(tsl)
        mask(~isnan(tsl{i}))=1;
    end
    mask(mask>=1)=1;
    topo(mask==0)=NaN;
    mask_topo=mask;
else % DEPTHslices -> one mask for each slice
    for i=1:length(tsl)
        mask{i}=zeros(size(tsl{i}));
        mask{i}(~isnan(tsl{i}))=1;
    end
    
    % one for topo:
    mask_topo=zeros(size(tsl{1}));
    for i=1:length(tsl)
        mask_topo(~isnan(tsl{i}))=1;
    end
    mask_topo(mask_topo>=1)=1;
    topo(mask_topo==0)=NaN;
end


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
        timerest=diff*length(tsl)-diff*i;
        if timerest>2
            disp(['    Approx. ',int2str(timerest),' minutes remaining'])
        else
            disp(['    Approx. ',int2str(timerest*60),' seconds remaining'])
        end
    end
    topo_interp=griddata(double(xgrid(~isnan(topo(:)))),double(ygrid(~isnan(topo(:)))),topo(~isnan(topo(:))),double(xgrid_interp),double(ygrid_interp));
    topo_interp(~mask_topointerp)=NaN;   % apply mask_interp

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

if normTsl
    for i=1:length(tsl)
        normTsl = tsl{i} - min(tsl{i}(:));
        tsl{i} = normTsl ./ max(normTsl(:));

        normTsl_interp = tsl_interp{i} - min(tsl_interp{i}(:));
        tsl_interp{i} = normTsl_interp ./ max(normTsl_interp(:));
    end
end

%%% Create folder and save data
disp('Saving data and info files.')
if proc==0
    mkdir(fullfile(pfad,'Timeslices'));
    mkdir(fullfile(pfad,'Timeslices','interpolated'));
    % save original timeslices
    save(fullfile(pfad,'Timeslices','tsl.mat'),'tsl','-v7.3');
    save(fullfile(pfad,'Timeslices','xgrid.mat'),'xgrid','-v7.3');
    save(fullfile(pfad,'Timeslices','ygrid.mat'),'ygrid','-v7.3');
    save(fullfile(pfad,'Timeslices','topo.mat'),'topo','-v7.3');
    save(fullfile(pfad,'Timeslices','mask.mat'),'mask','-v7.3');
    save(fullfile(pfad,'Timeslices','t_startende.mat'),'t_startende','-v7.3');
    if dsl
        depth=t;
        save(fullfile(pfad,'Timeslices','depth.mat'),'depth','maxElevation','followTopo','-v7.3'); % aldo save depth vector/maxElevation and followTopo-flag
    end
    % save interpolated timeslices
    save(fullfile(pfad,'Timeslices','interpolated','tsl_interp.mat'),'tsl_interp','-v7.3');
    save(fullfile(pfad,'Timeslices','interpolated','xgrid_interp.mat'),'xgrid_interp','-v7.3');
    save(fullfile(pfad,'Timeslices','interpolated','ygrid_interp.mat'),'ygrid_interp','-v7.3');
    save(fullfile(pfad,'Timeslices','interpolated','topo_interp.mat'),'topo_interp','-v7.3');
    save(fullfile(pfad,'Timeslices','interpolated','mask_interp.mat'),'mask_interp','-v7.3');
    save(fullfile(pfad,'Timeslices','interpolated','t_startende.mat'),'t_startende','-v7.3');
    if dsl
        depth=t;
        save(fullfile(pfad,'Timeslices','interpolated','depth.mat'),'depth','maxElevation','followTopo','-v7.3');
    end
    
    % write info-files
    fid=fopen(fullfile(pfad,'Timeslices','tslinfo.txt'),'wt');
    fprintf(fid,'Numbers of used rectangles: %d\n',rectangles);
    if dsl
        fprintf(fid,'Thickness of timeslices: %4.2f m\n',thick);
    else
        fprintf(fid,'Thickness of timeslices: %4.2f ns\n',thick);
    end
    if method==1
        fprintf(fid,'Method 1: sum of absolute amplitudes\n');
    elseif method==2
        fprintf(fid,'Method 2: rms of absolute amplitudes\n');
    elseif method==3
        fprintf(fid,'Method 3: amplitudes\n');
    elseif method==4
        fprintf(fid,'Method 4: absolute amplitudes\n');
    end
    fprintf(fid,'Grid increment dx: %4.2f m\n',dx);
    if dsl
        fprintf(fid,'Maximum Elevation (=0 m depth) is: %4.2f m\n',maxElevation);
        fprintf(fid,'followTopo = %d\n',followTopo);
    end
    fclose(fid);
    
    fid=fopen(fullfile(pfad,'Timeslices','interpolated','tslinfo.txt'),'wt');
    fprintf(fid,'Numbers of used rectangles: %d\n',rectangles);
    if dsl
        fprintf(fid,'Thickness of timeslices: %4.2f m\n',thick);
    else
        fprintf(fid,'Thickness of timeslices: %4.2f ns\n',thick);
    end
    if method==1
        fprintf(fid,'Method 1: sum of absolute amplitudes\n');
    elseif method==2
        fprintf(fid,'Method 2: rms of absolute amplitudes\n');
    elseif method==3
        fprintf(fid,'Method 3: amplitudes\n');
    elseif method==4
        fprintf(fid,'Method 4: absolute amplitudes\n');
    end
    if dsl
        fprintf(fid,'Maximum Elevation (=0 m depth) is: %4.2f m\n',maxElevation);
        fprintf(fid,'followTopo = %d\n',followTopo);
    end
    fprintf(fid,'Grid increment dx of interpolated timeslices: %4.2f m\n',dx_tsl);
    if griding==1
        fprintf(fid,'Interpolated using griddata\n');
    elseif griding==2
        fprintf(fid,'Interpolated using Inverse Distance Weighting (IDW) with radius = %dm and power = %d\n',[radius power]);
    end
    fclose(fid);
else
    mkdir(fullfile(pfad,'Timeslices_proc'));
    mkdir(fullfile(pfad,'Timeslices_proc','interpolated'));
    % save original timeslices
    save(fullfile(pfad,'Timeslices_proc','tsl.mat'),'tsl','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','xgrid.mat'),'xgrid','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','ygrid.mat'),'ygrid','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','topo.mat'),'topo','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','mask.mat'),'mask','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','t_startende.mat'),'t_startende','-v7.3');
    if dsl
        depth=t;
        save(fullfile(pfad,'Timeslices_proc','depth.mat'),'depth','maxElevation','followTopo','-v7.3');
    end
    % save interpolated timeslices
    save(fullfile(pfad,'Timeslices_proc','interpolated','tsl_interp.mat'),'tsl_interp','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','interpolated','xgrid_interp.mat'),'xgrid_interp','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','interpolated','ygrid_interp.mat'),'ygrid_interp','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','interpolated','topo_interp.mat'),'topo_interp','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','interpolated','mask_interp.mat'),'mask_interp','-v7.3');
    save(fullfile(pfad,'Timeslices_proc','interpolated','t_startende.mat'),'t_startende','-v7.3');
    if dsl
        depth=t;
        save(fullfile(pfad,'Timeslices_proc','interpolated','depth.mat'),'depth','maxElevation','followTopo','-v7.3');
    end
    
    % write info-files
    fid=fopen(fullfile(pfad,'Timeslices_proc','tslinfo.txt'),'wt');
    fprintf(fid,'Numbers of used rectangles: %d\n',rectangles);
    if dsl
        fprintf(fid,'Thickness of timeslices: %4.2f m\n',thick);
    else
        fprintf(fid,'Thickness of timeslices: %4.2f ns\n',thick);
    end
    if method==1
        fprintf(fid,'Method 1: sum of absolute amplitudes\n');
    elseif method==2
        fprintf(fid,'Method 2: rms of absolute amplitudes\n');
    elseif method==3
        fprintf(fid,'Method 3: amplitudes\n');
    elseif method==4
        fprintf(fid,'Method 4: absolute amplitudes\n');
    end
    fprintf(fid,'Grid increment dx: %4.2f m\n',dx);
    if dsl
        fprintf(fid,'Maximum Elevation (=0 m depth) is: %4.2f m\n',maxElevation);
        fprintf(fid,'followTopo = %d\n',followTopo);
    end
    fclose(fid);
    
    fid=fopen(fullfile(pfad,'Timeslices_proc','interpolated','tslinfo.txt'),'wt');
    fprintf(fid,'Numbers of used rectangles: %d\n',rectangles);
    if dsl
        fprintf(fid,'Thickness of timeslices: %4.2f m\n',thick);
    else
        fprintf(fid,'Thickness of timeslices: %4.2f ns\n',thick);
    end

    if method==1
        fprintf(fid,'Method 1: sum of absolute amplitudes\n');
    elseif method==2
        fprintf(fid,'Method 2: rms of absolute amplitudes\n');
    elseif method==3
        fprintf(fid,'Method 3: amplitudes\n');
    elseif method==4
        fprintf(fid,'Method 4: absolute amplitudes\n');
    end
    if dsl
        fprintf(fid,'Maximum Elevation (=0 m depth) is: %4.2f m\n',maxElevation);
        fprintf(fid,'followTopo = %d\n',followTopo);
    end
    fprintf(fid,'Grid increment dx of interpolated timeslices: %4.2f m\n',dx_tsl);
    if griding==1
        fprintf(fid,'Interpolated using griddata\n');
    elseif griding==2
        fprintf(fid,'Interpolated using Inverse Distance Weighting (IDW) with radius = %dm and power = %d\n',[radius power]);
    end
    fclose(fid);
end


% Apply mask_interp on interpolated tsl for plotting
for i=1:length(tsl)
    if dsl==0
        maske=mask_interp;
    else
        maske=mask_interp{i};
    end
    tsl_interp{i}=tsl_interp{i}.*maske;
end

% Plot timeslices
disp('Plot timeslices');
% plot timeslices with limited GUI options
if proc==0
    if exist(fullfile(pfad,['3D_Grid_R',int2str(j)],'coordtrans.mat'),'file')
        load(fullfile(pfad,['3D_Grid_R',int2str(j)],'coordtrans.mat'));
        % save in Tsl folder:
        copyfile(fullfile(pfad,['3D_Grid_R',int2str(j)],'coordtrans.mat'),fullfile(pfad,'Timeslices','coordtrans.mat'));
        copyfile(fullfile(pfad,['3D_Grid_R',int2str(j)],'coordtrans.mat'),fullfile(pfad,'Timeslices','interpolated','coordtrans.mat'));
        % plot:
        Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices','interpolated'),dsl,maxElevation,followTopo,coordtrans);
    else
        Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices','interpolated'),dsl,maxElevation,followTopo);
    end
else
    if exist(fullfile(pfad,['3D_Grid_R',int2str(j)],'processed','coordtrans.mat'),'file')
        load(fullfile(pfad,['3D_Grid_R',int2str(j)],'processed','coordtrans.mat'));
        % save in Tsl folder:
        copyfile(fullfile(pfad,['3D_Grid_R',int2str(j)],'processed','coordtrans.mat'),fullfile(pfad,'Timeslices_proc','coordtrans.mat'));
        copyfile(fullfile(pfad,['3D_Grid_R',int2str(j)],'processed','coordtrans.mat'),fullfile(pfad,'Timeslices_proc','interpolated','coordtrans.mat'));
        % plot:
        Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices_mig','interpolated'),dsl,maxElevation,followTopo,coordtrans);
    else
        Tsl_slider_plot(xgrid_interp,ygrid_interp,tsl_interp,topo_interp,t_startende,fullfile(pfad,'Timeslices_mig','interpolated'),dsl,maxElevation,followTopo);
    end
end

waitfor(gcf);

% set original path
path(oldpath);
