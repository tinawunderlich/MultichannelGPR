clear all
close all
clc

% Script for extracting radargrams from 3D-datablocks (Select rSlicer
% folder!)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires binned 3D data (processed or unprocessed),list of start/end coordinates of planned radargramms is optional, location af radargramms can be picked using this script



% which rectangles?
rect_start=1;
rect_end=15;

% which timeslice for profile selection?
tsl_start=20;   % start time in ns (or m if dsl=1)
tsl_end=22;     % ending time in ns (or m if dsl=1)
dsl=0; % =1 if vertical axis is depth, =0 if vertical axis is time
% Tipp: If you have depth migrated data, tsl_start and tsl_end are in m
% according to the depths in vector t. If you get an error message, please
% check these numbers. tsl_start<tsl_end!

% list of profile coordinates in rSlicer-folder?
% (columns are: xstart ystart xend yend)
list='Radargrams.txt'; % if no list, leave empty [] -> interactive picking in plot
%list =[];
coordglobal=1; % if the 3Dbins-R* are in local coordinates: do coordinate transformation if ==1

% trace spacing of 3D block:
dx=0.04;    % dx in m

dxrad=0.02; % trace spacing of new radargram in m

% use processed blocks in 3D_Grid_R*/processed? (yes==1)
proc=0;

% if you want to extract more profiles later, you can save the
% interpolation object and reload it later:
saveInterpolation=1;
reload=0; % set to 1 if you want to reload the interpolation handle in a second run of the program


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pathname=uigetdir([],'Choose rSlicer folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pathname);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pathname=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pathname);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pathname=uigetdir([],'Choose rSlicer folder');
        end
    else
        pathname=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pathname);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));


% load time vector (the same for all areas)
if proc==0
    load(fullfile(pathname,['3D_Grid_R',int2str(rect_start)],'t.mat'));
else
    load(fullfile(pathname,['3D_Grid_R',int2str(rect_start)],'processed','t.mat'));
end
t_ind=find(tsl_start<=t & tsl_end>=t);  % indices of samples for timeslice

if coordglobal==1
    load(fullfile(pathname,['3D_Grid_R',int2str(rect_start)],'coordtrans.mat'));
end

dt=t(2)-t(1);   % sampleinterval
ns=length(t);   % number of samples

for i=rect_start:rect_end % for each rectangle...
    % load coordinates
    if proc==0
        temp=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'x.mat'));
        temp2=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'y.mat'));
        temp3=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'z.mat'));
        matFileObj=matfile(fullfile(pathname,['3D_Grid_R',int2str(i)],'data.mat'));
    else
        temp=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'processed','x.mat'));
        temp2=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'processed','y.mat'));
        temp3=load(fullfile(pathname,['3D_Grid_R',int2str(i)],'processed','z.mat'));
        matFileObj=matfile(fullfile(pathname,['3D_Grid_R',int2str(i)],'processed','data.mat'));
    end
    xr{i}=temp.x;
    yr{i}=temp2.y;
    zr{i}=temp3.z;

    xmin(i)=min(xr{i}(1,:));
    xmax(i)=max(xr{i}(1,:));
    ymin(i)=min(yr{i}(:,1));
    ymax(i)=max(yr{i}(:,1));

    % load data only in selected timeslice
    tsl{i}=sum(abs(matFileObj.data(:,:,t_ind)),3);

    % for topomigrated/corrected data only: get bottom of data
    if dsl==1
        zrb{i}=NaN(size(xr{i}));
        data=permute(matFileObj.data(:,:,:),[3 1 2]);
        [r,c]=find(~isnan(data)); % r is index of time vector, c is linear index for x-y-grids
        b=unique(c); % find bins with non-isnan data (b is linear index corresponding to x-y-grids)
        for kk=1:length(b) % for each bin with data...
            zs=t(r(c==b(kk))); % get depth values with data
            zrb{i}(b(kk))=min(zs); % set minimum depth of data in large matrix
        end
    end

end

%%% make new grid
[xgrid,ygrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)]);  % including border!
alle=NaN(length(ygrid(:,1)),length(xgrid(1,:)));  % initialize big data matrix for timeslice
rectnumber=NaN(length(ygrid(:,1)),length(xgrid(1,:)));  % matrix for number of rectangle
rownum=NaN(length(ygrid(:,1)),length(xgrid(1,:))); % matrix with row number in each rectangle
colnum=NaN(length(ygrid(:,1)),length(xgrid(1,:))); % matrix with column number in each rectangle
if dsl==1
    topoall=NaN(length(ygrid(:,1)),length(xgrid(1,:))); % initialize big matrix for topography
    bottomall=NaN(length(ygrid(:,1)),length(xgrid(1,:))); % initialize big matrix for bottom topography beneath topomigrated/corrected data
end

for i=rect_start:rect_end
    % find indices of rectangle in large area (and adjust if necessary)
    startind(i,:)=[find(ygrid(:,1)>=ymin(i),1) find(xgrid(1,:)>=xmin(i),1)];
    endind(i,:)=[max(find(ygrid(:,1)<=ymax(i))) max(find(xgrid(1,:)<=xmax(i)))];
    if length(startind(i,1):endind(i,1))<=length(tsl{i}(:,1))   % if number of rows different
        anzrows=length(tsl{i}(:,1))-length(startind(i,1):endind(i,1));  % difference in rows
        if endind(i,1)+anzrows<=length(tsl{i}(:,1))
            endind(i,1)=endind(i,1)+anzrows;
        else
            startind(i,1)=startind(i,1)-anzrows;
        end
    end
    if length(startind(i,1):endind(i,1))>length(tsl{i}(:,1))   % if number of rows different
        endind(i,1)=startind(i,1)+length(tsl{i}(:,1))-1;
    end
    if length(startind(i,2):endind(i,2))<=length(tsl{i}(1,:))     % if number of columns different
        anzcol=length(tsl{i}(1,:))-length(startind(i,2):endind(i,2));
        if endind(i,2)+anzcol<=length(tsl{i}(1,:))
            endind(i,2)=endind(i,2)+anzcol;
        else
            startind(i,2)=startind(i,2)-anzcol;
        end
    end
    % set data in big matrix
    alle(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=tsl{i};
    if dsl==1
        topoall(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=zr{i};
        bottomall(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=zrb{i};
    end
    % save indices for finding of data in rectangles
    rectnumber(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=i;
    [c,r]=meshgrid([1:length(tsl{i}(1,:))],[1:length(tsl{i}(:,1))]);
    rownum(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=r;
    colnum(startind(i,1):endind(i,1),startind(i,2):endind(i,2))=c;
end

fig1=figure('Name','Timeslice','OuterPosition',[0 500 560 500]);
imagesc(xgrid(1,:),ygrid(:,1),alle)
axis xy
cmap=flipud(gray);
colormap(cmap)
xlabel('x [m]')
ylabel('y[m]')
crange=max(alle(:))-min(alle(:));
thresh=3;
set(gca,'FontSize',20,'Dataaspectratio',[1 1 1],'CLim',[median(alle(~isnan(alle(:))))-crange/100*15 median(alle(~isnan(alle(:))))+crange/100*30])

if isempty(list)
    ok=0;
    anz=1;
    while ok==0
        disp('Please pick the starting and ending point of the profile.');
        % coordinates will be in same system as timeslices
        [xpoints,ypoints]=ginput(2);
        xstart(anz,1)=xpoints(1);
        xend(anz,1)=xpoints(2);
        ystart(anz,1)=ypoints(1);
        yend(anz,1)=ypoints(2);

        hold on
        plot([xstart(anz) xend(anz)],[ystart(anz) yend(anz)],'r','Linewidth',2) % plot line
        drawnow;

        anz=anz+1;
        reply=input('Do you want to pick another profile? Yes=1, No=0 ','s');
        if strcmp(reply,'0')
            ok=1;
        end
    end
else
    disp('Reading list of coordinates.')
    coords=load(fullfile(pathname,list));
    xstart=coords(:,1);
    ystart=coords(:,2);
    xend=coords(:,3);
    yend=coords(:,4);
    if coordglobal==1
        % transform into local system fitting to timeslices
        for i=1:length(xstart)
            new=helmert([xstart(i) ystart(i); xend(i) yend(i)],coordtrans(:,3:4),coordtrans(:,1:2));
            xstart(i)=new(1,1);
            ystart(i)=new(1,2);
            xend(i)=new(2,1);
            yend(i)=new(2,2);
        end
    end

    hold on
    for i=1:length(xstart)
        plot([xstart(i) xend(i)],[ystart(i) yend(i)],'r','Linewidth',2) % plot line
    end
    drawnow;
end

disp('Please wait...')
for k=1:length(xstart)  % for every profile
    %%% interpolate onto selected profile
    dist=sqrt((xend(k)-xstart(k))^2+(ystart(k)-yend(k))^2); % length of profile
    xrad{k}=0:dxrad:dist; % x values along profile
    yy{k}=interp1([0 dist],[ystart(k) yend(k)],xrad{k}); % Koordinaten entlang Profil
    xx{k}=interp1([0 dist],[xstart(k) xend(k)],xrad{k});
end

mask_all=NaN(size(alle));
dat_t=NaN(size(alle));
disp('Creating all radargrams...')
for k=1:length(xstart)
    radargrams{k}=NaN(length(t),length(xrad{k}));
end
for i=1:length(t) % for every time step
    if ~mod(i,50)
        disp(['Time step ',int2str(i),'/',int2str(length(t))]);
    end
    % put together all rectangles for this time step:
    for j=rect_start:rect_end
        if proc==0
            matFileObj=matfile(fullfile(pathname,['3D_Grid_R',int2str(j)],'data.mat'));
        else
            matFileObj=matfile(fullfile(pathname,['3D_Grid_R',int2str(j)],'processed','data.mat'));
        end
%         if i==1 % create mask for all area
%             if proc==0
%                 load(fullfile(pathname,['3D_Grid_R',int2str(j)],'mask.mat'));
%             else
%                 load(fullfile(pathname,['3D_Grid_R',int2str(j)],'processed','mask.mat'));
%             end
%             mask_all(startind(j,1):endind(j,1),startind(j,2):endind(j,2))=mask;
%         end

        dat_t(startind(j,1):endind(j,1),startind(j,2):endind(j,2))=matFileObj.data(:,:,i); % complete area for this time step
    end

    % take data from this time step for every profile
    if i==1 && reload==1
        load(fullfile(pathname,'Interpolant.mat'));
    end
    if reload~=1
        F{i}=scatteredInterpolant(double(xgrid(~isnan(dat_t))),double(ygrid(~isnan(dat_t))),dat_t(~isnan(dat_t)));
    end
    for k=1:length(xstart)  % for every profile
        test=smooth(F{i}(xx{k},yy{k}),15);
        if ~isempty(test)
            radargrams{k}(i,:)=test;
        end
    end

end
if saveInterpolation==1
    save(fullfile(pathname,'Interpolant.mat'),'F','-v7.3');
end

% combine coordinates:
clear x;
for i=1:length(xx)
    global_coords{i}=[xx{i}' yy{i}' zeros(size(xx{i}'))];
    if coordglobal==1
        % transform into global system again for radargrams
        global_coords{i}(:,1:2)=helmert(global_coords{i}(:,1:2),coordtrans(:,1:2),coordtrans(:,3:4));
    end
    x{i}=xrad{i};
end


% if dsl==1: for each new profile delete pixels above topography and below
% bottom-of-data-topography:
if dsl==1
    Ft=scatteredInterpolant(xgrid(~isnan(topoall)),ygrid(~isnan(topoall)),topoall(~isnan(topoall)));
    Fb=scatteredInterpolant(xgrid(~isnan(bottomall)),ygrid(~isnan(bottomall)),bottomall(~isnan(bottomall)));
    for i=1:length(radargrams)
        topo=Ft(xx{i},yy{i});
        global_coords{i}(:,3)=topo; % set topo in global_coords
        bottom=Fb(xx{i},yy{i});
        % set values above topo and below bottom to NaN
        for k=1:length(radargrams{i}(1,:)) % for each trace
            radargrams{i}(t>=topo(k) | t<=bottom(k),k)=NaN;
        end
    end
end




% Save radargrams for further processing
if proc==0
    if ~exist(fullfile(pathname,'Radargrams'))
        mkdir(fullfile(pathname,'Radargrams'));
    end
    save(fullfile(pathname,'Radargrams','t.mat'),'t','-v7.3');
    save(fullfile(pathname,'Radargrams','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pathname,'Radargrams','x.mat'),'x','-v7.3');
    save(fullfile(pathname,'Radargrams','global_coords.mat'),'global_coords','-v7.3');

    % Save figure with profile location and radargram
    print(fullfile(pathname,'Radargrams',['Radargram_Location.png']),'-dpng','-r0');

else
    if ~exist(fullfile(pathname,'Radargrams_proc'))
        mkdir(fullfile(pathname,'Radargrams_proc'));
    end
    save(fullfile(pathname,'Radargrams_proc','t.mat'),'t','-v7.3');
    save(fullfile(pathname,'Radargrams_proc','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(pathname,'Radargrams_proc','x.mat'),'x','-v7.3');
    save(fullfile(pathname,'Radargrams_proc','global_coords.mat'),'global_coords','-v7.3');

    % Save figure with profile location and radargram
    print(fullfile(pathname,'Radargrams_proc',['Radargram_Locations.png']),'-dpng','-r0');
end


% set original path
path(oldpath);
