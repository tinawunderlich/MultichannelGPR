clear all
close all
clc

% Joining of binned data (on same or similar grid) into big matrix (probably not possible
% for very large areas due to computer memory!)
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% Requires binned 3D data in folders 3D_Grid_R* or 3D_Grid_R*/processed


numbers=1:16;   % Numbers of datasets for joining
% (not the 3D_Grid_R*-number, but number of calls for choosing a folder!)

rectangles=1:12; % Number of rectangles (3D_Grid_R*) (has to be the same in all datasets)


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

for i=1:length(numbers) % for every folder to be chosen...

    % select folder
    % get folder name
    if ispc
        if exist('temp.temp') % read last opened folder from temp.temp
            fid=fopen('temp.temp','r');
            fn=textscan(fid,'%s');
            fclose(fid);
            if ~isempty(fn{1})
                pathname=uigetdir(fn{1}{1},['Choose ',int2str(i),'. folder with bins']);
            else
                pathname=uigetdir([],['Choose ',int2str(i),'. folder with bins']);
            end
            fileattrib('temp.temp','-h');
            fid=fopen('temp.temp','wt');
            fprintf(fid,'%s',pathname);
            fclose(fid);
            fileattrib('temp.temp','+h');
        else
            pathname=uigetdir([],['Choose ',int2str(i),'. folder with bins']); % path to radargram-folder

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
                pathname=uigetdir(fn{1}{1},['Choose ',int2str(i),'. folder with bins']);
            else
                pathname=uigetdir([],['Choose ',int2str(i),'. folder with bins']);
            end
        else
            pathname=uigetdir([],['Choose ',int2str(i),'. folder with bins']); % path to radargram-folder
        end

        fid=fopen('.temp.temp','wt');
        fprintf(fid,'%s',pathname);
        fclose(fid);
    end
    pn{i}=pathname;
end

for r=1:length(rectangles)
    disp(['RECTANGLE ',int2str(rectangles(r)),':'])

    for i=1:length(numbers)
        % read only x and y
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'x.mat'));
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'y.mat'));

        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'t.mat'));

        dx(i)=x(1,2)-x(1,1);   % grid interval in m
        dt(i)=abs(t(2)-t(1));   % time interval in ns
        ns(i)=length(t);   % number of samples

        xmin(i)=min(x(1,:));
        xmax(i)=max(x(1,:));
        ymin(i)=min(y(:,1));
        ymax(i)=max(y(:,1));

        % save grids
        xx{i}=x;
        yy{i}=y;

        disp(['Folder: ',pn{i},'/3D_Grid_R',int2str(rectangles(r))])
        disp(['dx = ',num2str(dx(i)),' m'])
        if i>1 && dx(i)~=dx(i-1)
            disp('Warning! Different dx!')
        end
        disp(['dt = ',num2str(dt(i)),' ns'])
        if i>1 && dt(i)~=dt(i-1)
            disp('Warning! Different dt!')
        end
        disp(['ns = ',int2str(ns(i))])
        if i>1 && ns(i)~=ns(i-1)
            disp('Warning! Different ns!')
        end

    end
    disp('------------------------------')
    disp('Reading data...')

    clear x y t;

    % make new t for all:
    t=0:min(dt):max(dt.*ns)-min(dt);

    %%% make new grid
    xrg=[min(xmin)-dx/2:dx:max(xmax)+dx/2]; % edges for binning
    yrg=[min(ymin)-dx/2:dx:max(ymax)+dx/2];

    zbin=[];
    dbin=[];
    xbin=[];
    ybin=[];
    for i=1:length(numbers)
        disp([int2str(i),' / ',int2str(max(numbers))])

        % load data and mask of current rectangle
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'data.mat'));
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'mask.mat'));
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'z.mat'));
        load(fullfile(pn{i},['3D_Grid_R',int2str(rectangles(r))],'coordtrans.mat'));

        % prepare vectors for binning:
        zbin=[zbin; z(mask==1)]; % topo values
        dbintemp=zeros(sum(sum(mask==1)),length(data(1,1,:))); % initialize matrix
        for j=1:length(data(1,1,:)) % for every time sample
            temp=data(:,:,j);
            dbintemp(:,j)=temp(mask==1); % amplitude values
        end
        xbin=[xbin; xx{i}(mask==1)];
        ybin=[ybin; yy{i}(mask==1)];

        % interpolate dbin-matrix for new t-vector:
        dbin=[dbin; interp1([0:dt(i):(ns(i)-1)*dt(i)],dbintemp.',t)'];
    end

    %% binning
    disp('Start binning...')
    % delete all values outside xrg/yrg
    weg=(xbin>max(xrg) | xbin<min(xrg) | ybin>max(yrg) | ybin<min(yrg));
    if any(weg)
        x(weg)=[];
        y(weg)=[];
        z(weg)=[];
    end

    % finding the x and y bins for each coordinate
    [~,whichbinx] = histc(xbin,xrg);   % returns the bin number in x direction (left edge is included in bin)
    [~,whichbiny] = histc(ybin,yrg);   % returns the bin number in y direction (left edge is included in bin)

    % corrected bin numbers (if there is a coordinate == last right bin edge it is set into last bin)
    binsx = min(max(whichbinx,1),length(xrg)-1);
    binsy = min(max(whichbiny,1),length(yrg)-1);

    bins = (binsy-1).*(length(xrg)-1)+binsx; % linear index of bin in matrix

    xpos = ones(size(bins,1),1);
    ns = sparse(bins,xpos,1,(length(xrg)-1)*(length(yrg)-1),1); % sparse matrix with number of counts in each bin
    % topo:
    disp('    - Topography')
    ysum = sparse(bins,xpos,zbin,(length(xrg)-1)*(length(yrg)-1),1);   % sparse matrix with sum of amplitude values in each bin
    zm = full(ysum)./(full(ns));    % vector with mean values (sum/number) for each bin in vector bins (if no hit, then NaN)
    z = reshape(zm,length(xrg)-1,length(yrg)-1)';  % output matrix -> topo
    % amplitudes:
    disp('    - Amplitudes')
    for i=1:length(dbin(1,:)) % for each sample
        if ~mod(i,100)
            disp(['       Sample ',int2str(i),' / ',int2str(length(t))])
        end
        dsum = sparse(bins,xpos,dbin(:,i),(length(xrg)-1)*(length(yrg)-1),1);   % sparse matrix with sum of amplitude values in each bin
        dmtemp = full(dsum)./(full(ns));    % vector with mean values (sum/number) for each bin in vector bins (if no hit, then NaN)
        dm(:,:,i) = reshape(dmtemp,length(xrg)-1,length(yrg)-1)';  % output matrix
    end

    [x,y]=meshgrid([xrg(1)+dx/2:dx:xrg(end)-dx/2],[yrg(1)+dx/2:dx:yrg(end)-dx/2]);

    mask_all(~isnan(dm(:,:,1)))=1;

    % renaming variables
    data=dm;
    mask=mask_all;

    disp('-----------------------')
    disp(['Saving data in new folder: ',pn{1},['/3D_Grid_R',int2str(rectangles(r))],'/AllBins'])


    if ~exist(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins'),'dir')
        mkdir(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins'))
    end
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','x.mat'),'x','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','y.mat'),'y','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','data.mat'),'data','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','mask.mat'),'mask','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','z.mat'),'z','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','t.mat'),'t','-v7.3');
    save(fullfile(pn{1},['3D_Grid_R',int2str(rectangles(r))],'AllBins','coordtrans.mat'),'coordtrans','-v7.3');
end