clear all
close all
clc

%%% Read sampleslices and create a video of them or export as geotif/geopng
%
% created with the help of claude.ai (free version), carefully checked and
% edited by Tina Wunderlich, CAU Kiel, 2026
%
% When asked, give the folder with all sampleslices.

%%% Video settings:
createVideo=1; % yes=1, no=0
videoname='BorgsumburgSuedtor'; % name of file, will be appended with information on the creation of slices
fps = 1;  % frames per second (speed of video)

%%% geopng/geotif settings:
save_geopng=1; % 1=yes
save_geotif=1; % 1=yes
epsg=25832; % epsg code of CRS

%%% other settings:
colperc=3; % Colorscale clipping in percent (if =0: autoscale min-max)
removeBorder=1; % =1: remove border artifacts from interpolation, =0: leave as it is
pix=6; % if removeBorder==1: how many pixels are removed from border around area
medianFilter=1; % do you want to apply a 2D-median filter (1=yes, 0=no)
msize=3; % filter size in pixel
% use squareroot of amplitudes for visualization?
sq=0; % 1=yes, 0=no






%% ---- DO NOT CHANGE BELOW THIS LINE ------------------------
% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));



% get folder name
if ispc
    if exist('temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with sampleslices');
        else
            foldername=uigetdir([],'Choose folder with sampleslices');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Choose folder with sampleslices'); % path to sampleslices-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    end
else
    if exist('.temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with sampleslices');
        else
            foldername=uigetdir([],'Choose folder with sampleslices');
        end
    else
        foldername=uigetdir([],'Choose folder with sampleslices'); % path to sampleslices-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% read data
load(fullfile(foldername,'xgrid.mat'));
load(fullfile(foldername,'ygrid.mat'));
load(fullfile(foldername,'t.mat'));
load(fullfile(foldername,'mask_interp.mat'));
load(fullfile(foldername,'mask_interp_topo.mat'));
load(fullfile(foldername,'topo_interp.mat'));

if exist(fullfile(foldername,'coordtrans.mat'))
    load(fullfile(foldername,'coordtrans.mat'));
else
    coordtrans=[1 1 1 1; 2 2 2 2];
end

% topo:
if removeBorder==1 % remove interpolation artifacts around area
    disttopo=chamfer_DT(mask_interp_topo);
    topo_interp(disttopo<=pix)=NaN;
end
if medianFilter==1
    topo_interp=medianfilt2(topo_interp,[msize msize]);
end


% create datastore with sample slices:
files = dir(fullfile(foldername, 'slice_*.mat'));
% sort numerically:
nums = arrayfun(@(f) sscanf(f.name, 'slice_%d.mat'), files,'UniformOutput',false);
files(cellfun(@(x) isempty(x),nums))=[];
nums(cellfun(@(x) isempty(x),nums))=[];
nums=cell2mat(nums);
[~, sortIdx] = sort(nums);
files = files(sortIdx);
% with sorted list: 1, 2, 3, .... (instead of 1, 10, 11, ..., 19, 2, 20, ...)
filePaths = fullfile(foldername, {files.name});
ds = fileDatastore(filePaths, "ReadFcn", @load);
numFrames=numel(ds.Files); % number of slices


i=1;
disp('Loading data...')
while hasdata(ds)
    data=read(ds);

    fprintf('   %d/%d\n',i,numFrames);

    if i==1
        slices=NaN(size(data.slice,1),size(data.slice,2),numFrames);
    end
    % apply mask
    if isscalar(mask_interp)
        mask=mask_interp{1};
    else
        mask=mask_interp{i};
    end
    temp=data.slice;

    if removeBorder==1
        dist{i} = chamfer_DT(mask);
        temp(dist{i}<=pix)=NaN;
    else
        dist{i}=[];
    end

    if medianFilter==1
        temp=medianfilt2(temp,[msize msize]);
    end

    if sq==1
        temp=sqrt(abs(temp));
    end

    slices(:,:,i)=temp;

    i=i+1;
end


% read info file:
info = readConfig(fullfile(foldername,'configuration.txt'));


%% VIDEO:
if createVideo==1
    disp('-- Creating video --')
    % create video name:
    if sq==0
        name=fullfile(foldername,[videoname,'_tz_flag',int2str(info.tz_flag),'_followTopo',int2str(info.followTopo),'.mp4']);
    else
        name=fullfile(foldername,[videoname,'_tz_flag',int2str(info.tz_flag),'_followTopo',int2str(info.followTopo),'_sqrt.mp4']);
    end

    % colorscale:
    if colperc==0
        cmin  = min(slices(~isnan(slices)));
        cmax = max(slices(~isnan(slices)));
    else
        coldata=sort(slices(~isnan(slices)));
        if ~isempty(coldata) && length(coldata)>2
            cmin=coldata(round(length(coldata)/100*colperc));
            cmax=coldata(end-round(length(coldata)/100*colperc));
        end
    end

    v = VideoWriter(name, 'MPEG-4');
    v.FrameRate = fps;
    open(v);

    fig = figure('Color', 'w');
    for k = 1:numFrames
        imagesc(xgrid(1,:),ygrid(:,1),slices(:,:,k))
        set(gca,'CLim',[cmin cmax])
        colormap(flipud(gray));
        axis xy
        axis equal

        if info.tz_flag==1 && info.followTopo==0
            title(sprintf('t = %.2f ns', t(k)));
        else
            title(sprintf('z = %.2f m', t(k)));
        end

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame.cdata);
    end

    close(v);
    close(fig);
end


%% GEOPNG
if save_geopng==1
    disp('-- Exporting as geopng/pgw --')

    if info.tz_flag==1 && info.followTopo==0
        saveallslices_geopng(xgrid,ygrid,slices,topo_interp,t,foldername,colperc,coordtrans,1,epsg,sq);
    else
        saveallslices_geopng(xgrid,ygrid,slices,topo_interp,t,foldername,colperc,coordtrans,2,epsg,sq);
    end
end


%% GEOTIF
if save_geotif==1
    disp('-- Exporting as geotif --')
    if ~exist(fullfile(foldername,'georef_tif'),'dir')
        mkdir(fullfile(foldername,'georef_tif'));
    end

    % get UTM coords of corners auf area:
    corners=helmert([xgrid(end,1) ygrid(end,1); xgrid(end,end) ygrid(end,end); xgrid(1,1) ygrid(1,1); xgrid(1,end) ygrid(1,end)],coordtrans(:,1:2),coordtrans(:,3:4));
    
    % topo:
    maxZ=ceil(max(topo_interp(:))*10)/10; % max topo
    minZ=floor(min(topo_interp(:))*10)/10; % min topo
    writeGeoTIFF(fullfile(foldername,'georef_tif','Topo_interp.tif'), flipud(topo_interp), epsg, corners, 'jet', [minZ maxZ]);
    % sampleslices:
    if info.tz_flag==1 && info.followTopo==0
        saveallslices_geotif(corners,epsg,slices,t,fullfile(foldername,'georef_tif'),colperc,1,sq);
    else
        saveallslices_geotif(corners,epsg,slices,t,fullfile(foldername,'georef_tif'),colperc,2,sq);
    end
end


% set original path
path(oldpath);


%%
function cfg = readConfig(filename)
txt = fileread(filename);

cfg = struct();

cfg.binSize_m           = extractNum(txt, 'Bin size in m:\s*([-\d.]+)');
cfg.origNumSamples      = extractNum(txt, 'Original number of samples:\s*([-\d.]+)');
cfg.origSamplingInterval_ns = extractNum(txt, 'Original sampling interval:\s*([-\d.]+)');
cfg.origRange_ns        = extractNum(txt, 'Original range:\s*([-\d.]+)');

downBlock = regexp(txt, 'Downsampling:(.*?)(Area rotated|$)', 'tokens', 'once');
downTxt = downBlock{1};
cfg.downNumSamples      = extractNum(downTxt, 'Number of samples:\s*([-\d.]+)');
cfg.downSamplingInterval_ns = extractNum(downTxt, 'Sampling interval:\s*([-\d.]+)');
cfg.downRange_ns        = extractNum(downTxt, 'Range:\s*([-\d.]+)');

cfg.rotation_deg = extractNum(txt, 'Area rotated by\s*([-\d.]+)\s*degree');

shift = regexp(txt, 'Area shifted by\s*([-\d.]+)\s*m in x-direction and\s*([-\d.]+)\s*m in y-direction', 'tokens', 'once');
cfg.shiftX_m = str2double(shift{1});
cfg.shiftY_m = str2double(shift{2});

cfg.tz_flag     = extractNum(txt, 'tz_flag:\s*([-\d.]+)');
cfg.followTopo  = extractNum(txt, 'followTopo:\s*([-\d.]+)');
cfg.constV_m_ns = extractNum(txt, 'constV:\s*([-\d.]+)\s*m/ns');
end

function val = extractNum(txt, pattern)
tok = regexp(txt, pattern, 'tokens', 'once');
if isempty(tok)
    val = NaN;
else
    val = str2double(tok{1});
end
end


function saveallslices_geotif(corners,epsg,slice,t,foldername,colperc,tz,sq)

    for numtsl=1:size(slice,3)
        disp(['   ',int2str(numtsl),'/',int2str(size(slice,3))])

        cdata=slice(:,:,numtsl);
    
        cmin=min(cdata(:));
        cmax=max(cdata(:));
    
        tslname = fullfile(foldername,make_fname(numtsl,'.tif',t,tz,sq));
    
        if colperc==0
            writeGeoTIFF(tslname, flipud(cdata), epsg, corners, flipud(gray(16)), [cmin cmax]);
        else
            coldata=sort(cdata(~isnan(cdata)));
            if ~isempty(coldata) && length(coldata)>2
                cmin=coldata(round(length(coldata)/100*colperc));
                cmax=coldata(end-round(length(coldata)/100*colperc));
            end
            writeGeoTIFF(tslname, flipud(cdata), epsg, corners, flipud(gray(16)), [cmin cmax]);
        end
    end
end



function saveallslices_geopng(xgrid,ygrid,slice,topo,t,pfad,colperc,coordtrans,tz,epsg,sq)

dx=abs(xgrid(1,1)-xgrid(1,2));

if ~exist(fullfile(pfad,['georef_png_epsg',int2str(epsg)]),'dir')
    mkdir(fullfile(pfad,['georef_png_epsg',int2str(epsg)]));
end

% topo:
cdata=topo;
cmin=min(cdata(:));
cmax=max(cdata(:));
cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
cdata(isnan(cdata))=0;  % set nan to 0
imwrite(flipud(cdata).*256,parula(256),fullfile(pfad,['georef_png_epsg',int2str(epsg)],'Topo_interp.png'),'Transparency',0);
write_geoPNGW(xgrid,ygrid,coordtrans,fullfile(pfad,['georef_png_epsg',int2str(epsg)],'Topo_interp.pgw'));

% slices
for numtsl=1:size(slice,3)
    disp(['   ',int2str(numtsl),'/',int2str(size(slice,3))])

    cdata=slice(:,:,numtsl);
    cmin=min(cdata(:));
    cmax=max(cdata(:));

    tslname = fullfile(pfad,['georef_png_epsg',int2str(epsg)],make_fname(numtsl,'.png',t,tz,sq));

    if colperc==0
        cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
        cdata(isnan(cdata))=0;  % set nan to 0
        imwrite(flipud(cdata).*256,flipud(gray(256)),tslname,'Transparency',0);
    else
        coldata=sort(cdata(~isnan(cdata)));
        if ~isempty(coldata) && length(coldata)>2
            cmin=coldata(round(length(coldata)/100*colperc));
            cmax=coldata(end-round(length(coldata)/100*colperc));
            range=cmax-cmin;
            cdata=(cdata-cmin)/range;
            cdata(cdata<=0)=0;
            cdata(cdata>=1)=1;
            m=ones(size(cdata));
            m(isnan(cdata))=0;
        end
        im=cdata.*256;
        im(im<=2)=2;
        im(isnan(cdata))=0;  % set nan to 0
        imwrite(flipud(im),flipud(gray(256)),tslname,'Transparency',0);
    end


    % write pngw
    fname = make_fname(numtsl,'.pgw',t,tz,sq);
    if ~exist('coordtrans','var')    % local
        fid=fopen(fullfile(pfad,['georef_png_epsg',int2str(epsg)],fname),'wt');
        fprintf(fid,[num2str(dx),'\n0\n0\n',num2str(-dx),'\n',num2str(min(xgrid(:))),'\n',num2str(max(ygrid(:)))]);
        fclose(fid);
    else % global
        write_geoPNGW(xgrid,ygrid,coordtrans,fullfile(pfad,['georef_png_epsg',int2str(epsg)],fname));
    end

end
end

%%
function fnameStr = make_fname(numtsl,extension,t,tz,sq)
% create filename
if sq==0
    if tz==2 % depth
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'m',extension];
    elseif tz==1 % time
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'ns',extension];
    end
else
    if tz==2 % depth
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'m_sqrt',extension];
    elseif tz==1 % time
        fnameStr = ['Tsl','_',num2str(numtsl,'%3d'),'_t',num2str(t(numtsl),4),'ns_sqrt',extension];
    end
end
end

%%
function []=write_geoPNGW(x,y,coordtrans,filename)

dx=abs(x(1,2)-x(1,1));
dy=abs(y(2,1)-y(1,1));

% determine global coords of upper left and upper right pixel
pix_ol=helmert([min(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
pix_or=helmert([max(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
alpha=atand(abs(pix_ol(2)-pix_or(2))/abs(pix_ol(1)-pix_or(1))); % angle against west
if pix_ol(2)>pix_or(2)
    alpha=-alpha;
end
% determine pixel-lengths in all directions
A=dx*cosd(alpha);
D=dx*sind(alpha);
E=dy*cosd(alpha);
B=dy*sind(alpha);
    
% write pngw
fid=fopen(filename,'wt');
fprintf(fid,[num2str(A),'\n',num2str(D),'\n',num2str(B),'\n',num2str(-E),'\n',num2str(pix_ol(1)),'\n',num2str(pix_ol(2))]);
fclose(fid);
end
