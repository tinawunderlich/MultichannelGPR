clear all
close all
clc


% Read Sampleslices and save as georeferenced figures
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% requires Sampleslices in MultichannelGPR-format 


colperc=3; % Colorscale clipping in percent (if =0: autoscale min-max)

removeBorder=0; % =1: remove border artifacts from interpolation, =0: leave as it is
pix=6; % if removeBorder==1: how many pixels are removed from border around area

medianFilter=1; % do you want to apply a 2D-median filter (1=yes, 0=no)
msize=3; % filter size in pixel

% use squareroot of amplitudes for visualization?
sq=0; % 1=yes, 0=no

%% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

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
            pfad=uigetdir(fn{1}{1},'Choose sampleslices folder');
        else
            pfad=uigetdir([],'Choose sampleslices folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose sampleslices folder'); % path to sampleslices-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose sampleslices folder');
        else
            pfad=uigetdir([],'Choose sampleslices folder');
        end
    else
        pfad=uigetdir([],'Choose sampleslices folder'); % path to timeslices-folder
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


%% Load data
disp('Loading data...')

% load sampleslices:
if exist(fullfile(pfad,'mask.mat'),'file')
    temp=load(fullfile(pfad,'mask.mat'));
    bla=struct2cell(temp);
    mask=bla{1};
    temp=load(fullfile(pfad,'mask_interp.mat'));
    bla=struct2cell(temp);
    mask_interp=bla{1};
    if exist(fullfile(pfad,'coordtrans.mat'),'file')
    temp=load(fullfile(pfad,'coordtrans.mat'));
    bla=struct2cell(temp);
    coordtrans=bla{1};
    else
        coordtrans=[1 1 1 1; 2 2 2 2];
    end
    temp=load(fullfile(pfad,'xgrid.mat'));
    bla=struct2cell(temp);
    xgrid=bla{1};
    temp=load(fullfile(pfad,'ygrid.mat'));
    bla=struct2cell(temp);
    ygrid=bla{1};
    temp=load(fullfile(pfad,'topo_interp.mat'));
    bla=struct2cell(temp);
    topo_interp=bla{1};
    temp=load(fullfile(pfad,'slice_channelnum.mat'));
    bla=struct2cell(temp);
    slice_channelnum=bla{1};
    temp=load(fullfile(pfad,'slice_profilenum.mat'));
    bla=struct2cell(temp);
    slice_profilenum=bla{1};
    load(fullfile(pfad,'t.mat'));

    i=1;
    while exist(fullfile(pfad,['slice_',int2str(i),'.mat']),'file')
        temp=load(fullfile(pfad,['slice_',int2str(i),'.mat']));
        bla=struct2cell(temp);
        slice{i}=bla{1};
        i=i+1;
    end
else
    disp('No sampleslices found.')
    return;
end





%% Apply mask on tsl for plotting
if removeBorder==1 % remove interpolation artifacts around area
    disp('Remove interpolation border around area...')
    dist = chamfer_DT(mask_interp);
end
for i=1:length(slice)
    disp(['   ',int2str(i),'/',int2str(length(slice))])

    slice{i}=slice{i}.*mask_interp;

    if medianFilter==1
        slice{i}=medianfilt2(slice{i},[msize msize]);
    end

    if removeBorder==1 % remove interpolation artifacts around area
        slice{i}(dist<=pix)=NaN;
    end
end

%% Save timeslices
disp('Save sampleslice figures:');

saveallslices(xgrid,ygrid,slice,topo_interp,t,fullfile(pfad),colperc,coordtrans,sq);

% set original path
path(oldpath);

% End of script.


%%
function saveallslices(xgrid,ygrid,slice,topo,t,pfad,colperc,coordtrans,sq)
dx=abs(xgrid(1,1)-xgrid(1,2));
for numtsl=1:length(slice)
    disp(['   ',int2str(numtsl),'/',int2str(length(slice))])

    % Georeferenced png:
    if sq==1
        cdata=sqrt(slice{numtsl});
    else
        cdata=slice{numtsl};
    end
    cmin=min(cdata(:));
    cmax=max(cdata(:));

    if ~exist(fullfile(pfad,'georef'),'dir')
        mkdir(fullfile(pfad,'georef'));
    end

    tslname = fullfile(pfad,'georef',make_fname(numtsl,'.png',t));

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
    fname = make_fname(numtsl,'.pgw',t);
    if ~exist('coordtrans','var')    % local
        fid=fopen(fullfile(pfad,fname),'wt');
        fprintf(fid,[num2str(dx),'\n0\n0\n',num2str(-dx),'\n',num2str(min(xgrid(:))),'\n',num2str(max(ygrid(:)))]);
        fclose(fid);
    else % global
        write_geoPNGW(xgrid,ygrid,coordtrans,fullfile(pfad,'georef',fname));
    end

end
end

%%
function fnameStr = make_fname(numtsl,extension,t)
    % create filename
    fnameStr = ['Tsl','_',num2str(numtsl,'%2d'),'_t',num2str(t(numtsl),2),'ns',extension];
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