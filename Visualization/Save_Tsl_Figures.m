clear all
close all
clc


% Read Timeslices and save as georeferenced figures
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% requires Timeslices in MultichannelGPR-format (created with
% make_Timeslices.m)


colperc=1; % Colorscale clipping in percent (if =0: autoscale min-max)

removeBorder=1; % =1: remove border artifacts from interpolation, =0: leave as it is
pix=6; % if removeBorder==1: how many pixels are removed from border around area

medianFilter=1; % do you want to apply a 2D-median filter (1=yes, 0=no)
msize=3; % filter size in pixel

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
            pfad=uigetdir(fn{1}{1},'Choose timeslices folder');
        else
            pfad=uigetdir([],'Choose timeslices folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose timeslices folder'); % path to timeslices-folder

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
            pfad=uigetdir(fn{1}{1},'Choose timeslices folder');
        else
            pfad=uigetdir([],'Choose timeslices folder');
        end
    else
        pfad=uigetdir([],'Choose timeslices folder'); % path to timeslices-folder
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

% load timeslices
if exist(fullfile(pfad,'tsl.mat'),'file')
    temp=load(fullfile(pfad,'tsl.mat'));
    bla=struct2cell(temp);
    tsl=bla{1};
    temp=load(fullfile(pfad,'xgrid.mat'));
    bla=struct2cell(temp);
    xgrid=bla{1};
    temp=load(fullfile(pfad,'ygrid.mat'));
    bla=struct2cell(temp);
    ygrid=bla{1};
    temp=load(fullfile(pfad,'topo.mat'));
    bla=struct2cell(temp);
    topo=bla{1};
    temp=load(fullfile(pfad,'mask.mat'));
    bla=struct2cell(temp);
    mask=bla{1};
    load(fullfile(pfad,'t_startende.mat'));
else
    temp=load(fullfile(pfad,'tsl_interp.mat'));
    bla=struct2cell(temp);
    tsl=bla{1};
    temp=load(fullfile(pfad,'xgrid_interp.mat'));
    bla=struct2cell(temp);
    xgrid=bla{1};
    temp=load(fullfile(pfad,'ygrid_interp.mat'));
    bla=struct2cell(temp);
    ygrid=bla{1};
    temp=load(fullfile(pfad,'topo_interp.mat'));
    bla=struct2cell(temp);
    topo=bla{1};
    temp=load(fullfile(pfad,'mask_interp.mat'));
    bla=struct2cell(temp);
    mask=bla{1};
    load(fullfile(pfad,'t_startende.mat'));
end
if exist(fullfile(pfad,'depth.mat'),'file')
    load(fullfile(pfad,'depth.mat'));
    dsl=1;
else
    dsl=0;
end


%% read info-files
disp('Read info file:')
temp=readlines(fullfile(pfad,'tslinfo.txt'));
for i=1:length(temp)
    if startsWith(temp{i},'Thickness')
        temp2=strsplit(temp{i},' ');
        thickness=str2num(temp2{end-1}); % thickness Tsl
        unit=temp2{end}; % unit (m or ns)
        disp(['   Thickness of Tsl: ',num2str(thickness),' ',unit])
    elseif startsWith(temp{i},'Maximum')
        temp2=strsplit(temp{i},' ');
        maxElevation=str2num(temp2{end-1}); % maximum Elevation
        disp(['   Maximum Elevation: ',num2str(maxElevation),' m'])
    end
end


%% Apply mask on tsl for plotting
for i=1:length(tsl)
    if dsl==0
        maske=mask;
    else
        maske=mask{i}; % one mask for every depth
    end

    tsl{i}=tsl{i}.*maske;
    if medianFilter==1
        tsl{i}=medfilt2(tsl{i},[msize msize]);
    end

    if removeBorder==1 % remove interpolation artifacts around area
        if i==1
            disp('Remove interpolation border around area...')
        end
        dist = chamfer_DT(maske);
        tsl{i}(dist<=pix)=NaN;
    end
end

%% Save timeslices
disp('Save timeslice figures:');

if dsl==0
    maxElevation=[];
end

if exist(fullfile(pfad,'coordtrans.mat'),'file')
    load(fullfile(pfad,'coordtrans.mat'));
    savealltsl(xgrid,ygrid,tsl,topo,t_startende,fullfile(pfad),colperc,dsl,maxElevation,coordtrans);
else
    savealltsl(xgrid,ygrid,tsl,topo,t_startende,fullfile(pfad),colperc,dsl,maxElevation,[1 1 1 1; 2 2 2 2]);
end


% set original path
path(oldpath);

% End of script.


%%
function savealltsl(xgrid,ygrid,tsl,topo,t_startende,pfad,colperc,dsl,maxElevation,coordtrans)
dx=abs(xgrid(1,1)-xgrid(1,2));
for numtsl=1:length(tsl)
    disp(['   ',int2str(numtsl),'/',int2str(length(tsl))])

    % Georeferenced png:
    cdata=tsl{numtsl};
    cmin=min(cdata(:));
    cmax=max(cdata(:));

    if ~exist(fullfile(pfad,'georef'),'dir')
        mkdir(fullfile(pfad,'georef'));
    end


    tslname = fullfile(pfad,'georef',make_fname(numtsl,'.png',dsl,maxElevation,t_startende));

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
    fname = make_fname(numtsl,'.pgw',dsl,maxElevation,t_startende);
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
function fnameStr = make_fname(numtsl,extension,dsl,maxElevation,t)
% depth slice or a time slice ?
if dsl
    dslStr = ['_-_z',num2str(maxElevation - t(numtsl,2),'%5.2f'), '-',...
        num2str(maxElevation - t(numtsl,1),'%5.2f'),'m'];
    % create filename
    fnameStr = ['Tsl','_',num2str(numtsl,'%2d'),'_d',num2str(t(numtsl,1),2),...
    '-',num2str(t(numtsl,2),2),'m',dslStr,extension];
else
    dslStr = [];
    % create filename
    fnameStr = ['Tsl','_',num2str(numtsl,'%2d'),'_t',num2str(t(numtsl,1),2),...
    '-',num2str(t(numtsl,2),2),'ns',dslStr,extension];
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