clear all
clc
close all

% Sampleslices are used to create a 3D data block and this block is
% exported as vtk-file for Paraview
%
% Dr. Tina Wunderlich, 2026, tina.wunderlich@ifg.uni-kiel.de
% with the help of Claude.ai
%
% Select the folder with the sampleslices!



% Computer system
platform=2; % Linux=1, Mac=2, Windows=3

time_depth_flag=1;  % 1: time domain, 2: depth domain
maxTZ=40;  % maximum time [ns] or depth [m], if all slices set to Inf

followTopo=0;   % 0: make cube out of horizontal slices starting from t=0 or z=0, ignoring topography
% 1: use topography to create a 3D cube while bending the sampleslices along the topography
constantV=0.1; % if followTopo=1 and time_depth_flag=1: use this constant velocity [m/ns] for converting time to depth

force_division=1; % if =1: a division is forced although the computer memory is large enough for all data

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
disp('Loading data info...')

% load sampleslices:
if exist(fullfile(pfad,'mask_interp.mat'),'file')
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
else
    disp('No sampleslices found.')
    return;
end


name = 'Normalized Amplitude';   % Name in ParaView


% get memory information of system
if platform==1 % Linux
    system('less /proc/meminfo >meminfo.txt');
    fid=fopen('meminfo.txt','r');
    temp=textscan(fid,'%s','Headerlines',2);
    fclose(fid);
    memsize=str2double(temp{1}{2})*1e3;    % bytes
elseif platform==2 % Mac
    [status,cmdout]=system('sysctl hw.memsize');
    memsize=str2double(cmdout(13:end)); % memory size in bytes
elseif platform==3 % Windows
    [user,sys]=memory;
    memsize=str2double(sys.PhysicalMemory.Available);
end

% get size of complete 3d cube + 1000 bytes extra
% info from slices:
[ny,nx] = size(xgrid);
nz=numel(t(t<=maxTZ));

dx = abs(diff(xgrid(1,1:2)));
dy = abs(diff(ygrid(1:2,1)));
dz = abs(diff(t(1:2)));

% origin:
x0 = min(xgrid(:));
y0 = min(ygrid(:));
z0 = min(t);

x = x0 + (0:nx-1) * dx;
y = y0 + (0:ny-1) * dy;
z = z0 + (0:nz-1) * dz;

disp(sprintf('Data cube size is nx=%d / ny=%d / nz=%d',nx,ny,nz))
areasize=nx*ny*nz*8 + 1000; % (each element of the array * bytes in that element (8 for double)) + some overhead (1000)
if memsize/3*2>=areasize
    if force_division==0
        disp('Memory size of computer is larger than required memory! -> Proceed with complete area!')
        num_xrect=1;
        num_yrect=1;
    else
        disp(['Forcing division. Please divide area into rectangles:'])
        num_xrect=str2double(input('How many rectangles in x-direction?  ','s'));
        num_yrect=str2double(input('How many rectangles in y-direction?  ','s'));

        if isempty(num_xrect)
            num_xrect=1;
        end
        if isempty(num_yrect)
            num_yrect=1;
        end
    end
else
    numreq=2;
    while memsize/3*2<=areasize/numreq
        numreq=numreq+1;
    end
    disp(['Required memory size exceeds computer memory! Please divide area into minimum ',int2str(numreq),' rectangles:'])
    num_xrect=str2double(input('How many rectangles in x-direction?  ','s'));
    num_yrect=str2double(input('How many rectangles in y-direction?  ','s'));

    if isempty(num_xrect)
        num_xrect=1;
    end
    if isempty(num_yrect)
        num_yrect=1;
    end
end

if ~exist(fullfile(pfad,'VTK_Export'),'dir')
    mkdir(fullfile(pfad,'VTK_Export'));
end


disp('Creating data cubes...')
% Put data into cube:
if time_depth_flag==1 && followTopo==0
    % horizontal slices and ignoring topography

    disp(sprintf('Dividing area into %d rectangles:',num_xrect*num_yrect))

    % find min/max of coordinates
    minx=min(xgrid(:));
    maxx=max(xgrid(:));
    miny=min(ygrid(:));
    maxy=max(ygrid(:));

    %%% divide area into rectangles
    wid=round((maxx-minx)/num_xrect);  % width of rectangles in x-direction
    hei=round((maxy-miny)/num_yrect);  % height of rectangles in y-direction
    numx=round(nx/num_xrect); % number of samples in one rectangle in x-direction
    numy=round(ny/num_yrect);

    disp(['For each rectangle: Area width ',num2str(wid,4),' m, area height ',num2str(hei,4),' m']);

    disp('--------------------------------------')
    disp('Start reading data in rectangles and creating cubes')
    disp('  ')

    slice=cell(numel(t(t<=maxTZ)),1);
    anz=1;
    for i=1:num_xrect
        for j=1:num_yrect

            tic

            % make grid for this rectangle
            [X,Y,Z]=meshgrid(y((j-1)*numy+1:min(numy*j,ny)),x((i-1)*numx+1:min(numx*i,nx)),z);  % x and y interchanged!!!

            disp(['Reading data in rectangle no. ',int2str(anz),' ...'])

            % initialize data grid (3D) for current rectangle
            data=NaN(size(X));
            % for this rectangle:
            [nx_R, ny_R, nz_R]=size(data);
            x0_R=min(Y(:));
            y0_R=min(X(:));
            z0_R=min(Z(:));

            for ii=1:numel(t(t<=maxTZ)) % for each slice
                if mod(ii,10)==0; disp(['   ',int2str(ii),'/',int2str(length(t(t<=maxTZ)))]); end

                if anz==1 % for first rectangle, read all data and store in variable
                    temp=load(fullfile(pfad,['slice_',int2str(ii),'.mat']));
                    bla=struct2cell(temp);
                    slice{ii}=bla{1};
                    % put into 3D cube of this rectangle:
                    part=slice{ii}((j-1)*numy+1:min(numy*j,ny),(i-1)*numx+1:min(numx*i,nx));
                    data(:,:,ii)=permute(part,[2 1 3]); % change x and y!
                else
                    % put into 3D cube of this rectangle:
                    part=slice{ii}((j-1)*numy+1:min(numy*j,ny),(i-1)*numx+1:min(numx*i,nx));
                    data(:,:,ii)=permute(part,[2 1 3]); % change x and y!
                end
            end

            % make new folder for each rectangle
            if ~exist(fullfile(pfad,'VTK_Export',['Cube_3D_R',int2str(anz)]),'dir')
                mkdir(fullfile(pfad,'VTK_Export',['Cube_3D_R',int2str(anz)]));
            end

            % output name:
            filename = fullfile(pfad,'VTK_Export',['Cube_3D_R',int2str(anz)],'AmplitudeCube.vtk');

            write_vtk_structured_points(filename, ...
                x0_R, y0_R, z0_R, dx, dy, dz, nx_R, ny_R, nz_R, ...
                data, name,[],[],[],[]);

            fprintf('Data cube saved in: %s\n', filename);

            % save figure with rectangle
            hf=figure('Name','CurrentRectangle','visible','off');
            hold off
            imagesc(xgrid(1,:),ygrid(:,1),mask_interp)
            hold on
            plot([y0_R max(Y(:)) max(Y(:)) y0_R y0_R],[x0_R x0_R max(X(:)) max(X(:)) x0_R],'r','Linewidth',2)
            axis xy
            set(gca,'Dataaspectratio',[1 1 1])
            xlabel('x [m]')
            ylabel('y [m]')
            colormap(flipud(gray))
            print(fullfile(pfad,'VTK_Export',['Cube_3D_R',int2str(anz)],'Rectangle_Location.jpg'),'-djpeg');
            close(hf);

            toc
            disp('--------------------------------------')

            % increase number of rectangle
            anz=anz+1;
        end
    end


elseif time_depth_flag==1 && followTopo==1
    % convert time to depth and bend slices along topo

elseif time_depth_flag==2 && followTopo==0
    % horizontal slices in depth and igonring topography

elseif time_depth_flag==2 && followTopo==1
    % bend slices along topography


end

% exporting info file:
fid=fopen(fullfile(pfad,'VTK_Export','Info_vtk.txt'),'wt');
fprintf(fid,'Settings:\ntime_depth_flag=%d\nmaxTZ=%.2f\nfollowTopo=%d\nconstantV=%.2f\n',time_depth_flag,maxTZ,followTopo,constantV);

fclose(fid);


fprintf('Open the file in ParaView with: File -> Open.\n');


% set original path
path(oldpath);

% End of script.

%% =========================================================

function write_vtk_structured_points(fname, ...
    x0, y0, z0, dx, dy, dz, nx, ny, nz, ...
    scalar_data, scalar_name, ...
    Vx, Vy, Vz, vector_name)
%WRITE_VTK_STRUCTURED_POINTS  Schreibt VTK Legacy ASCII Datei.
%
%  Das Gitter ist ein STRUCTURED_POINTS Gitter (aequidistant, kartesisch).
%  ParaView liest dieses Format direkt ohne Plugin.
%
%  scalar_data : [ny x nx x nz] double  (leer [] = kein Skalar)
%  Vx,Vy,Vz   : [ny x nx x nz] double  (leer [] = kein Vektor)

npts = nx * ny * nz;

fid = fopen(fname, 'w');
if fid == -1
    error('Cannot open file: %s', fname);
end

% ---------- HEADER ----------
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'MATLAB VTK Export - %s\n', datestr(now));
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
fprintf(fid, 'ORIGIN %g %g %g\n', x0, y0, z0);
fprintf(fid, 'SPACING %g %g %g\n', dx, dy, dz);
fprintf(fid, 'POINT_DATA %d\n', npts);

% ---------- SKALARDATEN ----------
if ~isempty(scalar_data)
    % VTK erwartet Reihenfolge: x aendert sich zuerst, dann y, dann z
    % meshgrid liefert (y,x,z) -> permutieren auf (x,y,z)
    S = permute(scalar_data, [1 2 3]);   % nx x ny x nz %%% [2 1 3] falsch!
    S = S(:);                            % spaltenweise = x-zuerst

    fprintf(fid, 'SCALARS %s float 1\n', strrep(scalar_name,' ','_'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%g\n', S);
end

% ---------- VEKTORDATEN ----------
if ~isempty(Vx)
    Ux = permute(Vx, [2 1 3]); Ux = Ux(:);
    Uy = permute(Vy, [2 1 3]); Uy = Uy(:);
    Uz = permute(Vz, [2 1 3]); Uz = Uz(:);

    fprintf(fid, 'VECTORS %s float\n', strrep(vector_name,' ','_'));
    fprintf(fid, '%g %g %g\n', [Ux, Uy, Uz]');
end

fclose(fid);
end