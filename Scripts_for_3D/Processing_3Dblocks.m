clear all
close all
clc


% Processing of 3D-Blocks
% Dr. Tina Wunderlich, CAU 2020-2022, tina.wunderlich@ifg.uni-kiel.de
%
% requires binned data in rectangles 3D_Grid_R*
%
% AVAILABLE PROCESSING STEPS:
%%% Trace attributes:
% Coherence (in+crossline) (Trinks & Hinterleitner, 2020)
% Envelope
%
%%% Migration:
% Isochrone migration 3D with v(z)
%
%%% Others:
% SemblanceSmoothing (modified from Wilken et al. 2019)
% InverseDistanceWeighting-Interpolation
% Linear Interpolation



% Number of rectangles with binned data
rectangles=1; % e.g. 1:3

% is it topography corrected/migrated data?
dsl_topo=0; % =1: with topography correction/migration, =0: time or depth without topography (all sample slices have the same sample points)

% Name of settings.txt, which contains the possible processing steps
% if settings.txt does not exist, a default file will be written
settings='settings_proc3D.txt';

% -------------------------------------------------------------------------
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
addpath(fullfile(curFold,'Processing'),fullfile(curFold,'Migration'),fullfile(curFold,'Subfunctions'));


% read time vector
t=load(fullfile(pfad,['3D_Grid_R',int2str(rectangles(1))],'t.mat'));
t=t.t';
dt=t(2)-t(1);

disp('Reading processing settings:')
if ~exist(fullfile(pfad,settings),'file') % if no settings-file is found: create default file
    disp(['No file ',settings,' found. Creating default file ',settings,'.']);
    
    fid=fopen(fullfile(pfad,settings),'wt');
    fprintf(fid,'do_IsochroneMigration3D_vz 0\naperture[degree] 30\npath_to_vrms.mat .\npath_to_trms.mat .\ninterp 1\n\n');
    fprintf(fid,'do_Coherence 1\nnwavelength 1\n\n');
    fprintf(fid,'do_Envelope 0\n\n');
    fprintf(fid,'do_IDW_Interpolation 0\nradius[m] 1\npower 1\n\n');
    fprintf(fid,'do_linearInterpolation 0\nradLI[m] 1\n\n');
    fprintf(fid,'do_SemblanceSmoothing 0\nsembwin[bins] 5\ndtmin[ns] 0.5\ndtmax[ns] 2\ndtinc[ns] 0.5\nsembexp 2\nsembflag 1\nnormalize 0\nqclip 0.98\n\n');
    fclose(fid);
    
    
    disp(['Edit ',settings,' and start script again.'])
    return;
else
    fid=fopen(fullfile(pfad,settings),'r');
    temp=textscan(fid,'%s%s');
    fclose(fid);
    % get order of processing steps:
    steps=[{'do_IsochroneMigration3D_vz'} {'do_linearInterpolation'} {'do_IDW_Interpolation'} {'do_Coherence'} {'do_Envelope'} {'do_SemblanceSmoothing'}];
    order=zeros(1,length(steps));
    for i=1:length(steps)
        for j=1:length(temp{1})
            if strcmp(temp{1}{j},steps{i})
                order(i)=str2num(temp{2}{j});
            end
        end
    end
    if length(unique(order))~=length(find(order>0))+1
        disp(['Please check order of steps in ',settings,'!']);
        return;
    end
    % get settings for parameters:
    for i=1:length(temp{1})
        if strcmp(temp{1}(i),'aperture[degree]')
            aperture=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'interp')
            interp=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'path_to_vrms.mat')
            p_vrms=temp{2}{i};
        elseif strcmp(temp{1}(i),'path_to_trms.mat')
            p_trms=temp{2}{i};
        elseif strcmp(temp{1}(i),'nwavelength')
            nwavelength=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'sembwin[bins]')
            sembwin=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'dtmin[ns]')
            dtmin=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'dtmax[ns]')
            dtmax=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'dtinc[ns]')
            dtinc=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'sembexp')
            sembexp=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'sembflag')
            sembflag=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'normalize')
            normalize=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'qclip')
            qclip=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'radius[m]')
            radius=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'radLI[m]')
            radLI=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'power')
            power=str2num(temp{2}{i});
        end
    end
    for i=1:length(order(order>0))
        disp([int2str(i),'. ',steps{order==i}])
        if strcmp('do_IsochroneMigration3D_vz',steps{order==i})
            disp(['  aperture[degree] = ',num2str(aperture)])
            disp(['  path_to_vrms.mat = ',p_vrms])
            disp(['  path_to_trms.mat = ',p_trms])
            disp(['  interp = ',num2str(interp)])
        elseif strcmp('do_Coherence',steps{order==i})
            disp(['  nwavelength = ',num2str(nwavelength)])
        elseif strcmp('do_IDW_Interpolation',steps{order==i})
            disp(['  radius[m] = ',num2str(radius)])
            disp(['  power = ',num2str(power)])
        elseif strcmp('do_linearInterpolation',steps{order==i})
            disp(['  radLI[m] = ',num2str(radLI)])
        elseif strcmp('do_SemblanceSmoothing',steps{order==i})
            disp(['  sembwin = ',num2str(sembwin),' bins'])
            disp(['  dtmin = ',num2str(dtmin),' ns'])
            disp(['  dtmax = ',num2str(dtmax),' ns'])
            disp(['  dtinc = ',num2str(dtinc),' ns'])
            disp(['  sembexp = ',num2str(sembexp)])
            disp(['  sembflag = ',num2str(sembflag)])
            disp(['  normalize = ',num2str(normalize)])
            disp(['  qclip = ',num2str(qclip)])
        end
    end
end
% find overall number of processing steps:
numProc_all=length(order(order>0));


%%----------------------------------------------------------------------

% do processing for each rectangle individually
disp('Starting processing for each rectangle...')

% read coordinates
for i=1:length(rectangles)
    disp(['R',int2str(rectangles(i))])
    
    p=['3D_Grid_R',int2str(rectangles(i))];
    temp=load(fullfile(pfad,p,'x.mat'));
    x=temp.x;
    temp=load(fullfile(pfad,p,'y.mat'));
    y=temp.y;
    temp=load(fullfile(pfad,p,'z.mat'));
    z=temp.z;
    temp=load(fullfile(pfad,p,'mask.mat'));
    mask=temp.mask;
    
    if i==1
        dx=x(1,2)-x(1,1);   % dx of 3D-block
    end

    % load data
    temp=load(fullfile(pfad,p,'data.mat'));
    data=temp.data;
    
    %% ---------------------------
    % do processing
    for k=1:length(order(order>0))  % for all processing steps in right order
        
        if strcmp(steps{order==k},'do_IsochroneMigration3D_vz')
            % load vrms and trms:
            temp=load(fullfile(pfad,p_vrms));
            temp1=fieldnames(temp);
            v=getfield(temp,temp1{1});
            temp=load(fullfile(pfad,p_trms));
            temp1=fieldnames(temp);
            vt=getfield(temp,temp1{1});
            [data,zmig]=isochrone_mig_3d_varV(data,x,y,t,v,vt,aperture,interp,1);
        end
        
        if strcmp(steps{order==k},'do_Envelope')
            data=envelope(data);
        end

        if strcmp(steps{order==k},'do_linearInterpolation')
            data=linearInterpolation_3Dcube(data,x,y,radLI);
        end
        
        if strcmp(steps{order==k},'do_IDW_Interpolation')
            if dsl_topo==0
                data=idw3dblock(data,x,y,radius,power);
            else
                data=idw3dblock_topo(data,x,y,radius,power);
            end
        end
        
        if strcmp(steps{order==k},'do_Coherence')
            data=coherence(data,t,nwavelength);
        end
        
        if strcmp(steps{order==k},'do_SemblanceSmoothing')
            data=semblanceSmoothing(data,dt,sembwin,dtmin,dtmax,dtinc,sembexp,sembflag,normalize,qclip);
        end
    end
    
    disp('    Finished!')
    
    % save data
    if ~exist(fullfile(pfad,p,'processed'),'dir')
        mkdir(fullfile(pfad,p,'processed'))
    end
    % copy settings file
    copyfile(fullfile(pfad,settings),fullfile(pfad,p,'processed','ProcessingReadMe.txt'));
    copyfile(fullfile(pfad,p,'coordtrans.mat'),fullfile(pfad,p,'processed','coordtrans.mat'));
    
    save(fullfile(pfad,p,'processed','data.mat'),'data','-v7.3');
    save(fullfile(pfad,p,'processed','mask.mat'),'mask','-v7.3');
    save(fullfile(pfad,p,'processed','x.mat'),'x','-v7.3');
    save(fullfile(pfad,p,'processed','y.mat'),'y','-v7.3');
    save(fullfile(pfad,p,'processed','z.mat'),'z','-v7.3');
    save(fullfile(pfad,p,'processed','t.mat'),'t','-v7.3');
    disp('     Saved!')
    
end
