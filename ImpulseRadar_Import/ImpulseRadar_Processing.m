% Script for reading ImpulseRadar-datafiles (*.iprb, *.iprh, *.cor, ...)
% and processing them (no binning!)
%
% Dr. Tina Wunderlich, CAU Kiel 2025, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files:
% *.iprb: Profile data file
% *.iprh: header file
% *.cor: coordinate file
% *.time: timing file
% *.mkr: marker file
%
% requires MATLAB-files in following folders (path will be temporarily
% set): Processing, Export_Import, Subfunctions, Migration


clear all
close all
clc

% Computer system
platform=2; % Linux=1, Mac=2, Windows=3

% Select number of profiles:
profile_min=1;  % minimum profile number
profile_max=2;  % maximum profile number

% Coordinate settings:
add_Yoffset=0;    % add a constant offset to the y-antenna-GNSS-offset (e.g. due to non-vertical GNSS-stick), set =0 if not applicable
% negative if GNSS is tilted towards front, positive if tilted towards back
GNSS_height=0; % Height of GNSS antenna above ground [m]
utmzone=0; % set UTM Zone for measurements with GNSS, set =0 for measurements with total station
filter_coords=1; % =1 if you want to delete coordinates at the same position, =0 no filtering
smooth_coords=25; % width of moving median window. No smoothing if =0. It is recommended to smooth the GNSS coordinates before applying the offsets for each channel

%%% Processing steps before binning in txt-file (a default file with this
%%% name is created, if not existing)
settings='settings.txt';

%%% Raw data is only read if userawdata=0! Otherwise mat-files of raw
%%% data in folder profiles2mat are used.
userawdata=0;  % if =1: use aready read in raw data and apply new processing steps

% save radargrams.mat of processed profile data in one variable (requires large memory -> probably not
% working for every data set)?
rad=1; % if =0, processed data will be saved in folder proc, but not all radargrams in one file radargrams.mat
% if =1, processed data will be saved in radargrams.mat...


%--------------------------------------------------------------------------
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
            foldername=uigetdir(fn{1}{1},'Choose data folder');
        else
            foldername=uigetdir([],'Choose data folder');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Choose data folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose data folder');
        else
            foldername=uigetdir([],'Choose data folder');
        end
    else
        foldername=uigetdir([],'Choose data folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% get name
temp=dir(fullfile(foldername,'/*.cor')); % GNSS file
if isempty(temp)
    temp=dir(fullfile(foldername,'/*.tsp')); % total station file
    utmzone=0;
else
    % GNSS: check if utmzone set
    if utmzone==0
        disp('Please set correct UTM Zone and start again.')
        return;
    end
end
tempname=strsplit(temp(end).name,'_'); % Name of data files without '_???.cor/tsp
name=[tempname{1}];
for i=2:length(tempname)-1
    name=[name,'_',tempname{i}];
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Processing'),fullfile(curFold,'Subfunctions'),fullfile(curFold,'Migration'));


% file name parts
[pathstr,fname,ext]=fileparts(name);    % divide name into parts


%% Reading settings from file
disp('--------------------------------------')
disp('Reading processing settings:')
if ~exist(fullfile(foldername,settings),'file') % if no settings-file is found: create default file
    disp(['No file ',settings,' found. Creating default file ',settings,'.']);

    fid=fopen(fullfile(foldername,settings),'wt');
    fprintf(fid,'do_amplitudeOffset 1\ntstart[ns] 0\ntend[ns] 100\n\n');
    fprintf(fid,'do_badTraceRemoval 0\nminfactor 3\nmaxfactor 3\n\n');
    fprintf(fid,'do_constTraceDist 4\ndx[m] 0.02\n\n');
    fprintf(fid,'do_traceInterpolation 5\ngap 30\n\n');
    fprintf(fid,'do_medfilt 0\nnumsamp 3\n\n');
    fprintf(fid,'do_medfilt_x 0\nnumsamp_x 3\ntstart_x 0\n\n');
    fprintf(fid,'do_sphericalDivergence 0\n\n');
    fprintf(fid,'do_attenuationCorrection 0\nsigma[S/m] 0.002\neps 15\n\n');
    fprintf(fid,'do_spectralWhitening 0\nfmin_sw[MHz] 100\nfmax_sw[MHz] 600\nalpha 0.01\n\n');
    fprintf(fid,'do_bandpass 6\nfstart[MHz] 100\nfend[MHz] 600\n\n');
    fprintf(fid,'do_migration2d_vz 0\nv-file[m/ns] 0.1\ntv-file[ns] 0\naperture_m[degree] 30\n\n');
    fprintf(fid,'do_topomig2d 0\nv[m/ns] 0.1\nflag 1\naperture_t[degree] 30\n\n');
    fprintf(fid,'do_cutTWT 0\ntmax[ns] 80\n\n');
    fprintf(fid,'do_khighpass 0\nkcutoff[1/m] 0.1\n\n');
    fprintf(fid,'do_applygain 7\ng1 -20\ng2 0\ng3 10\ng4 20\ng5 30\n\n');
    fprintf(fid,'do_normalization 0\nqclip 0.98\n\n');
    fprintf(fid,'do_removeMeanMedianTrace 0\nmeanmedian 1\nnumtraces 0\n\n');
    fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
    fprintf(fid,'do_t0CorrectionThreshold 0\nthreshold -1\n\n');
    fprintf(fid,['do_t0CorrectionReferencetrace 0\nprofilenumber ',int2str(profile_min),'\nchannelnumber 1\ntracenumber 1\nt0[ns] 0.1\n\n']);
    fprintf(fid,'do_crossprofileshift 3\nmaxshift_cps 30\n\n');
    fprintf(fid,['do_channelshift 2\nOneOrAll 2\nmaxshift 30\nrefprofile ',int2str(profile_min),'\nshiftsamples\n']);
    fclose(fid);

    disp('Edit settings.txt and start script again.')
    return;
else
    fid=fopen(fullfile(foldername,settings),'r');
    temp=textscan(fid,'%s%s');
    fclose(fid);
    % get order of processing steps:
    steps=[{'do_amplitudeOffset'} {'do_badTraceRemoval'} {'do_constTraceDist'} {'do_traceInterpolation'}...
        {'do_medfilt'} {'do_medfilt_x'} {'do_sphericalDivergence'} {'do_attenuationCorrection'}...
        {'do_spectralWhitening'} {'do_bandpass'} {'do_migration2d_vz'} {'do_topomig2d'}...
        {'do_cutTWT'} {'do_khighpass'}  {'do_applygain'} {'do_normalization'} {'do_removeMeanMedianTrace'}...
        {'do_t0shift'} {'do_t0CorrectionThreshold'} {'do_t0CorrectionReferencetrace'} {'do_crossprofileshift'} {'do_channelshift'}]; 
    order=zeros(1,length(steps));
    for i=1:length(steps)
        for j=1:length(temp{1})
            if strcmp(temp{1}{j},steps{i})
                temp3=temp{2}(j);
                order(i)=str2double(temp3{1});
            end
        end
    end
    if (length(unique(order))~=length(find(order>0))+1) || max(order)==length(unique(order))
        disp(['Please check order of steps in ',settings,'!']);
        return;
    end
    % get settings for parameters:
    for i=1:length(temp{1})
        if strcmp(temp{1}(i),'tstart[ns]')
            temp3=temp{2}(i);
            params.tstart=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'tend[ns]')
            temp3=temp{2}(i);
            params.tend=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'dx[m]')
            temp3=temp{2}(i);
            params.dist=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'t0s[ns]')
            temp3=temp{2}(i);
            params.t0s=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'gap')
            temp3=temp{2}(i);
            params.gap=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'numsamp')
            temp3=temp{2}(i);
            params.numsamp=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'numsamp_x')
            temp3=temp{2}(i);
            params.numsamp_x=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'tstart_x')
            temp3=temp{2}(i);
            params.tstart_x=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'maxshift_cps')
            temp3=temp{2}(i);
            params.maxshift_cps=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'minfactor')
            temp3=temp{2}(i);
            params.minfactor=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'maxfactor')
            temp3=temp{2}(i);
            params.maxfactor=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'fmin_sw[MHz]')
            temp3=temp{2}(i);
            params.fmin_sw=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'fmax_sw[MHz]')
            temp3=temp{2}(i);
            params.fmax_sw=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'alpha')
            temp3=temp{2}(i);
            params.alpha=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'aperture_m[degree]')
            temp3=temp{2}(i);
            params.aperture_m=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'aperture_t[degree]')
            temp3=temp{2}(i);
            params.aperture_t=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'kcutoff[1/m]')
            temp3=temp{2}(i);
            params.kcutoff=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'v[m/ns]')
            temp3=temp{2}(i);
            params.v=temp3{1};
            params.vname=params.v;
            if ~isempty(str2num(params.v)) % one number
                params.v=str2num(params.v);
            else % file name
                temp4=load(fullfile(foldername,params.v));
                temp2=fieldnames(temp4);
                params.v=getfield(temp4,temp2{1});
            end
        elseif strcmp(temp{1}(i),'flag')
            temp3=temp{2}(i);
            params.flagtopo=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'v-file[m/ns]')
            temp3=temp{2}(i);
            params.vfile=temp3{1};
        elseif strcmp(temp{1}(i),'tv-file[ns]')
            temp3=temp{2}(i);
            params.tfile=temp3{1};
            % read v(t)
            params.vfilename=params.vfile;
            params.tfilename=params.tfile;
            if ~isempty(str2num(params.vfile)) % one number
                params.vfile=[str2num(params.vfile); str2num(params.vfile)];
                params.tfile=[0; 100];
            else
                temp4=load(fullfile(foldername,params.vfile));
                temp2=fieldnames(temp4);
                params.vfile=getfield(temp4,temp2{1});
                temp4=load(fullfile(foldername,params.tfile));
                temp2=fieldnames(temp4);
                params.tfile=getfield(temp4,temp2{1});
                % make column vectors
                if length(params.vfile(:,1))<length(params.vfile(1,:))
                    params.vfile=params.vfile';
                end
                if length(params.tfile(:,1))<length(params.tfile(1,:))
                    params.tfile=params.tfile';
                end
            end
        elseif strcmp(temp{1}(i),'g1')
            temp3=temp{2}(i);
            params.g1=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'g2')
            temp3=temp{2}(i);
            params.g2=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'g3')
            temp3=temp{2}(i);
            params.g3=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'g4')
            temp3=temp{2}(i);
            params.g4=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'g5')
            temp3=temp{2}(i);
            params.g5=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'threshold')
            temp3=temp{2}(i);
            params.threshold=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'profilenumber')
            temp3=temp{2}(i);
            params.profilenum=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'channelnumber')
            temp3=temp{2}(i);
            params.channelnum=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'tracenumber')
            temp3=temp{2}(i);
            params.tracenum=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'t0[ns]')
            temp3=temp{2}(i);
            params.t0=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'maxshift')
            temp3=temp{2}(i);
            params.maxshift=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'OneOrAll')
            temp3=temp{2}(i);
            params.oneall=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'refprofile')
            temp3=temp{2}(i);
            params.refprofile=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'shiftsamples')
            if length(temp{2})>=i && ~isempty(temp{2}{i})
                if channels==8
                    params.shiftsamples=[str2double(temp{2}{i}) str2double(temp{1}{i+1}) str2double(temp{2}{i+1}) str2double(temp{1}{i+2}) str2double(temp{2}{i+2}) str2double(temp{1}{i+3}) str2double(temp{2}{i+3}) str2double(temp{1}{i+4})];
                else
                    params.shiftsamples=[str2double(temp{2}{i}) str2double(temp{1}{i+1}) str2double(temp{2}{i+1}) str2double(temp{1}{i+2}) str2double(temp{2}{i+2}) str2double(temp{1}{i+3}) str2double(temp{2}{i+3}) str2double(temp{1}{i+4})...
                        str2double(temp{2}{i+4}) str2double(temp{1}{i+5}) str2double(temp{2}{i+5}) str2double(temp{1}{i+6}) str2double(temp{2}{i+6}) str2double(temp{1}{i+7}) str2double(temp{2}{i+7}) str2double(temp{1}{i+8})];
                end
            else
                params.shiftsamples=[];
            end
        elseif strcmp(temp{1}(i),'sigma[S/m]')
            temp3=temp{2}(i);
            params.sigma=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'eps')
            temp3=temp{2}(i);
            params.eps=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'tmax[ns]')
            temp3=temp{2}(i);
            params.tmax=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'meanmedian')
            temp3=temp{2}(i);
            params.meanmedian=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'numtraces')
            temp3=temp{2}(i);
            params.numtraces=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'fstart[MHz]')
            temp3=temp{2}(i);
            params.fstart=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'fend[MHz]')
            temp3=temp{2}(i);
            params.fend=str2double(temp3{1});
        elseif strcmp(temp{1}(i),'qclip')
            temp3=temp{2}(i);
            params.qclip=str2double(temp3{1});
        end
    end
    prepare_shift=0;
    prepare_t0=0;
    for i=1:length(order(order>0))
        disp([int2str(i),'. ',steps{order==i}])
        if strcmp('do_amplitudeOffset',steps{order==i})
            disp(['  tstart = ',num2str(params.tstart),' ns'])
            disp(['  tend = ',num2str(params.tend),' ns'])
        elseif strcmp('do_badTraceRemoval',steps{order==i})
            disp(['  minfactor = ',num2str(params.minfactor)])
            disp(['  maxfactor = ',num2str(params.maxfactor)])
        elseif strcmp('do_traceInterpolation',steps{order==i})
            disp(['  gap = ',num2str(params.gap)])
        elseif strcmp('do_t0CorrectionThreshold',steps{order==i})
            disp(['  threshhold = ',num2str(params.threshold)])
        elseif strcmp('do_medfilt',steps{order==i})
            disp(['  numsamp = ',num2str(params.numsamp)])
        elseif strcmp('do_medfilt_x',steps{order==i})
            disp(['  numsamp_x = ',num2str(params.numsamp_x)])
            disp(['  tstart_x = ',num2str(params.tstart_x)])
        elseif strcmp('do_constTraceDist',steps{order==i})
            disp(['  dist = ',num2str(params.dist)])
        elseif strcmp('do_spectralWhitening',steps{order==i})
            disp(['  fmin_sw = ',num2str(params.fmin_sw)])
            disp(['  fmax_sw = ',num2str(params.fmax_sw)])
            disp(['  alpha = ',num2str(params.alpha)])
        elseif strcmp('do_applygain',steps{order==i})
            disp(['  gain = [',num2str(params.g1),' ',num2str(params.g2),' ',num2str(params.g3),' ',num2str(params.g4),' ',num2str(params.g5),'] dB'])
        elseif strcmp('do_removeMeanMedianTrace',steps{order==i})
            disp(['  meanmedian = ',num2str(params.meanmedian)])
            disp(['  numtraces = ',num2str(params.numtraces)])
            meanmedianflag=1;
        elseif strcmp('do_t0CorrectionReferencetrace',steps{order==i})
            disp(['  profilenumber = ',num2str(params.profilenum)])
            disp(['  channelnumber = ',num2str(params.channelnum)])
            disp(['  tracenumber = ',num2str(params.tracenum)])
            disp(['  t0 = ',num2str(params.t0),' ns'])
            prepare_t0=1;
        elseif strcmp('do_cutTWT',steps{order==i})
            disp(['  tmax = ',num2str(params.tmax),' ns'])
        elseif strcmp('do_khighpass',steps{order==i})
            disp(['  kcutoff = ',num2str(params.kcutoff),' 1/m'])
        elseif strcmp('do_migration2d_vz',steps{order==i})
            disp(['  v-file[m/ns] = ',params.vfilename])
            disp(['  tv-file[ns] = ',params.tfilename])
            disp(['  aperture = ',num2str(params.aperture_m)])
        elseif strcmp('do_topomig2d',steps{order==i})
            disp(['  v = ',num2str(params.vname),' m/ns'])
            disp(['  flag = ',num2str(params.flagtopo)])
            disp(['  aperture = ',num2str(params.aperture_t)])
        elseif strcmp('do_t0shift',steps{order==i})
            disp(['  t0s = ',num2str(params.t0s),' ns'])
        elseif strcmp('do_crossprofileshift',steps{order==i})
            disp(['  maxshift_cps = ',num2str(params.maxshift_cps)])
        elseif strcmp('do_channelshift',steps{order==i})
            disp(['  maxshift = ',num2str(params.maxshift)])
            disp(['  OneOrAll = ',num2str(params.oneall)])
            disp(['  Reference profile = ',num2str(params.refprofile)])
            if ~isempty(params.shiftsamples)
                disp(['  shiftsamples = ',int2str(params.shiftsamples)])
            end
            prepare_shift=1;
        elseif strcmp('do_attenuationCorrection',steps{order==i})
            disp(['  sigma = ',num2str(params.sigma),' S/m'])
            disp(['  eps = ',num2str(params.eps)])
        elseif strcmp('do_bandpass',steps{order==i})
            disp(['  fstart = ',num2str(params.fstart),' MHz'])
            disp(['  fend = ',num2str(params.fend),' MHz'])
        elseif strcmp('do_normalization',steps{order==i})
            disp(['  qclip = ',num2str(params.qclip)])
        end
    end
    disp('--------------------------------------')
end


%% Load all profiles, do processing and save as mat-files

% Profile numbers
numbers=profile_min:profile_max;

if userawdata==0  % first run of program -> read all profiles
    disp('Read original profiles and save in folder profiles2mat -> raw data')
    % make folder to save raw data
    if ~exist(fullfile(foldername,'profiles2mat'),'dir')
        mkdir(fullfile(foldername,'profiles2mat'));
    end
    % read profile data and coordinates
    not=zeros(1,length(numbers));
    lnum=length(numbers);

    tempx=cell(lnum,1);
    tempy=cell(lnum,1);
    tempz=cell(lnum,1);
    for i=1:lnum
        % load data and coordinates
        [traces,dt,ns,tempx{i},tempy{i},tempz{i},channels]=readImpulseRadar(foldername,name,numbers(i),add_Yoffset,GNSS_height,utmzone,filter_coords,smooth_coords);
        traces=single(traces); % convert to single for saving memory
        % delete traces with NaN-coordinates
        del=find(isnan(tempx{i}));
        traces(:,del)=[];
        tempx{i}(del)=[];
        tempy{i}(del)=[];
        tempz{i}(del)=[];
        if ~isempty(dt)
            % resort into matrices (without cells)
            numtr=length(tempx{i}); % number of all traces per profile
            numtrch=numtr/channels; % number of traces per channel
            info=zeros(9,numtr);   % profilenum, tracenum per channel, channelnum, x, y, z, dt, ns (per trace), tracenum per profile
            info(1,:)=zeros(1,length(tempx{i}))+numbers(i); % profilenumber
            info(4:6,:)=[tempx{i}; tempy{i}; tempz{i}]; % x,y,z
            info(9,:)=1:numtr; % tracenumber per profile
            for ii=1:channels
                info(2:3,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[1:numtrch; zeros(1,numtrch)+ii]; % trace number per channel, channelnumber
                info(7:8,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[zeros(1,numtrch)+dt; zeros(1,numtrch)+ns]; % dt, ns
            end
            % save profile (raw data)
            infoname=fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info.mat']);
            trname=fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']);
            parsave(infoname,info);
            parsave(trname,traces);
            disp(['   ',int2str(numbers(i))])
        end
    end
    disp('All raw data read and saved in folder: profiles2mat.')
end
disp('--------------------------------------')

%% For some processing steps: Prepare reference data
if prepare_shift==1 && params.oneall==1
    if isempty(params.shiftsamples)
        disp('For Channelshift: Calculate shiftsamples from reference profile');
        temp2=load(fullfile(foldername,'profiles2mat',[name,'_',int2str(params.refprofile),'_info.mat']));
        info=temp2.info;
        temp=load(fullfile(foldername,'profiles2mat',[name,'_',int2str(params.refprofile),'.mat']));
        traces=temp.traces;
        clear temp temp2;
        t=0:info(7,1):info(7,1)*(info(8,1)-1);
        % do processing steps before:
        order2=zeros(size(order));
        for i=1:max(order(order>0))
            if ~strcmp('do_channelshift',steps{order==i})
                order2(order==i)=i;
            else
                break;
            end
        end
        traces=processing(steps,order2,traces,info,t,params);
        % do channelshift:
        test=cell(max(info(3,:)),1);
        for i=1:max(info(3,:))
            test{i}=traces(:,info(3,:)==i);
        end
        [~,params.shiftsamples]=channelshift(test,params.maxshift); % shiftsamples for reference profile

        % save shiftsamples in settings.txt
        fid=fopen(fullfile(foldername,settings),'a');
        %fprintf(fid,'shiftsamples\t');
        for ch=1:length(test)
            fprintf(fid,'%d\t',params.shiftsamples(ch));
        end
        fclose(fid);
        clear traces;
        clear test;
        disp('--------------------------------------')
    end
else
    params.shiftsamples=[];
end
if prepare_t0==1
    disp('For t0correction: Read reference trace');
    temp=load(fullfile(foldername,'profiles2mat',[name,'_',int2str(params.profilenum),'.mat']));
    traces=temp.traces;
    temp2=load(fullfile(foldername,'profiles2mat',[name,'_',int2str(params.profilenum),'_info.mat']));
    info=temp2.info;
    clear temp temp2;
    t=0:info(7,1):info(7,1)*(info(8,1)-1);
    % do processing steps before:
    order2=zeros(size(order));
    for i=1:max(order(order>0))
        if ~strcmp('do_t0CorrectionReferencetrace',steps{order==i})
            order2(order==i)=i;
        else
            break;
        end
    end
    % get reference trace
    params.reftrace=traces(:,info(2,:)==params.tracenum & info(3,:)==params.channelnum); % Reference trace
    clear traces;
    disp('--------------------------------------')
else
    params.reftrace=[];
end


%% Processing of data

disp('Start processing of raw data')
% make folder for processed data
if ~exist(fullfile(foldername,'profiles2mat','proc'),'dir')
    mkdir(fullfile(foldername,'profiles2mat','proc'));
end

% find maxz of complete data set
params.maxz=0;
params.dx=[];
params.minz=[];
for i=1:length(numbers) % (is not working as parfor!)
    % load raw data
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info.mat']));
    end
    if max(info(6,:))>params.maxz
        params.maxz=max(info(6,:));
    end
    if isempty(params.minz)
        params.minz=min(info(6,:));
    else
        if min(info(6,:))<params.minz
            params.minz=min(info(6,:));
        end
    end
    if isempty(params.dx)
        params.dx=mean(sqrt(diff(info(4,:)).^2+diff(info(5,:)).^2));
    end
end
params.maxz=params.maxz+0.5;
% get mean v:
if any(strcmp(steps,'do_migration2d_vz'))
    meanv=mean(params.vfile);
elseif any(strcmp(steps,'do_topomig2d'))
    meanv=mean(params.v);
else
    meanv=0.1;
end
params.minz=params.maxz-info(7,1)*info(8,1)/2*meanv-0.5-(params.maxz-params.minz);

% if do_crossprofileshift -> determine shiftsamples now!
if any(strcmp(steps,'do_crossprofileshift')) && order(strcmp(steps,'do_crossprofileshift'))>0
    if ~exist(fullfile(foldername,'shiftsamples_cps.txt'),'file')
        disp('Determine shiftsamples for crossprofileshift...')
        prof=NaN(length(numbers),1);
        for i=1:length(numbers) % (is not working as parfor!)
            % load raw data
            if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']),'file')
                mat=matfile(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']));
                load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info.mat']));

                cps{i}=mat.traces(:,find(info(3,:)==1)); % get all first channel traces
                prof(i,1)=numbers(i);
            end
        end

        % delete empty cells
        cps(isnan(prof))=[];
        prof(isnan(prof))=[];

        % determine shiftsamples_cps
        [~,params.shiftsamples_cps]=channelshift(cps,params.maxshift_cps);

        params.shiftsamples_cps=[prof params.shiftsamples_cps]; % profilenumber shiftsamples

        fid=fopen(fullfile(foldername,'shiftsamples_cps.txt'),'wt');
        fprintf(fid,'%d\t%d\n',params.shiftsamples_cps');
        fclose(fid);
        clear cps;
        disp('   Done and saved in shiftsamples_cps.txt')
    else
        disp('   Reading shiftsamples_cps from shiftsamples_cps.txt')
        params.shiftsamples_cps=load(fullfile(foldername,'shiftsamples_cps.txt'));
    end
end


anz=1;
for i=1:length(numbers)
    % load raw data
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']),'file')
        tic
        disp(['Profile #',int2str(numbers(i))]);
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'.mat']));
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info.mat']));

        params.dt=info(7,1);
        t=0:params.dt:params.dt*(info(8,1)-1);
        profileinfotemp{i}=[numbers(i) abs(t(2)-t(1)) info(8,1) max(info(3,:)) max(info(2,:))]; % profilnumber, dt, ns, channels, numtraces-per-channel


        % Plot figure in invisible mode
        b=1:profileinfotemp{i}(5):profileinfotemp{i}(5)*profileinfotemp{i}(4); % first trace of every channel
        fhProc=figure('Name','Processing control','Visible','off','PaperPosition',[0 0 50 50]);
        subplot(4,3,1)
        hold on
        plot(t,traces(:,b))
        grid on
        xlabel('t [ns]')
        title('Original data')

        % processing:
        [traces,ns,t,zmig,info]=processing(steps,order,traces,info,t,params,fhProc,b);

        % new profileinfo update
        profileinfotemp{i}=[numbers(i) abs(t(2)-t(1)) info(8,1) max(info(3,:)) max(info(2,:))]; % profilnumber, dt, ns, channels, numtraces-per-channel


        if rad==1
            % write all profiles in one variable
            for ii=1:max(info(3,:)) % for each channel
                radargrams{anz}=traces(:,info(3,:)==ii);
                global_coords{anz}=[info(4,info(3,:)==ii)' info(5,info(3,:)==ii)' info(6,info(3,:)==ii)'];
                x{anz}=[0; cumsum(sqrt(diff(global_coords{anz}(:,1)).^2+diff(global_coords{anz}(:,2)).^2))];
                anz=anz+1;
            end
        end

        % save processed data, t and zmig
        parsave(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(numbers(i)),'.mat']),traces);
        parsave(fullfile(foldername,'profiles2mat','proc','t.mat'),t);
        parsave(fullfile(foldername,'profiles2mat','proc','zmig.mat'),zmig);
        parsave(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat']),info);

        % save figure
        parprintfigure(fullfile(foldername,'profiles2mat','proc',['Processing_',int2str(numbers(i)),'.jpg']),fhProc);
        close(fhProc);
        disp('     ...saved')
        clear traces;
        clear info;
        timeiteration=toc;
        disp(['Elapsed time is ',num2str(timeiteration,5),' s.'])
    end
end
if rad==1
    % save all profiles in one variable
    save(fullfile(foldername,'profiles2mat','proc','radargrams.mat'),'radargrams','-v7.3');
    save(fullfile(foldername,'profiles2mat','proc','t.mat'),'t','-v7.3');
    save(fullfile(foldername,'profiles2mat','proc','x.mat'),'x','-v7.3');
    save(fullfile(foldername,'profiles2mat','proc','global_coords.mat'),'global_coords','-v7.3');
    if ~isempty(zmig)
        save(fullfile(foldername,'profiles2mat','proc','zmig.mat'),'zmig','-v7.3');
    end

    fid=fopen(fullfile(foldername,'profiles2mat','proc','radargrams.txt'),'wt');
    fprintf(fid,'Radargrams.mat contains channels\n');
    fprintf(fid,' %d\t',channels);
    fprintf(fid,'\nof profiles\n');
    fprintf(fid,' %d\n',numbers);
    fclose(fid);
end

% save settings in proc-folder
copyfile(fullfile(foldername,settings),fullfile(foldername,'profiles2mat','proc','ProcessingReadMe.txt'));

% update profile info
for i=1:length(numbers)
    if ~isempty(profileinfotemp{i})
        profileinfo(i,:)=profileinfotemp{i};
    end
end
profileinfo(profileinfo(:,2)==0,:)=[];
profileinfo(:,3)=ns;
if ~isempty(zmig)
    profileinfo(:,2)=abs(zmig(1)-zmig(2)); % replace dt[ns] by dz[m]
end
% append new profileinfo to old profileinfo
if exist(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'),'file')
    temp1=load(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'));
    temp1=temp1.profileinfo; % old info
    for ii=1:length(profileinfo(:,1))
        % check if already in old profileinfo
        if any(temp1(:,1)==profileinfo(ii,1))
            temp1(temp1(:,1)==profileinfo(ii,1),:)=profileinfo(ii,:); % replace old
        else % append
            temp1=[temp1; profileinfo(ii,:)];
        end
    end
    profileinfo=sortrows(temp1,1);
end
save(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'),'profileinfo','-v7.3');
disp('--------------------------------------')


% set original path
path(oldpath);

% End of script.

%--------------------------------------------------------------------------
%%% subfunctions:
function [datatraces,ns,t,zmig,info]=processing(steps,order,datatraces,info,t,params,fh,b)

ns=length(datatraces(:,1));
zmig=[];

for k=1:length(order(order>0))  % for all processing steps in right order


    if strcmp(steps{order==k},'do_constTraceDist')
        disp('  ...constant trace distance')
        channels=max(info(3,:)); % number of channels
        temp=cell(channels,1);
        x=cell(channels,1);
        gc=cell(channels,1);
        for ch=1:channels % for each channel
            [temp{ch},x{ch},gc{ch}]=constTraceDist(datatraces(:,info(3,:)==ch),params.dist,cumsum([0 sqrt((diff(info(4,info(3,:)==ch)).^2+diff(info(5,info(3,:)==ch)).^2))]),info(4:6,info(3,:)==ch)');
        end

        numtrch=min(cellfun(@(x) length(x),x)); % minimum number of traces per channel

        datatraces_temp=zeros(length(temp{1}(:,1)),numtrch*channels);
        info_temp=zeros(9,numtrch*channels);
        info_temp(1,:)=info(1,1); % profile number
        info_temp(9,:)=1:length(info_temp(9,:)); % overall trace number
        for ii=1:channels
            info_temp(2:3,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[1:numtrch; zeros(1,numtrch)+ii]; % trace number per channel, channelnumber
            info_temp(7:8,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=[zeros(1,numtrch)+info(7,1); zeros(1,numtrch)+info(8,1)]; % dt, ns
            info_temp(4:6,info_temp(3,:)==ii)=gc{ii}(1:numtrch,:)';
            datatraces_temp(:,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=temp{ii}(:,1:numtrch);
        end

        % replace datatraces and info with shorter profiles:
        datatraces=single(datatraces_temp);
        info=info_temp;
        clear datatraces_temp;
        clear info_temp;
        clear temp;

        b=1:numtrch:numtrch*max(info(3,:)); % first trace of every channel
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Constant trace distance')
        end
    end

    if strcmp(steps{order==k},'do_traceInterpolation')
        disp('  ...Trace interpolation')
        for ii=1:max(info(3,:))
            datatraces(:,info(3,:)==ii)=interpolation(datatraces(:,info(3,:)==ii),params.gap);
        end
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Trace interpolation')
        end
    end

    if strcmp(steps{order==k},'do_medfilt')
        disp('  ...medfilt')
        datatraces=medfilt(datatraces,params.numsamp);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Median Filter')
        end
    end

    if strcmp(steps{order==k},'do_medfilt_x')
        disp('  ...medfilt_x')
        datatraces=medfilt_x(datatraces,t,params.numsamp_x,params.tstart_x);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('x-wise Median Filter')
        end
    end

    if strcmp(steps{order==k},'do_badTraceRemoval')
        disp('  ...bad trace removal')
        datatraces=traceInterpolation(datatraces,params.minfactor,params.maxfactor);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Bad Trace Removal')
        end
    end

    if strcmp(steps{order==k},'do_t0shift')
        disp('  ...t0shift')
        datatraces=t0shift(datatraces,params.t0s,abs(t(2)-t(1)));
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('t0 shift')
        end
    end

    if strcmp(steps{order==k},'do_spectralWhitening')
        disp('  ...spectralWhitening')
        for ch=1:max(info(3,:)) % for each channel
            datatraces(:,info(3,:)==ch)=spectralWhitening(datatraces(:,info(3,:)==ch),params.dt,params.fmin_sw,params.fmax_sw,params.alpha);
        end
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Spectral Whitening')
        end
    end

    if strcmp(steps{order==k},'do_topomig2d')
        disp('  ...Topomigration/correction')
        % check for constant trace spacing:
        xtemp=[0:params.dx:params.dx*(length(datatraces(1,info(3,:)==1))-1)];
        dx=xtemp(2)-xtemp(1); %[m] Trace distance
        if round(dx*1000)~=round(mean(diff(xtemp))*1000)
            mode.Interpreter='tex';
            mode.WindowStyle='non-modal';
            msgbox('\fontsize{15}Constant trace spacing is neccessary before topomigration! Please add "constant trace distance" (and optionally "trace interpolation") before!','Error','warn',mode);
        else
            for ch=1:max(info(3,:)) % for each channel
                [datatemp(:,info(3,:)==ch),zmig]=topomig2d_varV(datatraces(:,info(3,:)==ch),[0:params.dx:params.dx*(length(datatraces(1,info(3,:)==ch))-1)],t,info(6,info(3,:)==ch),mean(params.v),params.aperture_t,params.flagtopo,0,params.minz,params.maxz);
            end
            datatraces=single(datatemp);
            ns=length(datatraces(:,1));
            clear datatemp;
            if exist('fh','var')
                subplot(4,3,k+1)
                hold on
                plot(zmig,datatraces(:,b))
                grid on
                xlabel('z [m]')
                title('Topomigration/correction')
            end
        end
    end

    if strcmp(steps{order==k},'do_migration2d_vz')
        disp('  ...Isochrone migration')
        tp=t;
        if length(tp(:,1))<length(tp(1,:))
            tp=tp';
        end
        for ch=1:max(info(3,:)) % for each channel
            % interpolate vgrid
            vgrid=repmat(interp1(params.tfile,params.vfile,tp),[1 length(datatraces(1,info(3,:)==ch))]);
            disp(['      channel #',int2str(ch)])
            [datatemp(:,info(3,:)==ch),zmig]=isochrone_mig_2d_varV(datatraces(:,info(3,:)==ch),[0:params.dx:params.dx*(length(datatraces(1,info(3,:)==ch))-1)],tp,vgrid,params.aperture_m,0);
        end
        datatraces=single(datatemp);
        ns=length(datatraces(:,1));
        clear datatemp;
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(zmig,datatraces(:,b))
            grid on
            xlabel('z [m]')
            title('Isochrone migration')
        end
    end

    if strcmp(steps{order==k},'do_khighpass')
        disp('  ...k_highpass')
        datatraces=k_highpass(datatraces,params.dx,params.kcutoff);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('k-highpass')
        end
    end

    if strcmp(steps{order==k},'do_amplitudeOffset')
        disp('  ...Removal of amplitude offset')
        datatraces=DCremoval(datatraces,t,params.tstart,params.tend);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Amplitude offset')
        end
    end

    if strcmp(steps{order==k},'do_cutTWT')
        disp('  ...Cut time axis')
        [datatraces,t,ns]=cutTWT(datatraces,t,params.tmax);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Cut TWT')
        end
    end

    if strcmp(steps{order==k},'do_removeMeanMedianTrace')
        disp('  ...Removal of mean/median trace')
        if params.meanmedian==1
            mm='mean';
        else
            mm='median';
        end
        for ch=1:max(info(3,:))
            datatraces(:,info(3,:)==ch)=removeHorizontalLines(datatraces(:,info(3,:)==ch),mm,params.numtraces);
        end
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Mean/median trace removal')
        end
    end


    if strcmp(steps{order==k},'do_t0CorrectionThreshold')
        disp('  ...t0 correction with threshold')
        [datatraces]=t0corr_thresh(datatraces,params.threshold);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('t0-threshold')
        end
    end


    if strcmp(steps{order==k},'do_t0CorrectionReferencetrace')
        disp('  ...t0 correction with reference trace')
        datatraces=t0correction(datatraces,params.reftrace,params.t0,params.dt);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('t0-reference trace')
        end
    end

    if strcmp(steps{order==k},'do_applygain')
        disp('  ...Apply gain')
        datatraces=applygain(datatraces,[params.g1 params.g2 params.g3 params.g4 params.g5]);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Gain')
        end
    end

    if strcmp(steps{order==k},'do_crossprofileshift')
        disp('  ...CrossProfileShift')
        tempch2{1}=datatraces;
        tempch=channelshift(tempch2,params.maxshift_cps,params.shiftsamples_cps(params.shiftsamples_cps(:,1)==info(1,1),2));
        datatraces=single(tempch2{1});
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('z [m]')
            title('CrossProfileShift')
        end
    end

    if strcmp(steps{order==k},'do_channelshift')
        disp('  ...Channelshift')
        % sort traces into channels
        tempch=cell(max(info(3,:)),1);
        for ch=1:length(tempch)
            tempch{ch}=datatraces(:,info(3,:)==ch);
        end
        if params.oneall==1 % apply shiftsamples from reference profile
            tempch=channelshift(tempch,params.maxshift,params.shiftsamples);
        elseif params.oneall==2 % apply shift individually for each profile
            tempch=channelshift(tempch,params.maxshift);
        end
        % sort traces back
        for ch=1:length(tempch)
            datatraces(:,info(3,:)==ch)=tempch{ch};
        end
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Channelshift')
        end
    end


    if strcmp(steps{order==k},'do_sphericalDivergence')
        disp('  ...correct for spherical divergence')
        datatraces=sphericalDivergence(datatraces,t);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Spherical divergence')
        end
    end

    if strcmp(steps{order==k},'do_attenuationCorrection')
        disp('  ...correct for attenuation')
        datatraces=attenuationcorrection(datatraces,t,params.sigma,params.eps);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Attenuation')
        end
    end


    if strcmp(steps{order==k},'do_bandpass')
        disp('  ...Bandpass')
        datatraces=bandpass_gpr(datatraces,params.dt,params.fstart,params.fend);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Bandpass')
        end
    end

    if strcmp(steps{order==k},'do_normalization')
        disp('  ...Normalization')
        datatraces=normalize2d(datatraces,params.qclip);
        if exist('fh','var')
            subplot(4,3,k+1)
            hold on
            plot(t,datatraces(:,b))
            grid on
            xlabel('t [ns]')
            title('Normalization')
        end
    end
end
if strcmp(class(datatraces),'double')
    datatraces=single(datatraces);
end
end


function parsave(name,var)
var_name=genvarname(inputname(2));
eval([var_name '=var;']);
save(name,var_name,'-v7.3');
end

function parprintfigure(name,fh)
print(fh,name,'-djpeg');
end
