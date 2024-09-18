clear all
close all
clc


% Read radargrams (radargrams.mat, global_coords.mat, x.mat and t.mat) from
% make_Radargrams or Bins2radargrams or other radargrams in same format
% (e.g. from DZT_Convert) and do processing of them -> save in new folder
% "processed"
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% requires content of Processing, Subfunctions and Migration folder



numbers=[]; % give numbers of processed radargrams or leave empty =[] for all

% Processing options:
settings='settings_khighpass.txt'; % give filename of settings-file in Radargram-folder (also give a filename, if you want to create a default file!)
plotflag=0; % =1 plot and show all radargrams during processing, =0 do not plot

% Plotting options
colorclip=3; % 0 is colorscale from min(data) to max(data), 1 is 1% clip value, 2 is 2% clip value and 3 is 3% clip value, ... (will not be saved, for plotting only!)
aspectratio_t=7;    % for time plots: give aspectratio for y-axis. If you want to plot over whole screen, set =0
aspectratio_z=1;    % for depth plots: give aspectratio for y-axis. If you want to plot over whole screen, set =0

% -------------------------------------------------------------------------
% Do not change the following part!

% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with radargrams','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
        fileattrib('radtemp.temp','-h');
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    end
else
    if exist('.radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('.radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad_rad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad_rad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end

    fid=fopen('.radtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end


% temporarily set path to required scripts
oldpath=path;
addpath('../Processing/','../Subfunctions/','../Migration/');


%%% Read data
disp('Reading data...')
temp=load(fullfile(pfad_rad,'global_coords.mat'));
gc=temp.global_coords; % global coordinates of starting end ending point
temp=load(fullfile(pfad_rad,'x.mat'));
xx=temp.x;   % profile coordinates
temp=load(fullfile(pfad_rad,'radargrams.mat'));
data=temp.radargrams;   % radargrams
temp=load(fullfile(pfad_rad,'t.mat'));
t=temp.t;   % time vector
dt=t(2)-t(1);

if isempty(numbers)
    numbers=1:length(data); % all radargrams
end

%-------------------------------------------------------------------------
disp('Reading processing settings:')
if ~exist(fullfile(pfad_rad,settings),'file') % if no settings-file is found: create default file
    disp(['No file ',settings,' found. Creating default file ',settings,'.']);
    
    fid=fopen(fullfile(pfad_rad,settings),'wt');
    fprintf(fid,'do_DCremoval 1\n\n');
    fprintf(fid,'do_bandpass 2\nfstart 100\nfend 600\n\n');
    fprintf(fid,'do_medfilt 0\nnumsamp 3\n\n');
    fprintf(fid,'do_medfilt_x 0\nnumsamp_x 3\ntstart_x 0\n\n');
    fprintf(fid,'do_constTraceDist 4\ndist 0.02\n\n');
    fprintf(fid,'do_reduceNumberOfSamples 0\nn 2\n\n');
    fprintf(fid,'do_cutTWT 3\ncutT 80\n\n');
    fprintf(fid,'do_t0correction 0\nt0 5\n\n');
    fprintf(fid,'do_t0shift 0\nt0s 5\n\n');
    fprintf(fid,'do_t0Threshold 0\nthreshold 1000\n\n');
    fprintf(fid,'do_sphericalDivergence 0\n\n');
    fprintf(fid,'do_attenuationCorrection 0\nsigma 0.002\neps 15\n\n');
    fprintf(fid,'do_normalization 0\nqclip 0.98\n\n');
    fprintf(fid,'do_removeMeanTrace 0\nmeanMedian 1\nnumberOfTraces 0\n\n');
    fprintf(fid,'do_makeAmpSpec 0\n\n');
    fprintf(fid,'do_turnProfiles 0\n\n');
    fprintf(fid,'do_exchange_x_y 0\n\n');
    fprintf(fid,'do_kHighpass 0\nkcutoff 0.1\n\n');
    fprintf(fid,'do_interpolation 0\ngap 3\n\n');
    fprintf(fid,'do_helmertTransformation 0\ncoordsfile coords.txt\n\n');
    fprintf(fid,'do_migration 0\nvfile_mig vgrid.mat\naperture_mig 30\nverbose_mig 1\n\n');
    fprintf(fid,'do_topomigration 0\ntopofile topo.mat\nvfile_topomig vgrid.mat\naperture_topomig 30\nflag 1\nverbose_topomig 1\nzmin -2\nzmax 0\n\n');
    fprintf(fid,'do_applyGain 5\ng -20 0 10 15 20\n');
    fclose(fid);
    
    disp(['Edit ',settings,' and start script again.'])
    return;
else
    fid=fopen(fullfile(pfad_rad,settings),'r');
    temp=textscan(fid,'%s%s');
    fclose(fid);
    % get order of processing steps:
    steps=[{'do_DCremoval'} {'do_medfilt'} {'do_medfilt_x'} {'do_helmertTransformation'} {'do_turnProfiles'} {'do_t0Threshold'} {'do_exchange_x_y'} {'do_migration'} {'do_reduceNumberOfSamples'} {'do_constTraceDist'} {'do_topomigration'} {'do_t0correction'} {'do_t0shift'} {'do_cutTWT'} {'do_makeAmpSpec'} {'do_sphericalDivergence'} {'do_attenuationCorrection'} {'do_bandpass'} {'do_normalization'} {'do_removeMeanTrace'} {'do_kHighpass'} {'do_applyGain'} {'do_interpolation'}];
    for i=1:length(steps)
        for j=1:length(temp{1})
            if strcmp(temp{1}{j},steps{i})
                order(i)=str2num(temp{2}{j});
            end
        end
    end
    % get settings for parameters:
    for i=1:length(temp{1})
        if strcmp(temp{1}(i),'sigma')
            sigma=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'eps')
            eps=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'numsamp')
            numsamp=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'numsamp_x')
            numsamp_x=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'tstart_x')
            tstart_x=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'t0')
            t0=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'t0s')
            t0s=str2num(temp{2}{i});
        elseif strcmp(temp{2}(i),'threshold')
            threshold=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'dist')
            dist=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'cutT')
            cutT=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'n')
            n=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'fstart')
            fstart=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'fend')
            fend=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'qclip')
            qclip=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'coordsfile')
            coordsfile=temp{2}{i};
        elseif strcmp(temp{1}(i),'vfile_mig')
            vfilem=temp{2}{i};
        elseif strcmp(temp{1}(i),'aperture_mig')
            aperturem=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'verbose_mig')
            verbosem=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'vfile_topomig')
            vfile=temp{2}{i};
        elseif strcmp(temp{1}(i),'topofile')
            topofile=temp{2}{i};
        elseif strcmp(temp{1}(i),'aperture_topomig')
            aperture=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'flag')
            flag=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'verbose_topomig')
            verbose=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'zmin')
            zmin=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'zmax')
            zmax=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'meanMedian')
            meanMedian=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'numberOfTraces')
            numtraces=str2num(temp{2}{i});
        elseif strcmp(temp{1}(i),'kcutoff')
            kcutoff=str2num(temp{2}{i});
        elseif strcmp(temp{1}{i},'g')
            g=[str2num(temp{2}{i}) str2num(temp{1}{i+1}) str2num(temp{2}{i+1}) str2num(temp{1}{i+2}) str2num(temp{2}{i+2})];
        elseif strcmp(temp{1}(i),'gap')
            gap=str2num(temp{2}{i});
        end
    end
    % display order of steps
    for i=1:length(order(order>0))
        disp([int2str(i),'. ',steps{order==i}])
        if strcmp('do_attenuationCorrection',steps{order==i})
            disp(['  sigma = ',num2str(sigma),' S/m'])
            disp(['  eps = ',num2str(eps)])
        elseif strcmp('do_bandpass',steps{order==i})
            disp(['  fstart = ',num2str(fstart),' MHz'])
            disp(['  fend = ',num2str(fend),' MHz'])
        elseif strcmp('do_constTraceDist',steps{order==i})
            disp(['  dist = ',num2str(dist),' m'])
        elseif strcmp('do_medfilt',steps{order==i})
            disp(['  numsamp = ',num2str(numsamp)])
        elseif strcmp('do_medfilt_x',steps{order==i})
            disp(['  numsamp_x = ',num2str(numsamp_x)])
            disp(['  tstart_x = ',num2str(tstart_x),' ns'])
        elseif strcmp('do_normalization',steps{order==i})
            disp(['  qclip = ',num2str(qclip)])
        elseif strcmp('do_cutTWT',steps{order==i})
            disp(['  cutT = ',num2str(cutT),' ns'])
        elseif strcmp('do_reduceNumberOfSamples',steps{order==i})
            disp(['  n = ',num2str(n)])
        elseif strcmp('do_t0correction',steps{order==i})
            disp(['  t0 = ',num2str(t0),' ns'])
        elseif strcmp('do_t0shift',steps{order==i})
            disp(['  t0 = ',num2str(t0s),' ns'])
        elseif strcmp('do_t0Threshold',steps{order==i})
            disp(['  threshold = ',num2str(threshold)])
        elseif strcmp('do_helmertTransformation',steps{order==i})
            disp(['  coordsfile = ',coordsfile])
        elseif strcmp('do_migration',steps{order==i})
            disp(['  vfile_mig = ',vfilem])
            disp(['  aperture_mig = ',num2str(aperturem),' degree'])
            disp(['  verbose_mig = ',num2str(verbosem)])
        elseif strcmp('do_topomigration',steps{order==i})
            disp(['  topofile = ',topofile])
            disp(['  vfile_topomig = ',vfile])
            disp(['  aperture_topomig = ',num2str(aperture),' degree'])
            if ~isempty(zmin)
                disp(['  zmin = ',num2str(zmin),' m'])
            else
                disp('zmin will be determined automatically.')
            end
            if ~isempty(zmax)
                disp(['  zmax = ',num2str(zmax),' m'])
            else
                disp('zmax will be determined automatically.')
            end
            disp(['  flag = ',num2str(flag)])
            disp(['  verbose_topomig = ',num2str(verbose)])
        elseif strcmp('do_removeMeanTrace',steps{order==i})
            if meanMedian==1
                meanMedian='mean';
            else
                meanMedian='median';
            end
            disp(['  meanMedian = ',meanMedian])
            if ~exist('numtraces','var')
                numtraces=0;
            end
            disp(['  numberOfTraces = ',numtraces])
        elseif strcmp('do_kHighpass',steps{order==i})
            disp(['  kcutoff = ',num2str(kcutoff)])
        elseif strcmp('do_applyGain',steps{order==i})
            disp(['  gain = ',num2str(g)])
        elseif strcmp('do_interpolation',steps{order==i})
            disp(['  gap = ',num2str(gap)])
        end
    end
end
disp('--------------------------------------')

%--------------------------------------------------------------------------
h=waitbar(0,'Processing, please wait!');
numsteps=length(order(order>0));    % number of steps per radargram
anz=length(numbers);
radargrams=cell(anz,1);
nn=0;
for kk=numbers % loop over radargrams
    disp(['   -> Profile ',int2str(kk),' of ',int2str(anz)])

    if ~isempty(data{kk}) && any(~isnan(data{kk}(:)))
        
        nn=nn+1;
        datatraces=data{kk};  % raw data of current radargram
        
        tcut_flag=0; % set to 0 for each new profile
        amp_flag=0; % if amplitude spectrum
        z_flag=0; % if depth migrated data
        
        if plotflag==1
            % Plot raw data
            figure(kk)
            subplot(2,1,1)
            imagesc(xx{kk},t,datatraces)
            grid on
            xlabel('x [m]')
            ylabel('t [ns]')
            title(['Original data of radargram ',int2str(kk)])
            colormap(flipud(gray));
            if aspectratio_t~=0
                set(gca,'Dataaspectratio',[1 aspectratio_t 1])
            end
            if colorclip~=0
                % determine color limits for plotting:
                coldata=sort(unique(datatraces));
                coldata(isnan(coldata))=[]; % delete nans
                cmin=coldata(round(length(coldata)/100*colorclip));
                cmax=coldata(end-round(length(coldata)/100*colorclip));
                set(gca,'CLim',[cmin cmax])
            else
                set(gca,'ClimMode','auto');
            end
        end
        
        
        for k=1:length(order(order>0))  % for all processing steps in right order
            
            %%% trace-wise median filter
            if strcmp(steps{order==k},'do_medfilt')
                [datatraces]=medfilt(datatraces,numsamp);
            end

            %%% x-wise median filter
            if strcmp(steps{order==k},'do_medfilt_x')
                [datatraces]=medfilt_x(datatraces,t,numsamp_x,tstart_x);
            end
            
            %%% Make constant trace distance
            if strcmp(steps{order==k},'do_constTraceDist')
                [datatraces,xx{kk},gc{kk}]=constTraceDist(datatraces,dist,xx{kk},gc{kk});
            end
            
            %%% Turn Profiles
            if strcmp(steps{order==k},'do_turnProfiles')
                [datatraces,xx{kk},gc{kk}]=turnProfiles(datatraces,xx{kk},gc{kk});
            end

            %%% Exchange x and y
            if strcmp(steps{order==k},'do_exchange_x_y')
                gc{kk}=exchange_x_y(gc{kk});
            end
            
            %%% Reduce Number of samples
            if strcmp(steps{order==k},'do_reduceNumberOfSamples')
                if tcut_flag==0
                    [datatraces,tnew]=reduceNumberOfSamples(datatraces,t,n);
                    nsnew=length(tnew);
                    tcut_flag=1;
                else
                    [datatraces,tnew]=reduceNumberOfSamples(datatraces,tnew,n);
                    nsnew=length(tnew);
                end
                dt=tnew(2)-tnew(1);
            end
                
            %%% Topomigration with 2d-v-model
            if strcmp(steps{order==k},'do_topomigration')
                % read required files and check format:
                top=load(fullfile(pfad_rad,topofile));
                name=fieldnames(top);
                eval(['topo=top.' name{1} ';']); % make variable topo from structure
                if length(xx{kk})~=length(topo{kk})
                    disp('Size of topography vector does not match number of traces! Please check and start again.')
                    return;
                end
                vf=load(fullfile(pfad_rad,vfile));
                name=fieldnames(vf);
                eval(['v=vf.' name{1} ';']); % make variable v from structure
                if length(v{kk})>1 && any(size(v{kk})~=size(datatraces))
                    disp('Size of vgrid does not match size of radargram! Trying to fix this...')
                    if length(v{kk}(1,:))~=length(datatraces(1,:))
                        disp('   Number of traces in vgrid and radargram does not match. Please check and start again.')
                        return;
                    elseif length(v{kk}(:,1))<length(datatraces(:,1))
                        v{kk}=[v{kk}; repmat(v{kk}(end,:),[length(datatraces(:,1))-length(v{kk}(:,1)),1])];
                        disp('   Extrapolated vgrid for all times in radargram. Problem solved.')
                    elseif length(v{kk}(:,1))>length(datatraces(:,1))
                        v{kk}=v{kk}(1:length(datatraces(:,1)));
                        disp('   Shortened vgrid corresponding to samples in radargram. Problem solved.')
                    end
                end
                
                if isempty(zmin) % if no zmin/zmax given -> determine one for all profiles
                    zm=[];
                    vm=[];
                    zmi=[];
                    for kk2=numbers
                        zm=[zm; max(topo{kk2})];
                        zmi=[zmi; min(topo{kk2})];
                        vm=[vm; max(v{kk2}(:))];
                    end
                    zmax=ceil(max(zm));
                    zmin=floor(min(zmi)-max(t)/2*max(vm));
                end
                    
                if tcut_flag==0
                    [datatraces,z]=topomig2d_varV(datatraces,xx{kk},t,topo{kk},v{kk},aperture,flag,verbose,zmin,zmax); % Output t is depth z in m!
                else % if t has been cut
                    [datatraces,znew]=topomig2d_varV(datatraces,xx{kk},tnew,topo{kk},v{kk},aperture,flag,verbose,zmin,zmax); % Output tnew is depth z in m!
                end
                z_flag=1;
            end
            
            %%% Isochrone migration with 2d-v-model
            if strcmp(steps{order==k},'do_migration')
                % read required files and check format:
                vf=load(fullfile(pfad_rad,vfilem));
                name=fieldnames(vf);
                eval(['v=vf.' name{1} ';']); % make variable topo from structure
                if length(v{kk})>1 && any(size(v{kk})~=size(datatraces))
                    disp('Size of vgrid does not match size of radargram! Trying to fix this...')
                    if length(v{kk}(1,:))~=length(datatraces(1,:))
                        disp('   Number of traces in vgrid and radargram does not match. Please check and start again.')
                        return;
                    elseif length(v{kk}(:,1))<length(datatraces(:,1))
                        v{kk}=[v{kk}; repmat(v{kk}(end,:),[length(datatraces(:,1))-length(v{kk}(:,1)),1])];
                        disp('   Extrapolated vgrid for all times in radargram. Problem solved.')
                    elseif length(v{kk}(:,1))>length(datatraces(:,1))
                        v{kk}=v{kk}(1:length(datatraces(:,1)));
                        disp('   Shortened vgrid corresponding to samples in radargram. Problem solved.')
                    end
                end
                if tcut_flag==0
                    [datatraces,z]=isochrone_mig_2d_varV(datatraces,xx{kk},t,v{kk},aperturem,verbosem); % Output t is depth z in m!
                else % if t has been cut
                    [datatraces,znew]=isochrone_mig_2d_varV(datatraces,xx{kk},tnew,v{kk},aperturem,verbosem); % Output tnew is depth z in m!
                end
                z_flag=1;
            end
            
            %%% Helmert Transformation
            if strcmp(steps{order==k},'do_helmertTransformation')
                coords=load(fullfile(pfad_rad,coordsfile));
                gc{kk}=helmert(gc{kk},coords(:,1:2),coords(:,3:4));
            end
            
            %%% t0 correction (correlation with first trace as reference)
            if strcmp(steps{order==k},'do_t0correction')
                [datatraces]=t0correction(datatraces,datatraces(:,1),t0,dt);
            end
            
            %%% t0 shift
            if strcmp(steps{order==k},'do_t0shift')
                [datatraces]=t0shift(datatraces,t0s,dt);
            end

            %%% t0 threshold
            if strcmp(steps{order==k},'do_t0Threshold')
                [datatraces]=t0corr_thresh(datatraces,threshold);
            end
            
            %%% DCremoval
            if strcmp(steps{order==k},'do_DCremoval')
                if tcut_flag==0
                    [datatraces]=DCremoval(datatraces,t);
                else
                    [datatraces]=DCremoval(datatraces,tnew);
                end
            end
            
            %%% cut TWT
            if strcmp(steps{order==k},'do_cutTWT')
                if tcut_flag==1 % if reduceNumberofSamples before
                    [datatraces,tnew,nsnew]=cutTWT(datatraces,tnew,cutT);
                elseif tcut_flag==0 && z_flag==0
                    [datatraces,tnew,nsnew]=cutTWT(datatraces,t,cutT);
                elseif z_flag==1 % if depth migration before
                    % make new depth vector starting from top of radargram
                    % pointing downwards (compatible with time direction)
                    z2=abs(z-max(z));
                    cutT2=abs(cutT-max(z));
                    [datatraces,znew2,nsnew]=cutTWT(datatraces,z2,cutT2);
                    % converting back to original z vector:
                    znew=-znew2+max(z);
                end
                tcut_flag=1;
            end
            
            %%% makeAmpSpec
            if strcmp(steps{order==k},'do_makeAmpSpec')
                if tcut_flag==0
                    [t,datatraces]=makeAmpspec(t,datatraces); % Output t is frequency and output datatraces is spectrum!
                else
                    [tnew,datatraces]=makeAmpspec(tnew,datatraces); % Output tnew is frequency and output datatraces is spectrum!
                end
                amp_flag=1;
            end
            
            %%% correct for spherical divergence
            if strcmp(steps{order==k},'do_sphericalDivergence')
                if tcut_flag==0
                    [datatraces]=sphericalDivergence(datatraces,t);
                else
                    [datatraces]=sphericalDivergence(datatraces,tnew);
                end
            end
            
            %%% attenuation correction
            if strcmp(steps{order==k},'do_attenuationCorrection')
                if tcut_flag==0
                    [datatraces]=attenuationcorrection(datatraces,t,sigma,eps);
                else
                    [datatraces]=attenuationcorrection(datatraces,tnew,sigma,eps);
                end
            end
            
            %%% apply bandpass
            if strcmp(steps{order==k},'do_bandpass')
                [datatraces]=bandpass_gpr(datatraces,dt,fstart,fend);
            end
            
            %%% trace normalization
            if strcmp(steps{order==k},'do_normalization')
                [datatraces]=normalize2d(datatraces,qclip);
            end
            
            %%% remove horizontal Lines (subtract mean trace)
            if strcmp(steps{order==k},'do_removeMeanTrace')
                [datatraces]=removeHorizontalLines(datatraces,meanMedian,numtraces);
            end
            
            %%% k-highpass for horizontal lines removal
            if strcmp(steps{order==k},'do_kHighpass')
                if length(datatraces(1,:))<24
                    disp('No k_highpass can be applied due to too less traces! Needs to have more than 24 traces. Continuing processing.');
                else
                    [datatraces]=k_highpass(datatraces,mean(diff(xx{kk})),kcutoff);
                end
            end
            
            %%% apply gain
            if strcmp(steps{order==k},'do_applyGain')
                [datatraces]=applygain(datatraces,g);
            end
            
            %%% interpolate gaps
            if strcmp(steps{order==k},'do_interpolation')
                [datatraces]=interpolation(datatraces,gap);
            end
            
            waitbar(((kk-1)*numsteps+k)/(numsteps*length(numbers)),h);
        end
       
        
        % set new coordinates (if changed during processing...)
        global_coords{nn}=gc{kk};
        x{nn}=xx{kk};
        
        if plotflag==1
            % Plot processed data
            figure(kk)
            subplot(2,1,2)
            if amp_flag==0
                % plot processed radargram
                if z_flag==0
                    if tcut_flag==0
                        imagesc(x{nn},t,datatraces)
                    else
                        imagesc(x{nn},tnew,datatraces)
                    end
                else
                    if tcut_flag==0
                        imagesc(x{nn},z,datatraces)
                    else
                        imagesc(x{nn},znew,datatraces)
                    end
                end
                grid on
                xlabel('x [m]')
                ylabel('t [ns]')
                if z_flag==1
                    xlabel('x [m]')
                    ylabel('z [m]')
                    axis xy
                end
                title(['Processed data of radargram ',int2str(kk)])
                colormap(flipud(gray));
                if z_flag==0
                    if aspectratio_t~=0
                        set(gca,'Dataaspectratio',[1 aspectratio_t 1])
                    end
                else
                    if aspectratio_z~=0
                        set(gca,'Dataaspectratio',[1 aspectratio_z 1])
                    end
                end
                if colorclip~=0
                    % determine color limits for plotting:
                    coldata=sort(unique(datatraces));
                    coldata(isnan(coldata))=[]; % delete nans
                    cmin=coldata(round(length(coldata)/100*colorclip));
                    cmax=coldata(end-round(length(coldata)/100*colorclip));
                    set(gca,'CLim',[cmin cmax])
                else
                    set(gca,'ClimMode','auto');
                end
            else
                % plot amplitude spectrum
                if tcut_flag==0
                    plot(t,datatraces)
                    hold on
                    plot(t,mean(datatraces,2),'k','Linewidth',2)
                else
                    plot(tnew,datatraces)
                    hold on
                    plot(tnew,mean(datatraces,2),'k','Linewidth',2)
                end
                xlabel('f [MHz]')
                ylabel('Amplitude')
                grid on
                title(['Amplitude spectrum of processed radargram ',int2str(kk)])
            end
        end
        
        
        radargrams{nn}=datatraces;
    else
        radargrams{nn}=data{kk};
    end
end
close(h);
if z_flag==1
    if tcut_flag==0
        t=z;
    else
        tnew=znew;
    end
end
%--------------------------------------------------------------------------
disp('Saving data...')
% completely processed radargrams are stored in variable "radargrams"
% save processed radargrams:
if ~exist(fullfile(pfad_rad,'processed'),'dir')
    mkdir(fullfile(pfad_rad,'processed'));
end
if tcut_flag==1
    t=tnew; % set new t
end
% save processed radargrams
save(fullfile(pfad_rad,'processed','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(pfad_rad,'processed','x.mat'),'x','-v7.3');
save(fullfile(pfad_rad,'processed','global_coords.mat'),'global_coords','-v7.3');
save(fullfile(pfad_rad,'processed','t.mat'),'t','-v7.3'); % can be time or depth!
% copy settings
copyfile(fullfile(pfad_rad,settings),fullfile(pfad_rad,'processed','ProcessingReadMe.txt'));

disp('Done!')

% restore original path
path(oldpath);