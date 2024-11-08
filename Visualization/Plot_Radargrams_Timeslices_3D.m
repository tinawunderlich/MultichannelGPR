clear all
close all
clc


%% Load radargrams and timeslices and plot them in 3D
%
% Dr. Tina Wunderlich, CAU Kiel, 2022, tina.wunderlich@ifg.uni-kiel.de
%
% Radargrams and timeslices need to be in the MultichannelGPR-format.


% time or depth? (Take care that both radargrams and timeslices are in the
% same domain!)
timedepth=1; % =1: time, =2: depth

% Radargrams:
rad_num=[]; % give number of radagrams to plot (if =[] all radargrams are used, be careful with large datasets!)

% Timeslices:
tsl_num=[]; % give number of timeslices to plot (or leave empty if no tsl)

% plot settings:
colorclip=5; % 0 is colorscale from min(data) to max(data), 1 is 1% clip value, 2 is 2% clip value and 3 is 3% clip value, ... (will not be saved, for plotting only!)
dxticks=[]; % give the increment for ticks along x-axis
dyticks=[]; % give the increment for ticks along y-axis
dzticks=[5]; % give the increment for ticks along t- or z-axis
aspectratio_t=5; % aspect ratio for time plots
aspectratio_z=1; % aspect ratio for depth plots


%% ------------- DO NOT CHANGE BELOW THIS LINE ----------------------------

% get folder name for radargrams
if ~ispc; menu('Choose folder with radargrams','OK'); end
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
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
        
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
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end
    
    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

% get folder name for timeslices
if ~isempty(tsl_num)
    if ~ispc; menu('Choose folder with timeslices','OK'); end
    if ispc
        if exist('temptsl.temp') % read last opened folder from temp.temp
            fid=fopen('temptsl.temp','r');
            fn=textscan(fid,'%s');
            fclose(fid);
            if ~isempty(fn{1})
                pfad_tsl=uigetdir(fn{1}{1},'Choose folder with timeslices');
            else
                pfad_tsl=uigetdir([],'Choose folder with timeslices');
            end
            fileattrib('temptsl.temp','-h');
            fid=fopen('temptsl.temp','wt');
            fprintf(fid,'%s',pfad_tsl);
            fclose(fid);
            fileattrib('temptsl.temp','+h');
        else
            pfad_tsl=uigetdir([],'Choose folder with timeslices'); % path to radargram-folder
            
            fid=fopen('temptsl.temp','wt');
            fprintf(fid,'%s',pfad_tsl);
            fclose(fid);
            fileattrib('temptsl.temp','+h');
        end
    else
        if exist('.temptsl.temp') % read last opened folder from temp.temp
            fid=fopen('.temptsl.temp','r');
            fn=textscan(fid,'%s');
            fclose(fid);
            if ~isempty(fn{1})
                pfad_tsl=uigetdir(fn{1}{1},'Choose folder with timeslices');
            else
                pfad_tsl=uigetdir([],'Choose folder with timeslices');
            end
        else
            pfad_tsl=uigetdir([],'Choose folder with timeslices'); % path to radargram-folder
        end
        
        fid=fopen('.temptsl.temp','wt');
        fprintf(fid,'%s',pfad_tsl);
        fclose(fid);
    end
end

% temporarily set path to required scripts
oldpath=path;
addpath('../Subfunctions/','../GUIs/');


% load radargrams:
load(fullfile(pfad,'radargrams.mat'));
load(fullfile(pfad,'global_coords.mat'));
load(fullfile(pfad,'t.mat'));
load(fullfile(pfad,'x.mat'));
dt=t(2)-t(1);
ns=length(t);

% load tsl:
if ~isempty(tsl_num)
    load(fullfile(pfad_tsl,'t_startende.mat'));
    if exist(fullfile(pfad_tsl,'coordtrans.mat'),'file')
        load(fullfile(pfad_tsl,'coordtrans.mat'));
    else
        coordtrans=[1 1 1 1; 2 2 2 2];
    end
    if exist(fullfile(pfad_tsl,'tsl_interp.mat'),'file')
        load(fullfile(pfad_tsl,'tsl_interp.mat'));
        load(fullfile(pfad_tsl,'topo_interp.mat'));
        load(fullfile(pfad_tsl,'xgrid_interp.mat'));
        load(fullfile(pfad_tsl,'ygrid_interp.mat'));
        if exist(fullfile(pfad_tsl,'depth.mat'),'file')
            temp=load(fullfile(pfad_tsl,'depth.mat'));
            dsl=1;
            maxElev=temp.maxElevation;
        else
            maxElev=[];
            dsl=0;
        end
        % rename for further internal use:
        tsl=tsl_interp;
        topo=topo_interp;
        xgrid=xgrid_interp;
        ygrid=ygrid_interp;
    elseif exist(fullfile(pfad_tsl,'tsl.mat'),'file')
        load(fullfile(pfad_tsl,'tsl.mat'));
        load(fullfile(pfad_tsl,'topo.mat'));
        load(fullfile(pfad_tsl,'xgrid.mat'));
        load(fullfile(pfad_tsl,'ygrid.mat'));
        if exist(fullfile(pfad_tsl,'depth.mat'),'file')
            temp=load(fullfile(pfad_tsl,'depth.mat'));
            dsl=1;
            maxElev=temp.maxElevation;
        else
            maxElev=[];
            dsl=0;
        end
    end
end

if isempty(rad_num)
    rad_num=1:length(radargrams);
end

%% plotting:
f=figure;
hold on
cdata=[];
for i=1:length(rad_num) % plot radargrams
    X=repmat(global_coords{rad_num(i)}(:,1)',[length(t),1]);
    Y=repmat(global_coords{rad_num(i)}(:,2)',[length(t),1]);
    Z=repmat(t',[1 length(X(1,:))]);
    C=normalize2d(radargrams{rad_num(i)});
    surf(X,Y,Z,C);
    cdata=[cdata; C(:)];
end

if ~isempty(tsl_num)
    for i=1:length(tsl_num) % plot timeslices
        if dsl==1
            Z=zeros(size(xgrid))+maxElev-mean(t_startende(tsl_num(i),:));
        else
            Z=zeros(size(xgrid))+mean(t_startende(tsl_num(i),:));
        end
        C=normalize2d(tsl{tsl_num(i)});
        % apply coordtrans
        temp=[xgrid(:) ygrid(:) C(:)]; % local coord system
        neu=helmert(temp,coordtrans(:,1:2),coordtrans(:,3:4)); % global coord system
        % check if same coordinate system between tsl and radargrams
        if abs(neu(1,1)-X(1,1))>abs(temp(1,1)-X(1,1))
            % local coordinates
            neu=temp; % use temp...
        end
        xgrid2=reshape(neu(:,1),size(xgrid));
        ygrid2=reshape(neu(:,2),size(xgrid));
       
        % plot
        surf(xgrid2,ygrid2,Z,C);
        cdata=[cdata; C(:)];
    end
end
shading interp
colormap(flipud(gray))
grid on
if timedepth==1 % time
    set(gca,'ZDir','reverse','DataAspectratio',[1 1 aspectratio_t]);
    zlabel('t [ns]')
else
    set(gca,'ZDir','normal','DataAspectratio',[1 1 aspectratio_z]);
    zlabel('z [m]')
end
if ~isempty(dxticks)
    xlim=get(gca,'XLim');
    set(gca,'XTick',[xlim(1):dxticks:xlim(2)])
end
if ~isempty(dyticks)
    ylim=get(gca,'YLim');
    set(gca,'YTick',[ylim(1):dyticks:ylim(2)])
end
if ~isempty(dzticks)
    zlim=get(gca,'ZLim');
    set(gca,'ZTick',[zlim(1):dzticks:zlim(2)])
end
xticks=get(gca,'XTick');
yticks=get(gca,'YTick');
for i=1:length(xticks)
    xticklabels{i}=sprintf('%d',xticks(i));
end
for i=1:length(yticks)
    yticklabels{i}=sprintf('%d',yticks(i));
end
xlabel('x [m]')
ylabel('y [m]')
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels)

if colorclip~=0
    % determine color limits for plotting:
    coldata=sort(unique(cdata));
    coldata(isnan(coldata))=[]; % delete nans
    cmin=coldata(round(length(coldata)/100*colorclip));
    cmax=coldata(end-round(length(coldata)/100*colorclip));
    set(gca,'CLim',[cmin cmax])
else
    set(gca,'ClimMode','auto');
end


% restore original path
path(oldpath);

%%
function [traces]=normalize2d(traces,qclip)

% Normalize traces of 2D-matrix [traces]=normalize2d(traces,qclip)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: 2d matrix with traces along columns
% qclip: cumulative probability in interval [0,1], default is 0.98
%
% Output:
% traces: normalized traces in matrix

if nargin==1
    qclip=0.98;
end

ns=length(traces(:,1));

% matrix with quantiles (constant for each trace):
%referscale=repmat(quantile(abs(traces),qclip,1),[ns 1]);
referscale=repmat(quanti(abs(traces(~isnan(traces))),qclip),[ns 1]);
% scale data:
traces=traces./referscale;
end