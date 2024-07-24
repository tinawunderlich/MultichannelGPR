clear all
close all
clc

%% kk-Pieslicefilter for Timeslices/Depthslices with selection of different areas (using a shape-file)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires interpolated Timeslices in MultichannelGPR-format (created with
% make_Timeslices.m)
% (requires scripts in "Subfunctions")
%
% When asked, choose folder with interpolated timeslices (containing tsl_interp.mat, xgrid_interp.mat, ...)


% name of shapefile for different parts of the area (in same folder as the timeslices):
shpfile='kkFilterShape.shp'; % needs to be a shape-file with polygons of the different areas for filtering 
all_area=0; % if =1 the complete area will be used with one stripe noise direction and shpfile will be ignored (for using shpfile set to 0)

% give name for new folder with filtered results and settings (will be
% created inside the interpolated timeslices folder)
name='kkfilt';

use_old_settings=1; % if =1: use old settings (angles of stripe noise from folder "name"), =0: make new angles

apply_to_all_tsl=0; % =1: if you want to appyl it to all timeslices and save in MultichannelGPR-format

apply_to_tslnum=1; % =1: if you want to apply it to tslnum only and save comparison figure
tslnum=14; % if you want to make new angles (use_old_settings=0), choose number of timeslice for plotting and interactive picking of stripe noise

opening=20; % opening angle of pieslice filter in degree

% plotting options:
colorclip=3;

plotflag=1; % do you want to plot original spectrum and filtered  spectrum (only for tslnum)? yes=1, no=0


%% ------------- DO NOT CHANGE FROM HERE ON! ----------------------------- %%
warning('off');

% get folder name
if ispc
    if exist('temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            folder=uigetdir(fn{1}{1},'Choose interpolated timeslices folder');
        else
            folder=uigetdir([],'Choose interpolated timeslices folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',folder);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        folder=uigetdir([],'Choose interpolated timeslices folder'); % path to timeslices-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',folder);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp','file') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            folder=uigetdir(fn{1}{1},'Choose interpolated timeslices folder');
        else
            folder=uigetdir([],'Choose interpolated timeslices folder');
        end
    else
        folder=uigetdir([],'Choose interpolated timeslices folder'); % path to timeslices-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',folder);
    fclose(fid);
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));

%%
disp('Loading data...')
load(fullfile(folder,'xgrid_interp.mat'));
xgrid=xgrid_interp;
load(fullfile(folder,'ygrid_interp.mat'));
ygrid=ygrid_interp;
load(fullfile(folder,'tsl_interp.mat'));
if exist(fullfile(folder,'coordtrans.mat'))
    load(fullfile(folder,'coordtrans.mat'));
else
    coordtrans=[1 1 1 1; 2 2 2 2];
end

dx=xgrid(1,2)-xgrid(1,1);

if ~exist(fullfile(folder,name))
    mkdir(fullfile(folder,name))
end

% prepare picking of angles:
data=tsl_interp{tslnum};
gesamtmaske=zeros(size(data));
gesamtmaske(~isnan(data))=1;
if use_old_settings==0
    if all_area==0
        disp('Reading shape-file and creating masks.')
        % read shape file
        S=shaperead(fullfile(folder,shpfile));
        % create masks for the areas
        for i=1:length(S)
            poly=helmert([S(i).X' S(i).Y'],coordtrans(:,3:4),coordtrans(:,1:2));
            mask{i}=scanline(data,xgrid(1,:),ygrid(:,1),poly,2); % sets poylgon to 1 and rest to 0
        end

        save(fullfile(folder,name,'mask.mat'),'mask','-v7.3');
    else
        mask{1}=ones(size(data));
        save(fullfile(folder,name,'mask.mat'),'mask','-v7.3');
    end
else
    disp('Loading masks.')
    load(fullfile(folder,name,'mask.mat'));
end

anz=length(mask); % number of different areas

if all_area==1
    % check if there is only one mask
    if anz~=1
        disp(['all_area==1, but there are ',int2str(anz),' masks. Please check settings and start again.'])
        return;
    end
end

if use_old_settings==0
    % pick stripe noise
    for i=1:anz % for every part
        part=data.*mask{i};
        part(isnan(part))=0;

        % calculate angle from picks
        figure
        imagesc(xgrid(1,:),ygrid(:,1),part)
        coldata=sort(unique(part(:)));
        coldata(isnan(coldata))=[]; % delete nans
        cmin=coldata(round(length(coldata)/100*colorclip));
        cmax=coldata(end-round(length(coldata)/100*colorclip));
        axis xy
        axis equal
        set(gca,'CLim',[cmin cmax],...
            'XLim',[min(xgrid(logical(mask{i})),[],'all') max(xgrid(logical(mask{i})),[],'all')],...
            'YLim',[min(ygrid(logical(mask{i})),[],'all') max(ygrid(logical(mask{i})),[],'all')])

        disp(['   Area ',int2str(i),': Pick two points in direction of stripe noise and press ENTER'])
        g=gline();
        pause;
        set(g,"Color","r")
        picks=[g.XData' g.YData'];
        close(gcf);

        angle(i)=-atand((picks(2,2)-picks(1,2))/(picks(2,1)-picks(1,1)));
        if abs(angle(i))>90
            angle(i)=90-atand((picks(2,2)-picks(1,2))/(picks(2,1)-picks(1,1)));
        end

        save(fullfile(folder,name,'angle.mat'),'angle','-v7.3');
    end
else
    disp('Loading angles.')
    load(fullfile(folder,name,'angle.mat'));
end

if apply_to_tslnum==1
 % apply filter to tslnum:
    disp('Apply filter to tslnum')
    tsl_filtered=data;
    for i=1:anz % for every part
        part=data.*mask{i};
        part(isnan(part))=0;
        
        sname=fullfile(folder,name,['Spectrum_Part_',int2str(i),'.png']);
        data_neu=kk_pieslicefilter(part,xgrid,ygrid,dx,angle(i),opening,plotflag,colorclip,sname);

        tsl_filtered(logical(mask{i}) & logical(gesamtmaske))=data_neu(logical(mask{i}) & logical(gesamtmaske));
    end

    % Plot figure and save
    fh=figure('Position',[0 0 1400 700]);
    subplot(1,2,1)
    imagesc(xgrid(1,:),ygrid(:,1),data)
    axis xy
    axis equal
    title(['Tsl ',int2str(tslnum),' original'])
    colormap(flipud(gray))
    coldata=sort(unique(data(:)));
    coldata(isnan(coldata))=[]; % delete nans
    cmin=coldata(round(length(coldata)/100*colorclip));
    cmax=coldata(end-round(length(coldata)/100*colorclip));
    set(gca,'CLim',[cmin cmax])

    subplot(1,2,2)
    imagesc(xgrid(1,:),ygrid(:,1),tsl_filtered)
    axis xy
    axis equal
    title(['Tsl ',int2str(tslnum),' filtered'])
    colormap(flipud(gray))
    coldata=sort(unique(tsl_filtered(:)));
    coldata(isnan(coldata))=[]; % delete nans
    cmin=coldata(round(length(coldata)/100*colorclip));
    cmax=coldata(end-round(length(coldata)/100*colorclip));
    set(gca,'CLim',[cmin cmax])

    saveas(fh,fullfile(folder,name,['Comparison_Tsl',int2str(tslnum),'.png']),'png');
    close(fh);
    disp(['Figures of spectrum and tslnum saved in folder ',name])
end

%% apply kk-filter to all timselices in different areas
if apply_to_all_tsl==1
    disp('Apply filter to all timeslices:')
    for j=1:length(tsl_interp)
        disp(['   ',int2str(j),'...'])
        data=tsl_interp{j};
        gesamtmaske=zeros(size(data));
        gesamtmaske(~isnan(data))=1;
        tsl_filtered=data;
        for i=1:anz % for every part
            part=data.*mask{i};
            part(isnan(part))=0;

            data_neu=kk_pieslicefilter(part,xgrid,ygrid,dx,angle(i),opening,0,colorclip,'');

            tsl_filtered(logical(mask{i}) & logical(gesamtmaske))=data_neu(logical(mask{i}) & logical(gesamtmaske));
        end
        % reset in tsl_interp:
        tsl_interp{j}=tsl_filtered;
    end

    % save readme
    disp('Saving ReadMe.')
    fid=fopen(fullfile(folder,name,'kk_filter_ReadMe.txt'),'wt');
    fprintf(fid,'kk-pieslice-filter with cosine taper\n');
    fprintf(fid,'Angles [°] for different parts:\n');
    fprintf(fid,'%.2f\n',angle);
    fprintf(fid,'Opening angle: %.1f°');
    if all_area==0
        fprintf(fid,'Parts are defined by polygons in shape-file:\n');
        fprintf(fid,'%s',fullfile(folder,shpfile));
    end
    fclose(fid);

    % save in MultichannelGPR-format (Use save_Tsl_Figures.m later for plotting)
    disp('Saving data in new folder...')
    save(fullfile(folder,name,'tsl_interp.mat'),'tsl_interp','-v7.3');
    copyfile(fullfile(folder,'coordtrans.mat'),fullfile(folder,name,'coordtrans.mat'));
    copyfile(fullfile(folder,'depth.mat'),fullfile(folder,name,'depth.mat'));
    copyfile(fullfile(folder,'mask_interp.mat'),fullfile(folder,name,'mask_interp.mat'));
    copyfile(fullfile(folder,'t_startende.mat'),fullfile(folder,name,'t_startende.mat'));
    copyfile(fullfile(folder,'topo_interp.mat'),fullfile(folder,name,'topo_interp.mat'));
    copyfile(fullfile(folder,'xgrid_interp.mat'),fullfile(folder,name,'xgrid_interp.mat'));
    copyfile(fullfile(folder,'ygrid_interp.mat'),fullfile(folder,name,'ygrid_interp.mat'));
    copyfile(fullfile(folder,'tslinfo.txt'),fullfile(folder,name,'tslinfo.txt'));

    disp('Done!')
end


% set original path
path(oldpath);

% End of script.
