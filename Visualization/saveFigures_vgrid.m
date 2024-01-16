clear all
close all
clc


% Read vgrid.mat (and global_coords.mat, x.mat and t.mat) and save them as
% figures (png and eps) and georeferenced PNGs.
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de


numbers=[1]; % give numbers of processed radargrams or leave empty =[] for all

% Plotting options for radargrams:
colorauto=0; % 0 is manual colorscale selection, 1 is automatic color range (min-max)
clim=[0.06 0.12]; % colorlimits in m/ns that will be used when colorauto=0
aspectratio_t=7;    % for time plots: give aspectratio for y-axis. If you want to plot over whole screen, set =0. If you want to make the taxis larger, make this number smaller.

dxtick=[]; % spacing between x-ticks in m (leave empty for automatic determination)
dytick=[]; % spacing between y-ticks in ns or m (leave empty for automatic determination)

% save as georeferenced png?
save_georef=1; % yes=1, no=0

% -------------------------------------------------------------------------
% Do not change the following part!

% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with vgrid.mat','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with vgrid.mat');
        else
            pfad_rad=uigetdir([],'Choose folder with vgrid.mat');
        end
        fileattrib('radtemp.temp','-h');
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    else
        pfad_rad=uigetdir([],'Choose folder with vgrid.mat'); % path to radargram-folder

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
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with vgrid.mat');
        else
            pfad_rad=uigetdir([],'Choose folder with vgrid.mat');
        end
    else
        pfad_rad=uigetdir([],'Choose folder with vgrid.mat'); % path to radargram-folder
    end

    fid=fopen('.radtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end


% temporarily set path to required scripts
oldpath=path;
addpath('../Subfunctions/');


%%% Read data
disp('Reading data...')
temp=load(fullfile(pfad_rad,'global_coords.mat'));
global_coords=temp.global_coords; % global coordinates of starting end ending point
temp=load(fullfile(pfad_rad,'x.mat'));
x=temp.x;   % profile coordinates
temp=load(fullfile(pfad_rad,'vgrid.mat'));
data=temp.vgrid;   % vgrid
temp=load(fullfile(pfad_rad,'t.mat'));
t=temp.t;   % time vector
dt=t(2)-t(1);

if isempty(numbers)
    numbers=1:length(data); % all radargrams
end

%-------------------------------------------------------------------------
%% plot as invisible figure and save as png and eps
disp('Saving figures...')
weg=[];
for kk=numbers % loop over radargrams
    if all(global_coords{kk}(:)==0) % if all coordinates ==0 -> do not use profile for map
        weg=[weg; kk];
    end
    if ~isempty(data{kk}) && any(~isnan(data{kk}(:))) && ~all(global_coords{kk}(:)==0)

        datatraces=data{kk}; % read radargrams

        f=figure('Visible','off');

        % plot v
        imagesc(x{kk},t,datatraces)
        grid on
        xlabel('x [m]')
        ylabel('t [ns]')
        axis ij
        colormap(jet);
        set(gca,'Dataaspectratio',[1 aspectratio_t 1])

        if colorauto==0
            set(gca,'CLim',clim)
        else
            set(gca,'ClimMode','auto');
        end
        xlim=get(gca,'Xlim');
        ylim=get(gca,'Ylim');
        if ~isempty(dxtick)
            set(gca,'XTick',[0:dxtick:xlim(2)])
        end
        if ~isempty(dytick)
            set(gca,'YTick',[round(ylim(1)):dytick:ylim(2)])
        end
        set(gca,'XLim',[0 x{kk}(end)],'YLim',[min([t(1) t(end)]) max([t(1) t(end)])])
        cb=colorbar;
        ylabel(cb,'Velocity [m/ns]');

        % Save figures
        if ~exist(fullfile(pfad_rad,'Figures_vmodel'),'dir')
            mkdir(fullfile(pfad_rad,'Figures_vmodel'));
        end
        saveas(f,fullfile(pfad_rad,'Figures_vmodel',['vgrid_',int2str(kk),'.eps']),'epsc')
        saveas(f,fullfile(pfad_rad,'Figures_vmodel',['vgrid_',int2str(kk),'.png']),'png')
    end
    disp(['  ',int2str(kk)])
end

% do not use profiles without coordinates
numbers(weg)=[];
if ~isempty(weg)
    disp('Following files do not have coordinates:')
    for i=1:length(weg)
        disp(int2str(weg))
    end
end


%% save as georef
if save_georef==1
    disp('Saving georeferenced v-models...')

    for kk=1:length(numbers) % loop over radargrams
        if ~isempty(data{numbers(kk)}) && any(~isnan(data{numbers(kk)}(:)))

            % Save figures
            if ~exist(fullfile(pfad_rad,'Figures_vmodel','georef'),'dir')
                mkdir(fullfile(pfad_rad,'Figures_vmodel','georef'));
            end

            cdata=data{numbers(kk)}; % read radargrams

            % turn around that radargram is running from west to east
            if global_coords{numbers(kk)}(1,1)>global_coords{numbers(kk)}(end,1)
                global_coords{numbers(kk)}=flipud(global_coords{numbers(kk)});
                cdata=fliplr(cdata);
            end

            % prepare data and save
            if colorauto==1
                cmin=min(cdata(:));
                cmax=max(cdata(:));
            else
                cmin=clim(1);
                cmax=clim(2);
            end
            range=cmax-cmin;
            cdata=(cdata-cmin)/range;
            cdata(cdata<=0)=0;
            cdata(cdata>=1)=1;
            m=ones(size(cdata));
            m(isnan(cdata))=0;
            im=cdata.*256;
            im(im<=2)=2;
            im(isnan(cdata))=0;  % set nan to 0 -> will be transparent
            imwrite(im,jet,fullfile(pfad_rad,'Figures_vmodel','georef',['vgrid_',int2str(numbers(kk)),'.png']),'Transparency',0);

            % write pngw
            coordtrans=[0 0 global_coords{numbers(kk)}(1,:); x{numbers(kk)}(end) 0 global_coords{numbers(kk)}(end,:)];
            if length(x{numbers(kk)}(:,1))>length(x{numbers(kk)}(1,:))
                x{numbers(kk)}=x{numbers(kk)}'; % make row vector
            end

            write_geoPNGW_radargrams(repmat(x{numbers(kk)},[length(t) 1]),repmat(t'./aspectratio_t,[1 length(x{numbers(kk)})]),coordtrans,fullfile(pfad_rad,'Figures_vmodel','georef',['vgrid_',int2str(numbers(kk)),'.pgw']));

        end
        disp(['  ',int2str(numbers(kk))])
    end

end

disp('Done!')

% restore original path
path(oldpath);