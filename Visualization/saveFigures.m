clear all
close all
clc


% Read radargrams (radargrams.mat, global_coords.mat, x.mat and t.mat) from
% make_Radargrams or Bins2radargrams or other radargrams in same format
% (e.g. from DZT_Convert) and save them as figures (png and eps) and
% georeferenced PNGs.
% Normal figures can also be plotted as wiggle plot, georeferenced pictures
% in grayscale only!
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de


numbers=[]; % give numbers of processed radargrams or leave empty =[] for all

% Plotting options for radargrams:
colorclip=3; % 0 is colorscale from min(data) to max(data), 1 is 1% clip value, 2 is 2% clip value and 3 is 3% clip value, ... (will not be saved, for plotting only!)
aspectratio_t=7;    % for time plots: give aspectratio for y-axis. If you want to plot over whole screen, set =0. If you want to make the taxis larger, make this number smaller.
aspectratio_z=0.3;    % for depth plots: give aspectratio for y-axis. If you want to plot over whole screen, set =0. If you want to make the zaxis larger, make this number smaller.
wiggleplot=0;  % =1: make wiggle plot instead of grayscale picture (=0)
wigglescale=1; % scaling factor for wiggle plot
tz_flag=2; % y-axis is 1=time [ns] or 2=depth [m]
ampSpec_flag=0; % 1: amplitude spectrum, 0: radargram
dxtick=[]; % spacing between x-ticks in m (leave empty for automatic determination)
dytick=[]; % spacing between y-ticks in ns or m (leave empty for automatic determination)

% plot picks on radargrams? (created with LayerPicking.m)
plot_picks=0; % =1: I want to plot picks on the radargrams, =0: no picks
flag=0; % =1: only lines, =2: only points, =0: all
linewidth=1; % linewidth
markersize=5; % markersize of points

% save as georeferenced png?
save_georef=1; % yes=1, no=0

% save profile coordinates as shape file?
save_shape=1; % yes=1, no=0
shapename='GPR_Profiles.shp'; % give name for shape file

% Plotting options for map:
plot_map=0; % plot map? yes=1, no=0
map_xlim=[]; % give x-limits of map [start end] (leave empty for automatic determination)
map_ylim=[]; % give y-limits of map [start end] (leave empty for automatic determination)
xtick=[]; % spacing between x-ticks for map in m (leave empty for automatic determination)
ytick=[]; % spacing between y-ticks for map in m (leave empty for automatic determination)

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

if plot_picks==1  %load picks.txt
    if ~ispc; menu('Choose *.txt-file with picks.','OK'); end
    if ispc
    if exist('picks.temp') % read last opened folder from temp.temp
        fid=fopen('picks.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
        fileattrib('picks.temp','-h');
        fid=fopen('picks.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('picks.temp','+h');
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',pfad_rad); % path to radargram-folder

        fid=fopen('picks.temp','wt');
        fprintf(fid,'%s',fullfile(folder,file));
        fclose(fid);
        fileattrib('picks.temp','+h');
    end
else
    if exist('.picks.temp') % read last opened folder from temp.temp
        fid=fopen('.picks.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',pfad_rad); % path to radargram-folder
    end

    fid=fopen('.picks.temp','wt');
    fprintf(fid,'%s',fullfile(folder,file));
    fclose(fid);
end
end


% temporarily set path to required scripts
oldpath=path;
addpath('../Subfunctions/');


%% read pick file
if plot_picks==1
    fid=fopen(fullfile(pfad_rad,file),'r');
    anz=fscanf(fid,'%d',1);
    for i=1:anz
        temp{i}=textscan(fid,'%s',4); % ID - name (1:line/2:points)
    end
    tempcol=textscan(fid,'%f%f%f',anz);
    col=[tempcol{1} tempcol{2} tempcol{3}];
    % make layer strings
    for i=1:anz
        layerlist{i}=temp{i}{1}{3};
        layerID(i)=str2num(temp{i}{1}{1});
    end
    temp=textscan(fid,'%f%f%f%f%d%d%d%d','Headerlines',2);
    % picks=[x z E N ID profnum linenum flag];
    picks=[temp{1} temp{2} temp{3} temp{4} double(temp{5}) double(temp{6}) double(temp{7}) double(temp{8})];
    fclose(fid);
end


%% Read data
disp('Reading data...')
temp=load(fullfile(pfad_rad,'global_coords.mat'));
global_coords=temp.global_coords; % global coordinates of starting end ending point
temp=load(fullfile(pfad_rad,'x.mat'));
x=temp.x;   % profile coordinates
temp=load(fullfile(pfad_rad,'radargrams.mat'));
data=temp.radargrams;   % radargrams
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

        if plot_picks==1
            p=unique(picks(picks(:,6)==kk,5)); % LayerIDs in this profile
            pi=picks(picks(:,6)==kk,:); % all picks in this profile
        end
        
        f=figure('Visible','off');
        
        if ampSpec_flag==0
            % plot processed radargram
            if wiggleplot==0
                imagesc(x{kk},t,datatraces)
            else
                wigglesc(datatraces,t,x{kk},wigglescale);
            end
            if plot_picks==1
                hold on
                for j=1:length(p)   % for all IDs
                    linenum=unique(pi(pi(:,5)==p(j),7));  % all line numbers for this ID
                    for k=1:length(linenum)
                        if pi(pi(:,5)==p(j),8)==1 & (flag==0 || flag==1) % line
                            plot(pi(pi(:,5)==p(j) & pi(:,7)==linenum(k),1),pi(pi(:,5)==p(j) & pi(:,7)==linenum(k),2),'linewidth',linewidth,'Color',col(p(j),:))
                        elseif pi(pi(:,5)==p(j),8)==2 & (flag==0 || flag==2) % points
                            plot(pi(pi(:,5)==p(j) & pi(:,7)==linenum(k),1),pi(pi(:,5)==p(j) & pi(:,7)==linenum(k),2),'*','Markersize',markersize,'linewidth',linewidth,'Color',col(p(j),:))
                        end
                    end
                end
            end
            grid on
            xlabel('x [m]')
            ylabel('t [ns]')
            axis ij
            if tz_flag==2
                xlabel('x [m]')
                ylabel('z [m]')
                axis xy
            end
            colormap(flipud(gray));
            if tz_flag==1
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
            xlim=get(gca,'Xlim');
            ylim=get(gca,'Ylim');
            if ~isempty(dxtick)
                set(gca,'XTick',[0:dxtick:xlim(2)])
            end
            if ~isempty(dytick)
                set(gca,'YTick',[round(ylim(1)):dytick:ylim(2)])
            end
            set(gca,'XLim',[0 x{kk}(end)],'YLim',[min([t(1) t(end)]) max([t(1) t(end)])])
        else
            % plot amplitude spectrum
            plot(t,datatraces)
            hold on
            plot(t,mean(datatraces,2),'k','Linewidth',2)
            xlabel('f [MHz]')
            ylabel('Amplitude')
            grid on
            xlim=get(gca,'Xlim');
            ylim=get(gca,'Ylim');
            if ~isempty(dxtick)
                set(gca,'XTick',[0:dxtick:xlim(2)])
            end
            if ~isempty(dytick)
                set(gca,'YTick',[0:dytick:ylim(2)])
            end
        end
        % Save figures
        if ~exist(fullfile(pfad_rad,'Figures'),'dir')
            mkdir(fullfile(pfad_rad,'Figures'));
        end
        saveas(f,fullfile(pfad_rad,'Figures',['Radargram_',int2str(kk),'.eps']),'epsc')
        saveas(f,fullfile(pfad_rad,'Figures',['Radargram_',int2str(kk),'.png']),'png')
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


%% plot map
if plot_map==1
    disp('Save map...')
    fm=figure('Visible','off');
    hold on
    led=[];
    for kk=1:length(numbers)
        quiver(global_coords{numbers(kk)}(1,1),global_coords{numbers(kk)}(1,2),global_coords{numbers(kk)}(end,1)-global_coords{numbers(kk)}(1,1),global_coords{numbers(kk)}(end,2)-global_coords{numbers(kk)}(1,2),0,'Linewidth',2);
        led=[led; {['Profile ',int2str(numbers(kk))]}];
        text(global_coords{numbers(kk)}(end,1)+1,global_coords{numbers(kk)}(end,2),int2str(numbers(kk)))
    end
    grid on
    legend(led,'Location','NorthEastOutside')
    if ~isempty(map_xlim)
        set(gca,'XLim',map_xlim)
    else
        map_xlim=get(gca,'XLim');
    end
    if ~isempty(map_ylim)
        set(gca,'yLim',map_ylim)
    else
        map_ylim=get(gca,'YLim');
    end
    if ~isempty(xtick)
        set(gca,'XTick',[map_xlim(1):xtick:map_xlim(2)])
    end
    if ~isempty(ytick)
        set(gca,'YTick',[map_ylim(1):ytick:map_ylim(2)])
    end
    xticks=get(gca,'XTick');
    yticks=get(gca,'YTick');
    for i=1:length(xticks)
        xticklabels{i}=sprintf('%d',xticks(i));
    end
    for i=1:length(yticks)
        yticklabels{i}=sprintf('%d',yticks(i));
    end
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0,'YTick',yticks,'YTickLabel',yticklabels,'DataAspectRatio',[1 1 1])
    
    saveas(fm,fullfile(pfad_rad,'Figures','Map.eps'),'epsc')
    saveas(fm,fullfile(pfad_rad,'Figures','Map.png'),'png')
end

%% save as shape file
if save_shape==1
    disp('Saving shape file...')
    
    for kk=1:length(numbers) % loop over radargrams
        S(kk).Geometry='Line';
        S(kk).BoundingBox=[min(global_coords{numbers(kk)}(:,1)) min(global_coords{numbers(kk)}(:,2)); max(global_coords{numbers(kk)}(:,1)) max(global_coords{numbers(kk)}(:,2))];
        S(kk).X=global_coords{numbers(kk)}(:,1);
        S(kk).Y=global_coords{numbers(kk)}(:,2);
        S(kk).id=numbers(kk);
    end
    
    % write shapefile with profile lines
    shapewrite(S,fullfile(pfad_rad,'Figures',shapename));
end


%% save as georef
if save_georef==1
    disp('Saving georeferenced radargrams...')
    if ampSpec_flag==0
        for kk=1:length(numbers) % loop over radargrams
            if ~isempty(data{numbers(kk)}) && any(~isnan(data{numbers(kk)}(:)))
                
                % Save figures
                if ~exist(fullfile(pfad_rad,'Figures','georef'),'dir')
                    mkdir(fullfile(pfad_rad,'Figures','georef'));
                end
                
                cdata=data{numbers(kk)}; % read radargrams
                
                % turn around that radargram is running from west to east
                if global_coords{numbers(kk)}(1,1)>global_coords{numbers(kk)}(end,1)
                    global_coords{numbers(kk)}=flipud(global_coords{numbers(kk)});
                    cdata=fliplr(cdata);
                end
                
                % prepare data and save
                cmin=min(cdata(:));
                cmax=max(cdata(:));
                coldata=sort(unique(cdata(~isnan(cdata(:)))));
                cmin=coldata(round(length(coldata)/100*colorclip));
                cmax=coldata(end-round(length(coldata)/100*colorclip));
                range=cmax-cmin;
                cdata=(cdata-cmin)/range;
                cdata(cdata<=0)=0;
                cdata(cdata>=1)=1;
                m=ones(size(cdata));
                m(isnan(cdata))=0;
                im=cdata.*256;
                im(im<=2)=2;
                im(isnan(cdata))=0;  % set nan to 0 -> will be transparent
                imwrite(im,flipud(gray(256)),fullfile(pfad_rad,'Figures','georef',['Radargram_',int2str(numbers(kk)),'.png']),'Transparency',0);
                
                % write pngw
                coordtrans=[0 0 global_coords{numbers(kk)}(1,:); x{numbers(kk)}(end) 0 global_coords{numbers(kk)}(end,:)];
                if length(x{numbers(kk)}(:,1))>length(x{numbers(kk)}(1,:))
                    x{numbers(kk)}=x{numbers(kk)}'; % make row vector
                end
                if tz_flag==1 % time
                    write_geoPNGW_radargrams(repmat(x{numbers(kk)},[length(t) 1]),repmat(t'./aspectratio_t,[1 length(x{numbers(kk)})]),coordtrans,fullfile(pfad_rad,'Figures','georef',['Radargram_',int2str(numbers(kk)),'.pgw']));
                else % depth
                    if aspectratio_z==0
                        aspectratio_z=1; % set to 1
                    end
                    ygrid=-1.*(repmat(t'.*aspectratio_z,[1 length(x{numbers(kk)})])-t(1)); % transform that z=0 is at surface and then positive down
                    write_geoPNGW_radargrams(repmat(x{numbers(kk)},[length(t) 1]),ygrid,coordtrans,fullfile(pfad_rad,'Figures','georef',['Radargram_',int2str(numbers(kk)),'.pgw']));
                end
                
            end
            disp(['  ',int2str(numbers(kk))])
        end
    else
        disp('Data are amplitude spectra, no referenced pngs are exported.')
    end
    
end

disp('Done!')

% restore original path
path(oldpath);