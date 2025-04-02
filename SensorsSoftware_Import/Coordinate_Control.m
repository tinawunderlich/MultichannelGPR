clear all
close all
clc

% Script for reading radargrams in radargrams.mat format and manipulating
% the start/end coordinates of the profiles
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files (please choose folder with these files):
% radargrams.mat: Radar data
% global_coords.mat: coordinates
% t.mat: time vector
%
% requires MATLAB-files in following folders (path will be temporarily
% set): Export_Import, Subfunctions


% Option 1: Read coordinates from mat-files and export to txt.
% Option 2: Read coordinates from txt-file and set into global_coords.mat.
% Option 3: Manipulate coordinates with settings below (no text-file written/required).
% Option 4: For skewed areas, when zigzag profile shift is linearly changing over the area
option=2;

% settings for option 3:
x_offset=50; % offset in x-direction [m] (e.g. x_offset=1: all profiles will be shifted 1 m towards positive x-axis)
                % if e.g. x_offset=[-0.3 0.3]: all odd profile numbers will
                % be shifted 0.3m towards negative x-axis and even profile
                % numbers 0.3m towards positive x-axis
y_offset=0; % offset in y-direction [m] (same explanation as for x)

% settings for option 4:
x_offset_start=0;   % offset for small y [m], skalar or 2 numbers (see above)
x_offset_end=0;     % offset for maximum y [m], skalar or 2 numbers (see above)
y_offset_start=[0 0];   % offset for small x [m], skalar or 2 numbers (see above)
y_offset_end=[1.2 -1.2];     % offset for maximum x [m], skalar or 2 numbers (see above) 


txtname='coords.txt'; % give name of coordinate text-file (saved in same folder as radargrams)

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
            foldername=uigetdir(fn{1}{1},'Choose folder with data');
        else
            foldername=uigetdir([],'Choose folder with data');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
    else
        foldername=uigetdir([],'Choose folder with data'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Choose folder with data');
        else
            foldername=uigetdir([],'Choose folder with data');
        end
    else
        foldername=uigetdir([],'Choose folder with data'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));



%% Load all profiles
load(fullfile(foldername,'radargrams.mat'));
load(fullfile(foldername,'global_coords.mat'));
load(fullfile(foldername,'t.mat'));
load(fullfile(foldername,'x.mat'));

dt=t(2)-t(1);
ns=length(t);
dx=mean(diff(x{1})); % just approx.

if option==1
    % plot coordinate map
    figure
    hold on
    for i=1:length(x)
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'b','Linewidth',2)
        textcoord(1)=interp1(1:length(x{i}),global_coords{i}(:,1),-1/dx,'linear','extrap')-0.2;
        textcoord(2)=interp1(1:length(x{i}),global_coords{i}(:,2),-1/dx,'linear','extrap');
        text(textcoord(1),textcoord(2),int2str(i),'Fontsize',12);
        minx(i,1)=min(global_coords{i}(:,1));
        maxx(i,1)=max(global_coords{i}(:,1));
        miny(i,1)=min(global_coords{i}(:,2));
        maxy(i,1)=max(global_coords{i}(:,2));
        startx(i,1)=global_coords{i}(1,1);
        starty(i,1)=global_coords{i}(1,2);
        endx(i,1)=global_coords{i}(end,1);
        endy(i,1)=global_coords{i}(end,2);
    end
    xlabel('x [m]')
    ylabel('y [m]')
    title('Profile plan (numbers at beginning of profile)')
    grid on
    set(gca,'Fontsize',15,'XLim',[min(minx)-2 max(maxx)+2],'YLim',[min(miny)-2 max(maxy)+2])
    saveas(gcf,fullfile(foldername,'ProfilePlan.png'),'png'); % save figure

    % write txt-file
    fid=fopen(fullfile(foldername,txtname),'wt');
    fprintf(fid,'ProfNr.\tx_start[m]\ty_start[m]\tx_end[m]\ty_end[m]\n');
    fprintf(fid,'%d\t%.2f\t%.2f\t%.2f\t%.2f\n',[[1:length(x)]' startx starty endx endy]');
    fclose(fid);

    disp(['File ',txtname,' saved. Please edit manually and restart script with option=2.'])

elseif option==2
    % Read txt-file and set new coords:
    if exist(fullfile(foldername,txtname),'file')
        fid=fopen(fullfile(foldername,txtname),'r');
        temp=textscan(fid,'%f%f%f%f%f','Headerlines',1);
        fclose(fid);

        coords=[temp{1} temp{2} temp{3} temp{4} temp{5}];
    else
        disp(['Error: No file ',txtname,' found.']);
        return;
    end

    % set coords:
    for i=1:length(x)
        n=length(x{i}); % number of traces
        global_coords{i}(:,1)=linspace(coords(i,2),coords(i,4),n);
        global_coords{i}(:,2)=linspace(coords(i,3),coords(i,5),n);
    end

    % create new folder:
    if ~exist(fullfile(foldername,'newCoords'),'dir')
        mkdir(fullfile(foldername,'newCoords'));
    end

    % plot coordinate map
    figure
    hold on
    for i=1:length(x)
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'b','Linewidth',2)
        textcoord(1)=interp1(1:length(x{i}),global_coords{i}(:,1),-1/dx,'linear','extrap')-0.2;
        textcoord(2)=interp1(1:length(x{i}),global_coords{i}(:,2),-1/dx,'linear','extrap');
        text(textcoord(1),textcoord(2),int2str(i),'Fontsize',12);
        minx(i,1)=min(global_coords{i}(:,1));
        maxx(i,1)=max(global_coords{i}(:,1));
        miny(i,1)=min(global_coords{i}(:,2));
        maxy(i,1)=max(global_coords{i}(:,2));
    end
    xlabel('x [m]')
    ylabel('y [m]')
    title('Profile plan (numbers at beginning of profile)')
    grid on
    set(gca,'Fontsize',15,'XLim',[min(minx)-2 max(maxx)+2],'YLim',[min(miny)-2 max(maxy)+2])
    saveas(gcf,fullfile(foldername,'newCoords','ProfilePlan.png'),'png'); % save figure

    % save data
    copyfile(fullfile(foldername,'radargrams.mat'),fullfile(foldername,'newCoords','radargrams.mat'));
    copyfile(fullfile(foldername,'t.mat'),fullfile(foldername,'newCoords','t.mat'));
    copyfile(fullfile(foldername,'x.mat'),fullfile(foldername,'newCoords','x.mat'));
    save(fullfile(foldername,'newCoords','global_coords.mat'),'global_coords','-v7.3');

    disp('Saved data with new coordinates in folder newCoords.')

elseif option==3  % manipulate coords
    % create new folder:
    if ~exist(fullfile(foldername,'newCoords'),'dir')
        mkdir(fullfile(foldername,'newCoords'));
    end

    % plot coordinate map
    figure
    hold on
    for i=1:length(x)
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'b','Linewidth',2)
        textcoord(1)=interp1(1:length(x{i}),global_coords{i}(:,1),-1/dx,'linear','extrap')-0.2;
        textcoord(2)=interp1(1:length(x{i}),global_coords{i}(:,2),-1/dx,'linear','extrap');
        text(textcoord(1),textcoord(2),int2str(i),'Fontsize',12);
        minx(i,1)=min(global_coords{i}(:,1));
        maxx(i,1)=max(global_coords{i}(:,1));
        miny(i,1)=min(global_coords{i}(:,2));
        maxy(i,1)=max(global_coords{i}(:,2));
    end
    xlabel('x [m]')
    ylabel('y [m]')
    title('Profile plan (numbers at beginning of profile) - blue=old/red=new')
    grid on
    
    % set coords:
    for i=1:length(x)
        if length(x_offset)==1
            xo=x_offset;
        elseif length(x_offset)==2
            xo=x_offset(2-mod(i,2));
        else
            disp('Wrong format of x_offset. Please give one or two numbers.')
            return;
        end
        if length(y_offset)==1
            yo=y_offset;
        elseif length(y_offset)==2
            yo=y_offset(2-mod(i,2));
        else
            disp('Wrong format of y_offset. Please give one or two numbers.')
            return;
        end
        global_coords{i}(:,1)=global_coords{i}(:,1)+xo;
        global_coords{i}(:,2)=global_coords{i}(:,2)+yo;
        minx(i+length(x),1)=min(global_coords{i}(:,1));
        maxx(i+length(x),1)=max(global_coords{i}(:,1));
        miny(i+length(x),1)=min(global_coords{i}(:,2));
        maxy(i+length(x),1)=max(global_coords{i}(:,2));
        % plot line
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'r','Linewidth',2)
    end
    set(gca,'Fontsize',15,'XLim',[min(minx)-2 max(maxx)+2],'YLim',[min(miny)-2 max(maxy)+2])
    saveas(gcf,fullfile(foldername,'newCoords','ProfilePlan.png'),'png'); % save figure

    % save data
    copyfile(fullfile(foldername,'radargrams.mat'),fullfile(foldername,'newCoords','radargrams.mat'));
    copyfile(fullfile(foldername,'t.mat'),fullfile(foldername,'newCoords','t.mat'));
    copyfile(fullfile(foldername,'x.mat'),fullfile(foldername,'newCoords','x.mat'));
    save(fullfile(foldername,'newCoords','global_coords.mat'),'global_coords','-v7.3');

    fid=fopen(fullfile(foldername,'newCoords','log.txt'),'wt');
    fprintf(fid,'Coordinate shift:\n');
    if length(x_offset)==1
        fprintf(fid,'x-offset: %.2f m\n',x_offset);
    elseif length(x_offset)==2
        fprintf(fid,'x-offset (odd profile numbers): %.2f m\nx-offset (even profile numbers): %.2f m\n',x_offset);
    end
    if length(y_offset)==1
        fprintf(fid,'y-offset: %.2f m\n',y_offset);
    elseif length(y_offset)==2
        fprintf(fid,'y-offset (odd profile numbers): %.2f m\ny-offset (even profile numbers): %.2f m\n',y_offset);
    end
    fclose(fid);

    disp('Saved data with new coordinates in folder newCoords.')

elseif option==4  % area skewed, linearly changing shifts
    % create new folder:
    if ~exist(fullfile(foldername,'newCoords'),'dir')
        mkdir(fullfile(foldername,'newCoords'));
    end

    % plot coordinate map
    figure
    hold on
    for i=1:length(x)
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'b','Linewidth',2)
        textcoord(1)=interp1(1:length(x{i}),global_coords{i}(:,1),-1/dx,'linear','extrap')-0.2;
        textcoord(2)=interp1(1:length(x{i}),global_coords{i}(:,2),-1/dx,'linear','extrap');
        text(textcoord(1),textcoord(2),int2str(i),'Fontsize',12);
        minx(i,1)=min(global_coords{i}(:,1));
        maxx(i,1)=max(global_coords{i}(:,1));
        miny(i,1)=min(global_coords{i}(:,2));
        maxy(i,1)=max(global_coords{i}(:,2));
    end
    xlabel('x [m]')
    ylabel('y [m]')
    title('Profile plan (numbers at beginning of profile) - blue=old/red=new')
    grid on
    
    % set coords:
    for i=1:length(x)
        % x start offset:
        if length(x_offset_start)==1
            xstart=x_offset_start;
        elseif length(x_offset_start)==2
            xstart=x_offset_start(2-mod(i,2));
        else
            disp('Wrong format of x_offset_start. Please give one or two numbers.')
            return;
        end
        % x end offset:
        if length(x_offset_end)==1
            xend=x_offset_end;
        elseif length(x_offset_end)==2
            xend=x_offset_end(2-mod(i,2));
        else
            disp('Wrong format of x_offset_end. Please give one or two numbers.')
            return;
        end
        % y start offset:
        if length(y_offset_start)==1
            ystart=y_offset_start;
        elseif length(y_offset_start)==2
            ystart=y_offset_start(2-mod(i,2));
        else
            disp('Wrong format of y_offset_start. Please give one or two numbers.')
            return;
        end
        % y end offset:
        if length(y_offset_end)==1
            yend=y_offset_end;
        elseif length(y_offset_end)==2
            yend=y_offset_end(2-mod(i,2));
        else
            disp('Wrong format of y_offset_end. Please give one or two numbers.')
            return;
        end
        % interpolate offset for this profile
        if xstart==xend
            xo=xstart;
        else
            xo=interp1([1 length(x)],[xstart xend],i); % linear interpolation of x-offset for this profile
        end
        if ystart==yend
            yo=ystart;
        else
            yo=interp1([1 length(x)],[ystart yend],i); % linear interpolation of x-offset for this profile
        end
        global_coords{i}(:,1)=global_coords{i}(:,1)+xo;
        global_coords{i}(:,2)=global_coords{i}(:,2)+yo;
        minx(i+length(x),1)=min(global_coords{i}(:,1));
        maxx(i+length(x),1)=max(global_coords{i}(:,1));
        miny(i+length(x),1)=min(global_coords{i}(:,2));
        maxy(i+length(x),1)=max(global_coords{i}(:,2));
        % plot line
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'r','Linewidth',2)
    end
    set(gca,'Fontsize',15,'XLim',[min(minx)-2 max(maxx)+2],'YLim',[min(miny)-2 max(maxy)+2])
    saveas(gcf,fullfile(foldername,'newCoords','ProfilePlan.png'),'png'); % save figure

    % save data
    copyfile(fullfile(foldername,'radargrams.mat'),fullfile(foldername,'newCoords','radargrams.mat'));
    copyfile(fullfile(foldername,'t.mat'),fullfile(foldername,'newCoords','t.mat'));
    copyfile(fullfile(foldername,'x.mat'),fullfile(foldername,'newCoords','x.mat'));
    save(fullfile(foldername,'newCoords','global_coords.mat'),'global_coords','-v7.3');

    fid=fopen(fullfile(foldername,'newCoords','log.txt'),'wt');
    fprintf(fid,'Coordinate shift:\n');
    if length(x_offset_start)==1
        fprintf(fid,'x-offset-start: %.2f m\n',x_offset_start);
    elseif length(x_offset_start)==2
        fprintf(fid,'x-offset-start (odd profile numbers): %.2f m\nx-offset-start (even profile numbers): %.2f m\n',x_offset_start);
    end
    if length(x_offset_end)==1
        fprintf(fid,'x-offset-end: %.2f m\n',x_offset_end);
    elseif length(x_offset_end)==2
        fprintf(fid,'x-offset-end (odd profile numbers): %.2f m\nx-offset-end (even profile numbers): %.2f m\n',x_offset_end);
    end
    if length(y_offset_start)==1
        fprintf(fid,'y-offset-start: %.2f m\n',y_offset_start);
    elseif length(y_offset_start)==2
        fprintf(fid,'y-offset-start (odd profile numbers): %.2f m\ny-offset-start (even profile numbers): %.2f m\n',y_offset_start);
    end
    if length(y_offset_end)==1
        fprintf(fid,'y-offset-end: %.2f m\n',y_offset_end);
    elseif length(y_offset_end)==2
        fprintf(fid,'y-offset-end (odd profile numbers): %.2f m\ny-offset-end (even profile numbers): %.2f m\n',y_offset_end);
    end
    fclose(fid);

    disp('Saved data with new coordinates in folder newCoords.')
end
