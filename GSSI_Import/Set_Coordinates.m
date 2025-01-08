clear all
close all
clc

%------------------------ Set coordinates to markers -------------------------------------
%
% Reads mat-files from DZT_Convert.m and sets coordinates to marker
% positions
%
% Dr. Tina Wunderlich, August 2020, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


% Mode of usage:
mode=3; % 1: plot profiles with markers and get marker numbers, export table
        % 2: set/delete markers (you can skip this step if not necessary)
        % 3: read table and set coordinates

mode2_options=[1 -1]; % for mode==2: list of profilenum and marker (+ means add marker, - means delete marker)

dx=0.02;    % make constant trace spacing in m (if not, then =0)

export2sgy=0;
export2mat=1;

%--------------------------------------------------------------------------
%%% DO NOT CHANGE FORM HERE ON!

warning off

test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));

if ~ispc; menu('Choose folder with mat-file(s)','OK'); end
if ispc
    if exist('radtemp') % read last opened folder from temp.temp
        fid=fopen('rad.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with mat-file(s)');
        else
            pfad=uigetdir([],'Choose folder with mat-file(s)');
        end
        fid=fopen('rad.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with mat-file(s)'); % path to radargram-folder
        fid=fopen('rad.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.rad.temp') % read last opened folder from temp.temp
        fid=fopen('.rad.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with mat-file(s)');
        else
            pfad=uigetdir([],'Choose folder with mat-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with mat-file(s)'); % path to radargram-folder
    end

    fid=fopen('.rad.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

%% read data
load(fullfile(pfad,'radargrams.mat'));
load(fullfile(pfad,'t.mat'));
load(fullfile(pfad,'global_coords.mat'));
load(fullfile(pfad,'x.mat'));
load(fullfile(pfad,'marker.mat'));
fid=fopen(fullfile(pfad,'radargrams.txt'),'r');
temp=textscan(fid,'%f%s%f','Headerlines',1);
fclose(fid);
list.num=temp{1};
list.name=temp{2};
list.chan=temp{3};
nametemp=split(list.name{1},'_');
name=nametemp{1}; % name of data files

%% different modes
if mode==1
    for i=1:length(x)
        % get marker
        marks{i}=find(marker{i}==1);
        
        % plot figure
        figure
        imagesc(1:length(radargrams{i}(1,:)),t,radargrams{i})
        hold on
        for j=1:length(marks{i})
            plot([marks{i}(j) marks{i}(j)],[0 max(t)/10],'k')
        end
        colorbar
        colormap(flipud(gray))
        title([int2str(list.num(i)),': ',list.name{i},' Chan ',int2str(list.chan(i))])
    end
    
    % export infos
    exportInfo(pfad,name,radargrams,global_coords,marks,list);
    
elseif mode==2
    % go through profiles and add/delete markers
    for i=1:length(mode2_options(:,1))
        prof=mode2_options(i,1);    % profile number
        add_del=mode2_options(i,2); % number of markers to add or delete
        
        if add_del<0
            disp(['Profile ',int2str(prof),': Delete ',int2str(abs(add_del)),' marker(s) by clicking.'])
        else
            disp(['Profile ',int2str(prof),': Add ',int2str(add_del),' marker(s) by clicking.'])
        end
            
        % get marker
        marks{prof}=find(marker{prof}==1);
        % plot figure
        figure
        imagesc(1:length(radargrams{prof}(1,:)),t,radargrams{prof})
        hold on
        for j=1:length(marks{prof})
            plot([marks{prof}(j) marks{prof}(j)],[0 max(t)/10],'k')
        end
        colorbar
        colormap(flipud(gray))
        title([int2str(list.num(prof)),': ',list.name{prof},' Chan ',int2str(list.chan(prof))])
        [xpick,ypick]=ginput(abs(add_del)); % picks are trace number and time
        xpick=round(xpick);
        
        if add_del<0
            % delete
            for j=1:length(xpick) % for every pick find nearest marker
                dist=abs(xpick(j)-marks{prof});
                marks{prof}(dist==min(dist))=[];
            end   
        else
            % add
            marks{prof}=sort([marks{prof} xpick]);
        end
        % update marker
        marker{prof}=zeros(size(marker{prof}));
        marker{prof}(marks{prof})=1;
        
        %plot new
        hold off
        imagesc(1:length(radargrams{prof}(1,:)),t,radargrams{prof})
        hold on
        for j=1:length(marks{prof})
            plot([marks{prof}(j) marks{prof}(j)],[0 max(t)/10],'k')
        end
        colorbar
        colormap(flipud(gray))
        title([int2str(list.num(prof)),': ',list.name{prof},' Chan ',int2str(list.chan(prof))])
    end
    
    % get marker for non-corrected profiles
    for i=1:length(list.num)
        if ~any(list.num(i)==mode2_options(:,1))
            % get marker
            marks{i}=find(marker{i}==1);
        end
    end
    
    if matlab==1
        save(fullfile(pfad,'marker.mat'),'marker','-v7.3');
    else
        save(fullfile(pfad,'marker.mat'),'marker');
    end
    
    disp('All markers corrected. Exporting new coordinate file.')
    
    % export infos
    exportInfo(pfad,name,radargrams,global_coords,marks,list);
    
    disp('Please add coordinates in this file and change to mode==3.')
    
elseif mode==3
    % read coordinate file and set coordinates
    fid=fopen(fullfile(pfad,[name,'_coordinates.txt']),'r');
    temp=textscan(fid,'%f%f%f%f%f%f','Headerlines',1);
    fclose(fid);
    coords=[temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}]; % profnum, xstart, ystart, xende, yende, anzMarker
    
    for i=1:length(x) % for all profiles
        % get marker
        marks{i}=find(marker{i}==1);
        
        % cut radargrams before first and after last marker
        radnew=radargrams{i}(:,marks{i}(1):marks{i}(end));
        
        % coordinates at markers
        xm=interp1([1 coords(i,6)],[coords(i,2) coords(i,4)],1:coords(i,6));
        ym=interp1([1 coords(i,6)],[coords(i,3) coords(i,5)],1:coords(i,6));
        
        % interpolate between markers
        xc{i}=interp1(marks{i},xm,marks{i}(1):marks{i}(end));
        yc{i}=interp1(marks{i},ym,marks{i}(1):marks{i}(end));

        % make *cm trace spacing
        if dx~=0
            pr=[0 cumsum(sqrt(diff(xc{i}).^2+diff(yc{i}).^2))]; % inline profile coordinate
            prnew{i}=[0:dx:pr(end)]; % new profile coordinate with dx trace spacing
            rad{i}=zeros(length(t),length(prnew{i}));
            for j=1:length(t) % for all samples...
                rad{i}(j,:)=interp1(pr,radnew(j,:),prnew{i});
            end
            xc{i}=interp1(pr,xc{i},prnew{i});
            yc{i}=interp1(pr,yc{i},prnew{i});
        else
            rad{i}=radnew;
            prnew{i}=[0 cumsum(sqrt(diff(xc{i}).^2+diff(yc{i}).^2))];
        end
        
        global_coords{i}=[xc{i}' yc{i}'];
    end
    radargrams=rad';
    x=prnew';
    
    
    % save data with coordinates
    if ~exist(fullfile(pfad,'setCoords'),'dir')
        mkdir(fullfile(pfad,'setCoords'));
    end
    if matlab==1
        save(fullfile(pfad,'setCoords','radargrams.mat'),'radargrams','-v7.3');
        save(fullfile(pfad,'setCoords','t.mat'),'t','-v7.3');
        save(fullfile(pfad,'setCoords','x.mat'),'x','-v7.3');
        save(fullfile(pfad,'setCoords','global_coords.mat'),'global_coords','-v7.3');
    else
        save(fullfile(pfad,'setCoords','radargrams.mat'),'radargrams');
        save(fullfile(pfad,'setCoords','t.mat'),'t');
        save(fullfile(pfad,'setCoords','x.mat'),'x');
        save(fullfile(pfad,'setCoords','global_coords.mat'),'global_coords');
    end
    
    disp('Data with coordinates saved in folder setCoords.')
    
    if export2sgy==1
        %%% exporting as sgy
        disp('Start saving as sgy...')
        % make new folder
        if ~exist(fullfile(pfad,'setCoords','sgy'),'dir')
            mkdir(fullfile(pfad,'setCoords','sgy'));
        end
        
        dt=t(2)-t(1);
        for ii=1:length(radargrams)
            disp([int2str(ii),' / ',int2str(length(radargrams))])
            export2sgy2D(radargrams{ii},dt,global_coords{ii}(:,1),global_coords{ii}(:,2),fullfile(pfad,'setCoords','sgy',[list.name{ii},'_Chan',int2str(list.chan(ii)),'.sgy']));
        end
        disp('   Finished!')
    end
end

% set original path
path(oldpath);


function exportInfo(pfad,name,radargrams,global_coords,marks,list)
for i=1:length(marks)
    anzmarks(i,1)=length(marks{i});
    
    startende(i,:)=[global_coords{i}(1,1:2) global_coords{i}(end,1:2)];
end

fid=fopen(fullfile(pfad,[name,'_coordinates.txt']),'wt');
fprintf(fid,'Profnum\txstart\tystart\txende\tyende\tanzMarker\n');
fprintf(fid,'%d\t%f\t%f\t%f\t%f\t%d\n',[list.num startende anzmarks]');
fclose(fid);
end

