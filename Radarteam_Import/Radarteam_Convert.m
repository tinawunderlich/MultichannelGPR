clear all
close all
clc

%------------------------ Radarteam_Convert -------------------------------------
%
% Converts sgy-files of Radarteam equipment to mat-files (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, September 2023, tina.wunderlich@ifg.uni-kiel.de, CAU
% Kiel
%
%--------------------------------------------------------------------------


dataplot=1; % plot radargram for controlling? 1=yes, 0=no

% Coordinate settings:
filter_coords=1; % if =1: coordinates are filtered and traces at
% approx. the same position will be removed (e.g. at the beginning or
% end of the profile): a local difference
% over m coordinates for x and y is calculated separately. If this
% absolute local difference is smaller than nstd*std(diff) it will be removed.
m=15; % number of samples for local difference calculation
nstd=1;  % nstd*standard deviation as valid range
show_filter_plot=1; % if =1: show plot for quality control of filtering
coords_smooth=0; % if=1: smooth coordinates along profile
nsmooth=15;  % number of traces for smoothing of coords
utmzone=33; % give UTM-zone for conversion from WGS84 to UTM

offsetGPS_X=0; % Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=0; % Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)
offsetGPS_Z=0.4; % Offset between GPS and ground surface (always positive) (use only if not corrected inside the GPS acquisition)

% Options for calculating inline coordinates for each trace:
coords_opt=2;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
% =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)



%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
%
% LOAD DATA

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));


% get file names
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
            pfad=uigetdir(fn{1}{1},'Choose folder with sgy-file(s)');
        else
            pfad=uigetdir([],'Choose folder with sgy-file(s)');
        end
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with sgy-file(s)'); % path to radargram-folder

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
            pfad=uigetdir(fn{1}{1},'Choose folder with sgy-file(s)');
        else
            pfad=uigetdir([],'Choose folder with sgy-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with sgy-file(s)'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*.sgy')); % Radarteam format
% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end

% go through list and read data
for i=1:length(list)
    disp(['File ',int2str(i),' of ',int2str(length(list))])

    filename=list(i).name;
    [~,name{i},ext]=fileparts(filename);

    [traces,coords,headers{i}]=readRadarteam(fullfile(pfad,filename),filter_coords,m,nstd,show_filter_plot,coords_smooth,nsmooth,utmzone);

    if ~isempty(traces)
        flag_ok(i)=1;
        % correct GPS-antenna offset x y:
        if offsetGPS_X~=0 || offsetGPS_Y~=0
            anz2=round(0.5/mean(sqrt(diff(coords(:,2)).^2+diff(coords(:,3)).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
            if anz2/2==round(anz2/2)
                anz2=anz2+1; % make odd
            end
            for ii=1:length(coords(:,2))-anz2
                dist=sqrt((coords(ii,2)-coords(ii+anz2,2))^2+(coords(ii,3)-coords(ii+anz2,3))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[coords(ii,2) coords(ii,3); coords(ii+anz2,2) coords(ii+anz2,3)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
            end
            anz1=anz2;
            % calculation for the last few traces
            anz2=anz2-1;
            for ii=length(coords(:,2))-anz1+1:length(coords(:,2))-1
                dist=sqrt((coords(ii,2)-coords(ii+anz2,2))^2+(coords(ii,3)-coords(ii+anz2,3))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[coords(ii,2) coords(ii,3); coords(ii+anz2,2) coords(ii+anz2,3)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
                anz2=anz2-1;
            end
            % extrapolate for last trace
            coords(end,2)=interp1([1 2],[coords(end-2,2) coords(end-1,2)],3,'linear','extrap');
            coords(end,3)=interp1([1 2],[coords(end-2,3) coords(end-1,3)],3,'linear','extrap');
        end
        % correct GPS height
        coords(:,4)=coords(:,4)-offsetGPS_Z;


        % set in variables:
        radargrams{i}=traces;
        global_coords{i}=coords(:,2:4);
        xtemp=coords(:,2);
        ytemp=coords(:,3);
        % calculate profile coordinates
        if coords_opt==1 % difference to beginning
            x{i}=[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)'];
        elseif coords_opt==2 % cumulative sum
            x{i}=[0; reshape(cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)'),length(xtemp)-1,1)];
        end
        t=double(0:headers{i}.dt:headers{i}.dt*(headers{i}.ns-1)).*1e-3; % in ns

        if dataplot==1
            figure
            imagesc(x{i},t,radargrams{i})
            title(name{i})
            xlabel('x [m]')
            ylabel('t [ns]')
            colormap(flipud(gray))
        end
    else
        flag_ok(i)=0; % empty traces because of bad coordinates
    end

end

%------------------------------ WRITE DATA --------------------------------


% delete empty cells
ok=~cellfun('isempty',radargrams);
radargrams=radargrams(ok);
x=x(ok);
global_coords=global_coords(ok);
headers=headers(ok);

disp('Start saving of mat-files...')


% make new folder
if ~exist(fullfile(pfad,'mat'),'dir')
    mkdir(fullfile(pfad,'mat'));
end
save(fullfile(pfad,'mat','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(pfad,'mat','t.mat'),'t','-v7.3');
save(fullfile(pfad,'mat','x.mat'),'x','-v7.3');
save(fullfile(pfad,'mat','global_coords.mat'),'global_coords','-v7.3');
save(fullfile(pfad,'mat','h.mat'),'headers','-v7.3');

fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
fprintf(fid,'Nr.\tName\n');
anz=1;
for i=1:length(name)
    if flag_ok(i)==1
        fprintf(fid,'%d\t',anz);
        fprintf(fid,'%s\n',[name{i},'.sgy']);
        anz=anz+1;
    end
end
fclose(fid);

disp('   Finished!')

% set original path
path(oldpath);
