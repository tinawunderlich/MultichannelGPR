clear all
close all
clc


% Script for reading radargrams in radargrams.mat format and sorting, no other processing
% is applied in this script!
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files (please choose folder with these files):
% radargrams.mat: Radar data
% x.mat: profile coordinates
% global_coords.mat: coordinates
% t.mat: time vector
%
% Assumes parallel profiles over an area and sorts them from south to north
% (changing the order of radarrgams in the file)
% and running from west to east (changing the direction of each radargram)

plot_fig=1; % if =1: plot figure for control

save_shape=1; % if =1: save shape file with profile locations
shape_name='sorted_Profiles.shp'; % name for shape-file (will be saved in the same folder as the data)


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');


% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose folder with radargrams.mat');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            foldername=uigetdir([],'Choose folder with radargrams.mat');
        end
    else
        foldername=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


%% Load all profiles
disp('Reading data...')
load(fullfile(foldername,'radargrams.mat'));
load(fullfile(foldername,'global_coords.mat'));
load(fullfile(foldername,'t.mat'));
load(fullfile(foldername,'x.mat'));

num_old=1:length(x); % old number of radargram

% get approximate direction (fit with straight line) and direction of
% profile
disp('Determine direction and order...')
minx=1e12; % initialize for determining min x value of area
for i=1:length(x)
    p(i,:)=polyfit(global_coords{i}(:,1),global_coords{i}(:,2),1); % fit straight line
    if global_coords{i}(1,1)<=global_coords{i}(end,1)
        flag(i,1)=0; % W -> E
    else
        flag(i,1)=1; % E -> W: shall be turned
        radargrams{i}=fliplr(radargrams{i});
        global_coords{i}=flipud(global_coords{i});
        x{i}=fliplr(abs(x{i}-max(x{i})));
    end
    if min(global_coords{i}(:,1))<minx
        minx=min(global_coords{i}(:,1));
    end
end
% determining y-value at minx-line (and reset y-cutoff-value in p)
p(:,2)=p(:,1).*minx+p(:,2);

% sort from south to north
[new,num_new]=sortrows(p,2);

%% write new files
disp('Write sorted and turned files...')
if ~exist(fullfile(foldername,'sorted'))
    mkdir(fullfile(foldername,'sorted'))
end

for i=1:length(num_new)
    r{i}=radargrams{num_new(i)};
    xx{i}=x{num_new(i)};
    gc{i}=global_coords{num_new(i)};
end
radargrams=r;
x=xx;
global_coords=gc;
save(fullfile(foldername,'sorted','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(foldername,'sorted','x.mat'),'x','-v7.3');
save(fullfile(foldername,'sorted','t.mat'),'t','-v7.3');
save(fullfile(foldername,'sorted','global_coords.mat'),'global_coords','-v7.3');

% write info-file
fid=fopen(fullfile(foldername,'sorted','info_sorted.txt'),'wt');
fprintf(fid,'#old\tturned?t#new\n');
fprintf(fid,'%d\t%d\t%d\n',[num_old' flag num_new]');
fclose(fid);

% plot figure
if plot_fig==1
    col=jet(length(x));
    figure
    hold on
    for i=1:length(x)
        plot(global_coords{i}(:,1),global_coords{i}(:,2),'Color',col(i,:),'Linewidth',2)
    end
    grid on
    axis xy
    axis equal
    colorbar
    colormap(jet)
end


%% write shape-file
if save_shape==1
    disp('Write shape-file...')
    anz=1;
    for kk=1:length(x) % loop over radargrams
        if ~isempty(radargrams{kk})
            S(anz).Geometry='Line';
            S(anz).BoundingBox=[min(global_coords{kk}(:,1)) min(global_coords{kk}(:,2)); max(global_coords{kk}(:,1)) max(global_coords{kk}(:,2))];
            S(anz).X=global_coords{kk}(:,1);
            S(anz).Y=global_coords{kk}(:,2);
            S(anz).id=kk;
            anz=anz+1;
        end
    end

    % write shapefile with profile lines
    shapewrite(S,fullfile(foldername,'sorted',shape_name));
end

disp('Done!')

