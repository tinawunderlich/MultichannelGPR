clear all
close all
clc

% Script for making new depth vector with constant dz and adjusting
% radargrams: When you applied topomigration with given zmin/zmax but
% different velocity models for each profile, the radargrams will all have
% different dz and thus different lengths. With this script the depth
% vector and radargrams are interpolated onto a unique depth vector with
% constant dz.
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% requires following files (please choose folder with these files):
% radargrams.mat: Radar data
% global_coords.mat: coordinates
% t.mat: depth vector after topomigration

% Options for choosing dz
dz_opt=3;   % =1: choose the smallest dz from radargrams for all,
            % =2: choose the largest dz from radargrams for all,
            % =3: choose dz below for all 
dz=0.003; % dz in m for dz_opt=3

%% -------------- DO NOT CHANGE BELOW THIS LINE! --------------------------

% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            foldername=uigetdir([],'Choose folder with radargrams');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            foldername=uigetdir([],'Choose folder with radargrams');
        end
    else
        foldername=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end

%% Load all profiles
load(fullfile(foldername,'radargrams.mat'));
load(fullfile(foldername,'t.mat'));

zmin=min(t);
zmax=max(t);

%% determine dz
for i=1:length(radargrams)
    dzi(i)=(zmax-zmin)/(length(radargrams{i}(:,1))-1);
end
disp(['min(dz) = ',num2str(min(dzi),4),' m / max(dz) = ',num2str(max(dzi),4),' m'])

if dz_opt==1 % smallest dzi
    disp(['Choosing dz = ',num2str(min(dzi),4),' m'])
    dz=min(dzi);
elseif dz_opt==2 % largest dzi
    disp(['Choosing dz = ',num2str(max(dzi),4),' m'])
    dz=max(dzi);
else
    disp(['Choosing dz = ',num2str(dz),' m'])
end


%% interpolate on new z vector
z=zmax:-dz:zmin;
for i=1:length(radargrams)
    if ~mod(i,10)
        disp([int2str(i),'/',int2str(length(radargrams))])
    end
    temp=zeros(length(z),length(radargrams{i}(1,:)));
    for j=1:length(radargrams{i}(1,:)) % for each trace
        temp(:,j)=interp1([zmax:-dzi(i):zmin]',radargrams{i}(:,j),z');
    end
    radargrams{i}=temp;
end

t=z;

%% save data
if ~exist(fullfile(foldername,'const_dz'),'dir')
    mkdir(fullfile(foldername,'const_dz'));
end
save(fullfile(foldername,'const_dz','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(foldername,'const_dz','t.mat'),'t','-v7.3');
copyfile(fullfile(foldername,'global_coords.mat'),fullfile(foldername,'const_dz','global_coords.mat'));
copyfile(fullfile(foldername,'x.mat'),fullfile(foldername,'const_dz','x.mat'));