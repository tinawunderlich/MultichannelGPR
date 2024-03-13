clear all
close all
clc

% Script for Plotting of profiles and picking of hyperbolas -> velocity
% determination
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% requires folders Export_Import, Subfunctions, Processing, Plotting
%
% You can either select the rSlicer-folder (=original processed files will be used, Mala3D had to be run before!) or a folder containing
% radargrams.mat and corresponding files (=any radargrams, could also be from DZT or RDT import).


%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!

warning off

% get folder name
if exist('.temp.temp') % read last opened folder from temp.temp
    fid=fopen('.temp.temp','r');
    fn=textscan(fid,'%s');
    fclose(fid);
    if ~isempty(fn{1})
        foldername=uigetdir(fn{1}{1},'Choose rSlicer folder or folder with radargrams.mat');
    else
        foldername=uigetdir('Choose rSlicer folder or folder with radargrams.mat');
    end
else
    foldername=uigetdir('Choose rSlicer folder or folder with radargrams.mat'); % path to rSlicer_folder
end

% save last selected folder in file
fid=fopen('.temp.temp','w');
fprintf(fid,foldername);
fclose(fid);

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Processing'),fullfile(curFold,'Subfunctions'),fullfile(curFold,'Export_Import'),fullfile(curFold,'GUIs'));

% get all available profile numbers
outtemp=dir(foldername);
test=cellfun(@(x) endsWith(x,'.rad'),{outtemp.name});
if sum(test)>0 % rad-files found -> rSlicer folder -> original processed profiles
    out=outtemp(test); % only rad-files
    fname=out(1).name; % complete name and number.rad
    name=fname(1:end-8); % only file name without _000.rad
    for i=1:length(out)
        n=out(i).name;  % name of file
        profilelist(i,1)=str2num(n(end-6:end-4)); % profile number
    end
    
    
    % plot with gui options
    Velocity_2D_picking_plot(foldername,name,profilelist);
    
else  % folder with radargrams.mat, ...
    if exist(fullfile(foldername,'radargrams.txt'),'file')
        fid=fopen(fullfile(foldername,'radargrams.txt'),'r');
        % test for which radargrams.txt-file (DZT_Convert or Mala)
        test=fscanf(fid,'%s',1);
        if strcmp(test,'Radargrams.mat') % Mala
            temp=textscan(fid,'%f','Headerlines',1);
            temp2=textscan(fid,'%f','Headerlines',1);
            chanlist=temp{1}';
            for i=1:length(temp2)
                profilelist=temp2{1}';
            end
            name='dummyname';
        else  % DZT_Convert
            temp=textscan(fid,'%f%s%f','Headerlines',1);
            chanlist=unique(temp{3}');
            plist=temp{2};
            for i=1:length(temp{2})
                c=strsplit(plist{i},'_','CollapseDelimiters',1);
                profilelist(i)=str2num(c{end});
            end
            name=c{1};
        end
        fclose(fid);
    else
        load(fullfile(foldername,'x.mat'));
        profilelist=[1:length(x)];
        chanlist=1;
        name='dummyname';
    end

    
    % plot with gui options
    Velocity_2D_picking_plot_rad(foldername,name,profilelist,chanlist);

end


% check if figure is still open
waitfor(gcf);

% set original path
path(oldpath);
