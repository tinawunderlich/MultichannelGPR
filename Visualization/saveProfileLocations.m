clear all
close all
clc


% Read *.pos-files and save a shape-file with the profile location for GIS.
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de


shapename='GPR_Profiles_alle.shp'; % give name for shape file

% -------------------------------------------------------------------------
% Do not change the following part!

% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with *.pos-files','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with *.pos-files');
        else
            pfad_rad=uigetdir([],'Choose folder with *.pos-files');
        end
        fileattrib('radtemp.temp','-h');
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
        fileattrib('radtemp.temp','+h');
    else
        pfad_rad=uigetdir([],'Choose folder with *.pos-files'); % path to *.pos-files-folder

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
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with *.pos-files');
        else
            pfad_rad=uigetdir([],'Choose folder with *.pos-files');
        end
    else
        pfad_rad=uigetdir([],'Choose folder with *.pos-files'); % path to *.pos-files-folder
    end

    fid=fopen('.radtemp.temp','wt');
    fprintf(fid,'%s',pfad_rad);
    fclose(fid);
end

%% Read data
disp('Reading data...')

list=dir([pfad_rad,'/*.pos']);
for i=1:length(list)
    % get profile number
    temp=extractBetween(list(i).name,'_','_');
    numbers(i,1)=str2num(temp{1});

    % get coordinates
    fid=fopen(fullfile(pfad_rad,list(i).name),'r');
    temp=textscan(fid,'%f%f%f%f','Headerlines',1);
    fclose(fid);
    data{i}=[temp{1} temp{3} temp{2} temp{4}]; % change order of x and y here!
end

%% save as shape file
disp('Saving shape file...')

anz=1;
for kk=1:length(numbers) % loop over radargrams
    if ~isempty(data{kk})
        S(anz).Geometry='Line';
        S(anz).BoundingBox=[min(data{kk}(:,2)) min(data{kk}(:,3)); max(data{kk}(:,2)) max(data{kk}(:,3))];
        S(anz).X=data{kk}(:,2);
        S(anz).Y=data{kk}(:,3);
        S(anz).id=numbers(kk);
        anz=anz+1;
    end
end

% write shapefile with profile lines
shapewrite(S,fullfile(pfad_rad,shapename));
