clear all
close all
clc


% Read Mala Multichannel *.pos-files and save a geoJSON-file with the profile location for GIS.
%
% Dr. Tina Wunderlich, CAU Kiel 2023-2025, tina.wunderlich@ifg.uni-kiel.de


% save profile coordinates as geoJSON file (-> profile lines for QGIS)
save_geoJSON=1; % yes=1, no=0
jsonname='GPR_Profiles.json'; % give name for geoJSON file
epsg=32631; % EPGS code for CRS

% -------------------------------------------------------------------------
% Do not change the following part!

% get folder name - RADARGRAMS
if ~ispc; menu('Choose folder with *.pos-files','OK'); end
if ispc
    if exist('radtemp.temp') % read last opened folder from temp.temp
        fid=fopen('radtemp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad_rad=uigetdir(fn{1}{1},'Choose folder with *.pos-files');
        else
            pfad_rad=uigetdir([],'Choose folder with *.pos-files');
        end
        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
    else
        pfad_rad=uigetdir([],'Choose folder with *.pos-files'); % path to *.pos-files-folder

        fid=fopen('radtemp.temp','wt');
        fprintf(fid,'%s',pfad_rad);
        fclose(fid);
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
    if ~startsWith(list(i).name,'.')
    % get profile number
    temp=extractBetween(list(i).name,'_','_');
    numbers(i,1)=str2num(temp{1});

    % get coordinates
    fid=fopen(fullfile(pfad_rad,list(i).name),'r');
    temp=textscan(fid,'%f%f%f%f','Headerlines',1);
    fclose(fid);
    data{i}=[temp{1} temp{3} temp{2} temp{4}]; % change order of x and y here!
    end
end

%% save as shape file
disp('Saving shape file...')

s.type='FeatureCollection';
s.features=[];
s.crs=[];
s.crs.type='name';
s.crs.properties.name=['EPSG:',num2str(epsg)];

for kk=1:length(numbers) % loop over radargrams
    s.features(kk).type='Feature';
    s.features(kk).geometry.type='LineString';
    s.features(kk).geometry.coordinates=[data{kk}(:,2) data{kk}(:,3)];
    s.features(kk).ids=numbers(kk); % line number
end

fid=fopen(fullfile(pfad_rad,'Figures',jsonname),'w');
fwrite(fid,jsonencode(s));
fclose(fid);
