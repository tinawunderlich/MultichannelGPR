clear all
close all
clc

% Script for preparing the txt-file for radargram-extraction from a
% shape-file for make_Radargram.m
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires a shape-file of lines

radargram_file='Radargrams.txt'; % name of radargram-file for make_Radargram.m (will be saved in same folder as shape-file)

% for plotting:
off=4; % [m]; Offset for plotting the numbers at the starting points

%% ------------ DO NOT CHANGE FROM HERE ON! -------------------------------
warning('off');

% get folder name
% get file name
if ispc
    if exist('stemp.temp') % read last opened folder from temp.temp
        fid=fopen('stemp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.shp','Choose *.shp file with lines',fn{1}{1});
        else
            [file,folder]=uigetfile('*.shp','Choose *.shp file with lines');
        end
    else
        [file,folder]=uigetfile('*.shp','Choose *.shp file with lines');
    end
    fid=fopen('stemp.temp','wt');
    fprintf(fid,'%s',fullfile(folder,file));
    fclose(fid);
else
    if exist('.stemp.temp') % read last opened folder from temp.temp
        fid=fopen('.stemp.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.shp','Choose *.shp file with lines',fn{1}{1});
        else
            [file,folder]=uigetfile('*.shp','Choose *.shp file with lines');
        end
    else
        [file,folder]=uigetfile('*.shp','Choose *.shp file with lines');
    end
    fid=fopen('.stemp.temp','wt');
    fprintf(fid,'%s',fullfile(folder,file));
    fclose(fid);
end

% read shape-file
data=shaperead(fullfile(folder,file));

% go through file
num=1; % counter for profiles
f=1; % counter for features/id
list=[];
for i=1:length(data)
    % check if line
    if strcmp(data(i).Geometry,'Line')
        if length(data(i).X)==3 % line from one point to the other (2 points + NaN)
            list=[list; data(i).X(1) data(i).Y(1) data(i).X(2) data(i).Y(2) f num];
            num=num+1;
        else % line with multiple knots -> divide into several short profiles with same f but different num
            for j=1:length(data(i).X)-2
                list=[list; data(i).X(j) data(i).Y(j) data(i).X(j+1) data(i).Y(j+1) f num];
                num=num+1;
            end
        end
        f=f+1;
    else
        disp(['Feature ',int2str(i),': Geometry needs to be LINE.'])
    end
end

col=jet(f-1); % colors for visualization

figure('Position',[0 0 1000 1000])
plot(list(:,1),list(:,2),'*k','Linewidth',2)
hold on
plot(list(:,3),list(:,4),'*k','Linewidth',2)
for i=1:f-1
    plot([list(list(:,5)==i,1) list(list(:,5)==i,3)],[list(list(:,5)==i,2) list(list(:,5)==i,4)],'Color',col(i,:),'Linewidth',2)
    text(list(list(:,5)==i,1)-off,list(list(:,5)==i,2)+off,int2str(i))
end
axis xy
axis equal
set(gca,'Fontsize',20)
xlabel('Easting [m]')
ylabel('Northing [m]')
% save figure
saveas(gcf,fullfile(folder,'Radargrams_Figure.png'),'png');

% write txt-file for make_Radargram.m
fid=fopen(fullfile(folder,radargram_file),'wt');
fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\n',list(:,1:4)');
fclose(fid);

% write log-file
fid=fopen(fullfile(folder,'Radargram_Log.txt'),'wt');
fprintf(fid,'Log for conversion of %s',file);
fprintf(fid,' to %s.\n\n',radargram_file);
fprintf(fid,'xstart\tystart\txend\tyend\tFeature/ID\tRadargram#\n');
fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n',list');
fclose(fid);

