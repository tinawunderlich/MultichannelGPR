%------------------------ Split Markergroups -------------------------------------
%
% When data was measured with marker groups, e.g. 2 zigzag lines in one
% file, split automatically between the groups. Data has to be in
% MultichannelGPR format with radargrams.mat, t.mat, x.mat,
% global_coords.mat and marker.mat
%
% Dr. Tina Wunderlich, 2024, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


clear all
close all
clc

dataplot=1; % plot radargram for controlling? 1=yes, 0=no


%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
%
% LOAD DATA

% get file names
if ispc
    if exist('dzt.temp') % read last opened folder from temp.temp
        fid=fopen('dzt.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            pfad=uigetdir([],'Choose folder with radargrams.mat');
        end
        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    else
        pfad=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder

        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
    end
else
    if exist('.dzt.temp') % read last opened folder from temp.temp
        fid=fopen('.dzt.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams.mat');
        else
            pfad=uigetdir([],'Choose folder with radargrams.mat');
        end
    else
        pfad=uigetdir([],'Choose folder with radargrams.mat'); % path to radargram-folder
    end

    fid=fopen('.dzt.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

if exist(fullfile(pfad,'radargrams.mat'))
    load(fullfile(pfad,'radargrams.mat'));
else
    disp('No data found. Please check path and try again.')
    return;
end
if exist(fullfile(pfad,'t.mat'))
    load(fullfile(pfad,'t.mat'));
else
    disp('No t.mat found. Please check path and try again.')
    return;
end
if exist(fullfile(pfad,'x.mat'))
    load(fullfile(pfad,'x.mat'));
else
    disp('No x.mat found. Please check path and try again.')
    return;
end
if exist(fullfile(pfad,'global_coords.mat'))
    load(fullfile(pfad,'global_coords.mat'));
else
    disp('No global_coords.mat found. Please check path and try again.')
    return;
end
if exist(fullfile(pfad,'marker.mat'))
    load(fullfile(pfad,'marker.mat'));
else
    disp('No marker.mat found. Please check path and try again.')
    return;
end
if exist(fullfile(pfad,'radargrams.txt'))
    fid=fopen(fullfile(pfad,'radargrams.txt'),'r');
    temp=textscan(fid,'%s%s%s','Headerlines',1);
    fclose(fid);
    name=temp{2}; % names of files
end

count=1;
for i=1:length(marker)
    if ~exist('name','var') || length(name)<i
        name{i}=int2str(i); % profile number (=name)
    end
    disp(name{i})

    % find marker indices
    marks=find(marker{i}==1);
    % differences between markers (in traces):
    d=diff(marks);
    % find large gaps between markers:
    gaps=find(d>=mean(d)+mean(d)/100*50); % index of marker

    if dataplot==1
        figure
        subplot(2,3,1:3)
        imagesc([1:length(x{i})],t,radargrams{i})
        hold on
        colorbar
        for ii=1:length(marks)
            plot([marks(ii) marks(ii)],[0 10],'k')
        end
        colormap(flipud(gray))
        xlabel('Trace number')
        ylabel('t [ns]')
    end

    if ~isempty(gaps) % marker groups found...
        if length(gaps)==1
            % first group
            data{count}=radargrams{i}(:,1:marks(gaps)+10);
            gc{count}=global_coords{i}(1:marks(gaps)+10,:);
            xx{count}=x{i}(1:marks(gaps)+10);
            ma{count}=marker{i}(1:marks(gaps)+10);

            if dataplot==1
                subplot(2,3,4)
                imagesc([1:length(xx{count})],t,data{count})
                hold on
                colorbar
                m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end
            count=count+1;

            % last group
            data{count}=radargrams{i}(:,marks(gaps)+11:end);
            gc{count}=global_coords{i}(marks(gaps)+11:end,:);
            xx{count}=x{i}(marks(gaps)+11:end);
            ma{count}=marker{i}(marks(gaps)+11:end);

            part{i}=1:2;

            if dataplot==1
                subplot(2,3,5)
                imagesc([1:length(xx{count})],t,data{count})
                hold on
                colorbar
                m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end
            count=count+1;

        elseif length(gaps)==2
            % first group
            data{count}=radargrams{i}(:,1:marks(gaps(1))+10);
            gc{count}=global_coords{i}(1:marks(gaps(1))+10,:);
            xx{count}=x{i}(1:marks(gaps(1))+10);
            ma{count}=marker{i}(1:marks(gaps(1))+10);

            if dataplot==1
                subplot(2,3,4)
                imagesc([1:length(xx{count})],t,data{count})
                hold on
                colorbar
                m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end
            count=count+1;

            % second group
            data{count}=radargrams{i}(:,marks(gaps(1))+11:marks(gaps(2))+10);
            gc{count}=global_coords{i}(marks(gaps(1))+11:marks(gaps(2))+10,:);
            xx{count}=x{i}(marks(gaps(1))+11:marks(gaps(2))+10);
            ma{count}=marker{i}(marks(gaps(1))+11:marks(gaps(2))+10);

            if dataplot==1
                subplot(2,3,5)
                imagesc([1:length(xx{count})],t,data{count})
                hold on
                colorbar
                m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end

            count=count+1;

            % last group
            data{count}=radargrams{i}(:,marks(gaps(2))+11:end);
            gc{count}=global_coords{i}(marks(gaps(2))+11:end,:);
            xx{count}=x{i}(marks(gaps(2))+11:end);
            ma{count}=marker{i}(marks(gaps(2))+11:end);

            part{i}=1:3;

            if dataplot==1
                subplot(2,3,6)
                imagesc([1:length(xx{count})],t,data{count})
                hold on
                colorbar
                m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
                colormap(flipud(gray))
                xlabel('Trace number')
                ylabel('t [ns]')
            end
            count=count+1;
        end
    else
        part{i}=1;

        data{count}=radargrams{i};
        gc{count}=global_coords{i};
        xx{count}=x{i};
        ma{count}=marker{i};

        if dataplot==1
            subplot(2,3,4:6)
            imagesc([1:length(xx{count})],t,data{count})
            hold on
            colorbar
            m=find(ma{count}==1);
                for ii=1:length(m)
                    plot([m(ii) m(ii)],[0 10],'k')
                end
            colormap(flipud(gray))
            xlabel('Trace number')
            ylabel('t [ns]')
        end
        count=count+1;
    end
end

% set data in correct variables:
radargrams=data;
global_coords=gc;
x=xx;
marker=ma;

% save data:
if ~exist(fullfile(pfad,'split_markergroups'),'dir')
    mkdir(fullfile(pfad,'split_markergroups'));
end
save(fullfile(pfad,'split_markergroups','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(pfad,'split_markergroups','t.mat'),'t','-v7.3');
save(fullfile(pfad,'split_markergroups','x.mat'),'x','-v7.3');
save(fullfile(pfad,'split_markergroups','global_coords.mat'),'global_coords','-v7.3');
save(fullfile(pfad,'split_markergroups','marker.mat'),'marker','-v7.3');

% save txt-file
fid=fopen(fullfile(pfad,'split_markergroups','radargrams_splitMarkergroups.txt'),'wt');
fprintf(fid,'Nr.\tName\tPart\n');
c=1;
for i=1:length(name)
    for j=1:length(part{i})
        fprintf(fid,'%d\t',c);
        fprintf(fid,'%s\t',name{i});
        fprintf(fid,'%d\n',part{i}(j));
        c=c+1;
    end
end
fclose(fid);