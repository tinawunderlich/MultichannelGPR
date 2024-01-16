clear all
close all
clc

%% Read all DZG-files in folder and produce profile plan
%
% Dr. Tina Wunderlich, CAU Kiel, 2022, tina.wunderlich@ifg.uni-kiel.de
%%

option=2;   % Option=1: plot arrow from start to endpoint of profile
% Option=2: plot original coordinates
% Option=3: plot smoothed profile coordinates (Matlab only! For Octave the same as option=2.)

convert2utm=1; % convert WGS Lat/Long to UTM (=1 if measured with e.g. Stonex-GPS)
zone=32; % if convert2utm==1 -> give UTM-zone

checkCoords=0; % normally use =0; checkCoords=1: Check for large coordinate jumps and remove them (occurs sometimes for measurements with total station)

legend_on=1; % =1 if legend, =0 no legend
map_xlim=[]; % give x-limits of map [start end] (leave empty for automatic determination)
map_ylim=[]; % give y-limits of map [start end] (leave empty for automatic determination)
xtick=[]; % spacing between x-ticks for map in m (leave empty for automatic determination)
ytick=[]; % spacing between y-ticks for map in m (leave empty for automatic determination)

% offsets between arrow head and profile number [m]
offx=0.5;
offy=0.2;

%-------------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));

test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end

% get file names
if ispc && matlab==1
    if exist('dzg.temp') % read last opened folder from temp.temp
        fid=fopen('dzg.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZG-file(s)');
        else
            pfad=uigetdir([],'Choose folder with DZG-file(s)');
        end
        fileattrib('dzg.temp','-h');
        fid=fopen('dzg.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('dzt.temp','+h');
    else
        pfad=uigetdir([],'Choose folder with DZG-file(s)'); % path to radargram-folder
        
        fid=fopen('dzg.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('dzg.temp','+h');
    end
else
    if exist('.dzg.temp') % read last opened folder from temp.temp
        fid=fopen('.dzg.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZG-file(s)');
        else
            if matlab==1
                pfad=uigetdir([],'Choose folder with DZG-file(s)'); % path to radargram-folder
            else
                pfad=uigetdir('','Choose folder with DZG-file(s)');
            end
        end
    else
        if matlab==1
            pfad=uigetdir([],'Choose folder with DZG-file(s)'); % path to radargram-folder
        else
            pfad=uigetdir('','Choose folder with DZG-file(s)');
        end
    end
    
    fid=fopen('.dzg.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*.DZG')); % GSSI format

% check for names starting with . and ending with _backup.DZG
ii=1;
while ii<=length(list)
    if strcmp(list(ii).name(1),'.') || strcmp(list(ii).name(end-10:end),'_backup.DZG')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end

fm=figure('Visible','on');
hold on

leg=[];
for i=1:length(list)
    % name zerlegen
    [pathstr,fname,ext] = fileparts(list(i).name);
    
    if ~exist([pfad,'/',fname,'_backup.DZG'])  % if backup already exists, some changes have been made before, so do not do the changes again.
        % PART 1: check DZG-file if corrupted (sometimes if measured with Stonex-GPS)
        copyfile([pfad,'/',fname,'.DZG'],[pfad,'/',fname,'_backup.DZG']); % Save backup file
        
        fileIDorig = fopen([pfad,'/',fname,'_backup.DZG'],'r');
        fileIDnew = fopen([pfad,'/',fname,'.DZG'],'wt');
        
        dataorig = textscan(fileIDorig,'%s','Delimiter','\n');
        dataorig = dataorig{1};
        ndata = length(dataorig);
        for j = 1:ndata
            
            tmpline = strsplit(dataorig{j},',','CollapseDelimiters',false);
            nentry = length(tmpline);
            if strcmp(tmpline{1},'$GSSIS') && nentry == 3 && ...
                    ~isempty(str2double(tmpline{2})) && ...
                    ~isempty(str2double(tmpline{3}))
                fprintf(fileIDnew,'%s\n',dataorig{j});
            elseif strcmp(tmpline{1},'$GPGGA') && nentry == 15
                fprintf(fileIDnew,'%s\n',dataorig{j});
            elseif isempty(tmpline{1})
                fprintf(fileIDnew,'%s','');
            end
        end
        fclose(fileIDorig);
        fclose(fileIDnew);
        s1=dir([pfad,'/',fname,'_backup.DZG']);
        s2=dir([pfad,'/',fname,'.DZG']);
        if s1.bytes==s2.bytes  % if same size, delete backup, because no changes
            delete([pfad,'/',fname,'_backup.DZG']);
        end
        % file.DZG is now without GPSerrors
    end
    
    % PART 2: Read DZG-file
    fid=fopen([pfad,'/',fname,'.DZG']);
    temp=textscan(fid,'%s');
    fclose(fid);
    
    if ~isempty(temp{1})
        
        count = 1;
        % loop over lines of file
        for j=1:length(temp{1})
            
            % the first line contains the trace number for the coordinate
            if strfind(temp{1}{j},'$GSSIS')

                if matlab==1
                    GSSILine=split(temp{1}{j},',');
                else
                    GSSILine=strsplit(temp{1}{j},',');
                end
                
                if j < length(temp{1})
                    % the coordinates are found in the second line
                    if strfind(temp{1}{j+1},'$GPGGA')
                        if matlab==1
                            GPGGALine=split(temp{1}{j+1},',');
                        else
                            GPGGALine=strsplit(temp{1}{j+1},',');
                        end
                        % extract the coordinate information
                        x(count) = str2double(GPGGALine{5});
                        y(count) = str2double(GPGGALine{3});
                        z(count) = str2double(GPGGALine{10});
                        
                        % extract the trace number
                        trn(count)=str2double(GSSILine{2});
                        count = count + 1;
                    end
                end
            end
            
        end
        [trnUnique, ia, ~] = unique(trn);
        
        xyz{i} = [ trnUnique', x(ia)', y(ia)', z(ia)', zeros(size(x(ia)'))];
        
        % check for large coordinate jumps:
        if checkCoords==1
            mdx=median(abs(diff(xyz{i}(:,2)))); % median difference in x-coordinates
            mdy=median(abs(diff(xyz{i}(:,3)))); % median difference in y-coordinates
            scal=5;
            for j=2:length(xyz{i}(:,1))
                if abs(diff(xyz{i}(j-1:j,2)))>=scal*mdx || abs(diff(xyz{i}(j-1:j,3)))>=scal*mdy
                    xyz{i}(j:end,5)=1;
                    break;
                end
            end
            xyz{i}(xyz{i}(:,5)==1,:)=[]; % delete these lines
        end
        % delete line if coordinate is NaN
        weg=isnan(xyz{i}(:,2)) | isnan(xyz{i}(:,3)) | isnan(xyz{i}(:,4));
        xyz{i}(weg,:)=[];
        
        clear trn;
        clear x;
        clear y;
        clear z;
        
        % Lat/Lon to UTM conversion
        if convert2utm==1
            lontemp=num2str(xyz{i}(:,2),'%.8f');
            lattemp=num2str(xyz{i}(:,3),'%.8f');
            for j=1:length(lattemp(:,1))
                temp=strsplit(lattemp(j,:),'.');
                lat(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
                temp=strsplit(lontemp(j,:),'.');
                lon(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
            end
            [xneu,yneu]=wgs2utm(lat,lon,zone,'N');
            xyz{i}(:,2)=xneu';
            xyz{i}(:,3)=yneu';
            
            clear lat;
            clear lon;
        end
        
        if matlab==0
            % subtract coordinate offset
            coordstrx=num2str(round(min(xyz{i}(:,2))));
            coordstry=num2str(round(min(xyz{i}(:,3))));
            coordstrx(end-2:end)='000';
            coordstry(end-2:end)='000';
            xyz{i}(:,2)=xyz{i}(:,2)-str2num(coordstrx);
            xyz{i}(:,3)=xyz{i}(:,3)-str2num(coordstry);
        end
        
        leg=[leg; {num2str(i)}];
        % plot profile in map
        if option==1
            if matlab==1
                quiver(xyz{i}(1,2),xyz{i}(1,3),xyz{i}(end,2)-xyz{i}(1,2),xyz{i}(end,3)-xyz{i}(1,3),0,'Linewidth',2);
            else
                quiver(xyz{i}(1,2),xyz{i}(1,3),xyz{i}(end,2)-xyz{i}(1,2),xyz{i}(end,3)-xyz{i}(1,3),0,'Linewidth',2,'MaxHeadSize',0.05);
            end
        elseif option==2
            p=plot(xyz{i}(:,2),xyz{i}(:,3),'Linewidth',2);
            pc{i}=get(p,'Color');
        else
            if matlab==0
                p=plot(xyz{i}(:,2),xyz{i}(:,3),'Linewidth',2);
            else
                p=plot(smooth(xyz{i}(:,2)),smooth(xyz{i}(:,3)),'Linewidth',2);
            end
            pc{i}=get(p,'Color');
        end
    end
end
for i=1:length(list)
    % plot additional things in map
    if option==2 || option==3
        if matlab==1
            quiver(xyz{i}(1,2),xyz{i}(1,3),xyz{i}(end,2)-xyz{i}(1,2),xyz{i}(end,3)-xyz{i}(1,3),0,'Color',pc{i},'Linewidth',1,'LineStyle','-','ShowArrowHead','on');
        else
            quiver(xyz{i}(1,2),xyz{i}(1,3),xyz{i}(end,2)-xyz{i}(1,2),xyz{i}(end,3)-xyz{i}(1,3),0,'Color',pc{i},'MaxHeadSize',0.05,'Linewidth',1,'LineStyle','-','ShowArrowHead','on');
        end
    end
    text(xyz{i}(end,2)+offx,xyz{i}(end,3)+offy,int2str(i))
end

if legend_on==1
    legend(leg,'Location','eastoutside')
end
grid on
if ~isempty(map_xlim)
    set(gca,'XLim',map_xlim)
else
    map_xlim=get(gca,'XLim');
end
if ~isempty(map_ylim)
    set(gca,'yLim',map_ylim)
else
    map_ylim=get(gca,'YLim');
end
if ~isempty(xtick)
    set(gca,'XTick',[map_xlim(1):xtick:map_xlim(2)])
end
if ~isempty(ytick)
    set(gca,'YTick',[map_ylim(1):ytick:map_ylim(2)])
end
xticks=get(gca,'XTick');
yticks=get(gca,'YTick');
for i=1:length(xticks)
    xticklabels{i}=sprintf('%d',xticks(i));
end
for i=1:length(yticks)
    yticklabels{i}=sprintf('%d',yticks(i));
end
xlabel('Easting [m]')
ylabel('Northing [m]')
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'XTickLabelRotation',0,'YTick',yticks,'YTickLabel',yticklabels,'DataAspectRatio',[1 1 1])

if matlab==0
    title(['Coordinate offset: x=',coordstrx,'m, y=',coordstry,'m'])
end

if ~exist(fullfile(pfad,'Figures'),'dir')
    mkdir(fullfile(pfad,'Figures'));
end

saveas(fm,fullfile(pfad,'Figures','Map.eps'),'epsc')
saveas(fm,fullfile(pfad,'Figures','Map.pdf'),'pdf')
saveas(fm,fullfile(pfad,'Figures','Map.png'),'png')



% set original path
path(oldpath);