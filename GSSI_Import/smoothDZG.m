clear all
close all
clc

%% Takes *.DZG-files from GSSI-GPS-Data and smoothes the profile if there are coordinate jumps due to bad GPS signal
% =If you pulled a straight profile like this: ...........................
% but your profile now looks like: .....       .............   ...........
%                                       .......             ...
% 
% A straight line will be fitted with RANSAC to most of the coordinate
% points and all points with distance<maxDist will be projected on this
% profile. You will receive a straight line in the end.
% 
% Attention: This script will read your *.DZG-files, move them to
% a folder DZG_old and write new *.DZG files. These new corrected files are then
% read in with DZT_Convert.m. So run this script before DZT_Convert!
%
% Dr. Tina Wunderlich, August 2022, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------

convert2utm=1; % if you measured with a stonex use =1 
% (Attention: Your corrected DZG-Files will be in UTM already, so use convert2utm=0 in DZT_Convert!)
zone=34; % UTM zone

% do you want to plot the coordinates for control?
plotflag=1; % =1: yes, =0: no

% Settings for RANSAC fit:
maxNum=1/3; % maxNum of all data points will be used as sample for fitting
maxDist=0.05; % max allowable distance for inliers [m]

%--------------------------------------------------------------------------
%% DO NOT CHANGE FROM HERE ON!

warning('off');

test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end

% get file names
if ispc
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
        fileattrib('dzg.temp','+h');
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
            pfad=uigetdir([],'Choose folder with DZG-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with DZG-file(s)'); % path to radargram-folder
    end
    
    fid=fopen('.dzg.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*.DZG')); % GSSI format
% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end


% RANSAC functions
fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2); % distance evaluation function

if plotflag==1
    figure(1)
    hold on
    figure(2)
    hold on
end
for i=1:length(list)
    % name of DZG-file
    name=list(i).name;
    
    fid=fopen(fullfile(pfad,name));
    temp=textscan(fid,'%s');
    fclose(fid);
    
    if ~isempty(temp{1})
        count = 1;
        % loop over lines of file
        for j=1:length(temp{1})
            
            % the first line contains the trace number for the coordinate
            if strfind(temp{1}{j},'$GSSIS')
                if matlab
                    GSSILine=split(temp{1}{j},',');
                else
                    GSSILine=strsplit(temp{1}{j},',');
                end
                
                
                if j < length(temp{1})
                    % the coordinates are found in the second line
                    if strfind(temp{1}{j+1},'$GPGGA')
                        if matlab
                            GPGGALine=split(temp{1}{j+1},',');
                        else
                            GPGGALine=strsplit(temp{1}{j+1},',');
                        end
                        % ectract the coordinate information
                        x(count) = str2double(GPGGALine{5});
                        y(count) = str2double(GPGGALine{3});
                        z(count) = str2double(GPGGALine{10});
                        
                        % extract the time of each sample
                        if matlab
                            t(count) = datetime(GPGGALine{2},'InputFormat','HHmmss.SS');
                        end
                        % extract the trace number
                        trn(count)=str2double(GSSILine{2});
                        count = count + 1;
                    end
                end
            end
        end
    end
    if ~exist(fullfile(pfad,'DZG_Orig'))
        mkdir(fullfile(pfad,'DZG_Orig'));
    end
    if ~exist(fullfile(pfad,'DZG_Orig',name))
        copyfile(fullfile(pfad,name),fullfile(pfad,'DZG_Orig',name)); % save Original file in new folder
    end
    
    % Lat/Lon to UTM conversion
    if convert2utm==1
        lontemp=num2str(x','%.8f');
        lattemp=num2str(y','%.8f');
        for j=1:length(lattemp(:,1))
            temp2=strsplit(lattemp(j,:),'.');
            lat(j)=str2num(temp2{1}(1:end-2))+str2num([temp2{1}(end-1:end),'.',temp2{2}])/60;
            temp2=strsplit(lontemp(j,:),'.');
            lon(j)=str2num(temp2{1}(1:end-2))+str2num([temp2{1}(end-1:end),'.',temp2{2}])/60;
        end
        [xneu,yneu]=wgs2utm(lat,lon,zone,'N');
        x=xneu;
        y=yneu;
        
        clear lat;
        clear lon;
    end
    
    % evtl. sorting
    data=sortrows([trn' x' y' z'],1); % sorting for trace number
    
    % RANSAC fit
    points=data(:,2:3);
    sampleSize=round(length(trn)*maxNum);  % number of points to sample per trial
    [modelRANSAC, inlierIdx, status] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDist);
    
    % project points onto line:
    perpSlope=-1/modelRANSAC(1); % Slope of perpendicular line
    for k=1:length(points(:,1))
        if inlierIdx(k)==1
            yInt = -perpSlope * points(k,1) + points(k,2);
            xIntersection(k) = (yInt - modelRANSAC(2)) / (modelRANSAC(1) - perpSlope);
            yIntersection(k) = perpSlope * xIntersection(k) + yInt;
        else
            xIntersection(k)=NaN;
            yIntersection(k)=NaN;
        end
    end
    
    if plotflag==1
%     figure
%     subplot(3,1,1)
%     plot(points(:,1),points(:,2),'b*')
%     hold on
%     plot(points(inlierIdx,1),points(inlierIdx,2),'r*')
%     plot(xIntersection(inlierIdx),yIntersection(inlierIdx),'*k')
%     plot(points(inlierIdx,1),modelRANSAC(1).*points(inlierIdx,1)+modelRANSAC(2),'r','Linewidth',1)
%     axis equal
%     legend('raw data','inlier','points on straight line','line fit')
%     title(['Profile ',int2str(i)])
%     subplot(3,1,2)
%     plot(data(:,1),points(:,1),'*')
%     hold on
%     plot(data(inlierIdx,1),points(inlierIdx,1),'r*')
%     title('x')
%     subplot(3,1,3)
%     plot(data(:,1),points(:,2),'*')
%     hold on
%     plot(data(inlierIdx,1),points(inlierIdx,2),'r*')
%     title('y')
    
        figure(1)
        plot(points(:,1),points(:,2),'*')

        figure(2)
        plot(xIntersection,yIntersection,'*')
    end
    
    clear trn
    clear x
    clear y
    clear z
    
    
    
    % write new DZG-file, only use inlier points projected on line
    fid=fopen(fullfile(pfad,name),'wt');
    count = 1;
    % loop over lines of file
    for j=1:length(temp{1})
        % the first line contains the trace number for the coordinate
        if strfind(temp{1}{j},'$GSSIS')
            if matlab
                GSSILine=split(temp{1}{j},',');
            else
                GSSILine=strsplit(temp{1}{j},',');
            end
            
            if j < length(temp{1})
                % the coordinates are found in the second line
                if strfind(temp{1}{j+1},'$GPGGA')  % if GGA-line
                    if matlab
                        GL=split(temp{1}{j+1},',');
                    else
                        GL=strsplit(temp{1}{j+1},',');
                    end
                    
                    % extract the trace number
                    trn=str2double(GSSILine{2});
                    
                    test=find(trn==data(inlierIdx,1),1,'first');
                    if ~isempty(test) && ~isnan(xIntersection(test))
                        % its a good trace (inlier)
                        fprintf(fid,'%s\n',temp{1}{j}); % GSSI line
                        fprintf(fid,'%s',[GL{1},',',GL{2},',']);
                        fprintf(fid,'%14.10f,N,%14.10f,E,',[yIntersection(test) xIntersection(test)]);
                        fprintf(fid,'%s',[GL{7},',', GL{8},',', GL{9},',', GL{10},',M,', GL{12},',M,', GL{14},',', GL{15}]);
                        fprintf(fid,'\n');
                    end
                end
            end
        end
    end
    fclose(fid);
    
    clear xIntersection
    clear yIntersection
end

if plotflag==1
    figure(1)
    axis equal
    grid on
    title('Raw data')
    
    figure(2)
    axis equal
    grid on
    title('Corrected data')
end