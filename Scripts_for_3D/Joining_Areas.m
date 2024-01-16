clear all
close all
clc

% Joining of rectangles with data into big matrix (probably not possible
% for very large areas due to computer memory!)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Requires binned 3D data in folders 3D_Grid_R*, if asked give rSlicer
% folder (folder with 3D_Grid_R*-folders inside)


numbers=1:2;   % Numbers of rectangles for joining

% use processed data (if =1, then use data in /processed folder)
proc=0;

% only use one of the following (or set both to 0!):
sub_mean=0;     % if ==1: subtract mean for each rectangle for each time sample (to reduce offsets)
sub_median=0;   % if ==1: subtract median for each rectangle for each time sample (to reduce offsets)

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');

% select rSlicer-folder
% get folder name
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pathname=uigetdir([],'Choose rSlicer folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pathname);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pathname=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pathname);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pathname=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            pathname=uigetdir([],'Choose rSlicer folder');
        end
    else
        pathname=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pathname);
    fclose(fid);
end


for i=1:length(numbers)
    % read only x and y
    if proc==0
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'x.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'y.mat'));
    else
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','x.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','y.mat'));
    end

    if i==1
        if proc==0
            load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'t.mat'));
        else
            load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','t.mat'));
        end
        dx=x(1,2)-x(1,1);   % grid interval in m
        dt=t(2)-t(1);   % time interval in ns
        ns=length(t);   % number of samples
    end
    
    xmin(i)=min(x(1,:));
    xmax(i)=max(x(1,:));
    ymin(i)=min(y(:,1));
    ymax(i)=max(y(:,1));
    
    % save grids
    xx{i}=x;
    yy{i}=y;    
end

clear x y;

%%% make new grid
[xgrid,ygrid]=meshgrid([min(xmin):dx:max(xmax)],[min(ymin):dx:max(ymax)]);  % including border!
alles=NaN(length(ygrid(:,1)),length(xgrid(1,:)),ns);  % initialize big data matrix
mask_all=zeros(length(ygrid(:,1)),length(xgrid(1,:)));  % initialize big mask matrix
topo=NaN(length(ygrid(:,1)),length(xgrid(1,:)));  % initialize big topo matrix

h=waitbar(0,'Read data');
a=1;
for i=1:length(numbers)
    % load data and mask of current rectangle
    if proc==1
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','data.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','mask.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'processed','z.mat'));
    else
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'data.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'mask.mat'));
        load(fullfile(pathname,['3D_Grid_R',int2str(numbers(i))],'z.mat'));
    end
    
    % find indices of rectangle in large area (and adjust if necessary)
    startind=[find(ygrid(:,1)>=ymin(i),1) find(xgrid(1,:)>=xmin(i),1)];
    endind=[max(find(ygrid(:,1)<=ymax(i))) max(find(xgrid(1,:)<=xmax(i)))];
    if length(startind(1):endind(1))<=length(data(:,1,1))   % if number of rows different
        anzrows=length(data(:,1,1))-length(startind(1):endind(1));  % difference in rows
        if endind(1)+anzrows<=length(data(:,1,1))
            endind(1)=endind(1)+anzrows;
        else
            startind(1)=startind(1)-anzrows;
        end
    end
    if length(startind(2):endind(2))<=length(data(1,:,1))     % if number of columns different
        anzcol=length(data(1,:,1))-length(startind(2):endind(2));
        if endind(2)+anzcol<=length(data(1,:,1))
            endind(2)=endind(2)+anzcol;
        else
            startind(2)=startind(2)-anzcol;
        end
    end 

    % set data for each time sample in big matrix
    for j=1:ns
        % apply mask on data
        temp=data(:,:,j);
        temp(mask==0)=NaN;
        data(:,:,j)=temp;
        
        % set data in big matrix
        if sub_mean==1
            alles(startind(1):endind(1), startind(2):endind(2),j)=data(:,:,j)-mean(mean(data(:,:,j)));
        elseif sub_median==1
            alles(startind(1):endind(1), startind(2):endind(2),j)=data(:,:,j)-median(median(data(:,:,j)));
        else
            alles(startind(1):endind(1), startind(2):endind(2),j)=data(:,:,j);
        end
        
        
        waitbar(a/(length(numbers)*ns));
        a=a+1;
    end
    
    % set topo-matrix
    topo(startind(1):endind(1),startind(2):endind(2))=z;
end
close(h);
mask_all(~isnan(alles(:,:,1)))=1;

% renaming variables
data=alles;
mask=mask_all;


disp('Saving data in new folder')

if proc==0
    if ~exist(fullfile(pathname,'CompleteArea'),'dir')
        mkdir(fullfile(pathname,'CompleteArea'))
    end
    save(fullfile(pathname,'CompleteArea','xgrid.mat'),'xgrid','-v7.3');
    save(fullfile(pathname,'CompleteArea','ygrid.mat'),'ygrid','-v7.3');
    save(fullfile(pathname,'CompleteArea','data.mat'),'data','-v7.3');
    save(fullfile(pathname,'CompleteArea','mask.mat'),'mask','-v7.3');
    save(fullfile(pathname,'CompleteArea','topo.mat'),'topo','-v7.3');
else
    if ~exist(fullfile(pathname,'CompleteArea'),'dir')
        mkdir(fullfile(pathname,'CompleteArea'))
    end
    if ~exist(fullfile(pathname,'CompleteArea','processed'),'dir')
        mkdir(fullfile(pathname,'CompleteArea','processed'))
    end
    save(fullfile(pathname,'CompleteArea','processed','xgrid.mat'),'xgrid','-v7.3');
    save(fullfile(pathname,'CompleteArea','processed','ygrid.mat'),'ygrid','-v7.3');
    save(fullfile(pathname,'CompleteArea','processed','data.mat'),'data','-v7.3');
    save(fullfile(pathname,'CompleteArea','processed','mask.mat'),'mask','-v7.3');
    save(fullfile(pathname,'CompleteArea','processed','topo.mat'),'topo','-v7.3');
end