% Script for reading processed Mala-datafiles in profiles2mat/proc and
% binning them onto a rectangular grid
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires MATLAB-files in following folders (path will be temporarily
% set):  Subfunctions


clear all
close all
clc

% Computer system
platform=2; % Linux=1, Mac=2, Windows=3

% Bin size of grid
dx=0.05; % [m]

% Automatic rotation of measurement area for minimum memory size
rotate_area=1;  % 1=yes, 0=no

% Division into rectangles for processing
force_division=0; % if =1: a division is forced although the computer memory is large enough for all data
border=0.5;   % overlapping border of rectangles in m

% Continue reading data for rectangles (if a former run of the program has
% been stopped after some rectangles, you can continue here using the same
% processed data and area configuration as before)
continuedata=0; % if new area: =0
continuerectangle=1;    % number of starting rectangle

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
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            foldername=uigetdir([],'Choose rSlicer folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Choose rSlicer folder');
        else
            foldername=uigetdir([],'Choose rSlicer folder');
        end
    else
        foldername=uigetdir([],'Choose rSlicer folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% get name
temp=dir(fullfile(foldername,'/*.rad'));
tempname=strsplit(temp(end).name,'_'); % Name of data files without '_???.rd3'
name=[tempname{1}];
for i=2:length(tempname)-1
    name=[name,'_',tempname{i}];
end


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Subfunctions'));

% get matlab version (to shut down parfor_wait)
rel=version('-release');
if str2num(rel(1:4))<2017
    pfw=0; % no parfor_wait
else
    pfw=1; % use parfor_wait
end

% get memory information of system
if platform==1 % Linux
    system('less /proc/meminfo >meminfo.txt');
    fid=fopen('meminfo.txt','r');
    temp=textscan(fid,'%s','Headerlines',2);
    fclose(fid);
    memsize=str2double(temp{1}{2})*1e3;    % bytes
elseif platform==2 % Mac
    [status,cmdout]=system('sysctl hw.memsize');
    memsize=str2double(cmdout(13:end)); % memory size in bytes
elseif platform==3 % Windows
    [user,sys]=memory;
    memsize=str2double(sys.PhysicalMemory.Available);
end

%% check if correct processed data available:
if ~exist(fullfile(foldername,'profiles2mat','proc'),'dir')
    disp('No processed data for binning. Stopping!')
    return;
end
% load profileinfo
temp=load(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'));
profileinfo=temp.profileinfo; % profilnumber, dt, ns, channels, numtraces-per-channel
t=0:profileinfo(1,2):profileinfo(1,2)*(profileinfo(1,3)-1);
dt=profileinfo(1,2);
numbers=profileinfo(:,1); % profile numbers
disp(['Found processed data of ',int2str(length(numbers)),' profiles (Min: ',int2str(min(numbers)),'; Max: ',int2str(max(numbers)),').'])


%% Size of area
% read info files for coordinates
disp('Reading coordinates for determining area size...')
xylist=zeros(sum(profileinfo(:,4).*profileinfo(:,5)),6);
anz=0;
for i=1:length(numbers)
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat']),'file')
        m=matfile(fullfile(foldername,'profiles2mat',[name,'_',int2str(numbers(i)),'_info_proc.mat']));
        if i==1
            xylist(1:profileinfo(i,4)*profileinfo(i,5),:)=[zeros(length(m.info(4,:)),1)+numbers(i) m.info(4,:)' m.info(5,:)' m.info(6,:)' m.info(3,:)' [1:profileinfo(i,4)*profileinfo(i,5)]']; % Number, x, y, z, channel of profile, tracenumber in profile
            anz=anz+profileinfo(i,4)*profileinfo(i,5);
        else
            xylist(anz+1:anz+length(m.info(4,:)),:)=[zeros(length(m.info(4,:)),1)+numbers(i) m.info(4,:)' m.info(5,:)' m.info(6,:)' m.info(3,:)' [1:length(m.info(3,:))]']; % Number, x, y, z, channel of profile, tracenumber in profile
            anz=anz+length(m.info(4,:));
        end
    end
end
channels=max(xylist(:,5));
ns=profileinfo(1,3);

%% Configuration file
if continuedata==1  % continue a former run of the program
    if ~exist(fullfile(foldername,'configuration.txt'))
        disp('No former run of program found. Please set continuedata=0 and start script again.')
        return;
    end
    % read configuration file
    fid=fopen(fullfile(foldername,'configuration.txt'),'r');
    temp=textscan(fid,'%f%*s%*s%*s%*s','Headerlines',4);
    num_xrect=temp{1}(1);   % number of rectangles in x direction
    num_yrect=temp{1}(2);   % number of rectangles in y direction
    temp=textscan(fid,'%*s%*s%*s%*s%f',2);
    border=temp{1}(1);
    dx=temp{1}(2);
    temp=textscan(fid,'%*s%*s%*s%f');
    rotbest=temp{1}(1);
    fclose(fid);
    disp('Configuration file found. Using these values:');
    disp(['Dividing area into ',int2str(num_xrect),' rectangles in x direction and ',int2str(num_yrect),' rectangles in y direction.']);
    disp(['Overlapping border of rectangles: ',num2str(border)]);
    disp(['Bin size: ',num2str(dx)]);
    if ~isempty(rotbest)
        disp(['Rotating area by ',num2str(rotbest),' degree.'])
        rmat=[cosd(rotbest) -sind(rotbest); sind(rotbest) cosd(rotbest)]; % rotational matrix
        new=zeros(length(xylist(:,1)),2);
        for rr=1:length(xylist(:,1))
            new(rr,:)=xylist(rr,2:3)*rmat;   % rotate coordinates
        end
        % move origin
        shiftx=floor(min(new(:,1)));
        shifty=floor(min(new(:,2)));
        new(:,1)=new(:,1)-shiftx;
        new(:,2)=new(:,2)-shifty;
        % overwrite coordinates in position
        xylist(:,2:3)=new;

        fig1=figure('Visible','off');
        plot(xylist(:,2),xylist(:,3),'k.')
        hold on
        set(gca,'Dataaspectratio',[1 1 1])
        axis xy
        xlabel('x [m]')
        ylabel('y [m]')
    end
    % load coordtrans of previous rectangle
    load(fullfile(foldername,['3D_Grid_R',int2str(continuerectangle-1)],'coordtrans.mat'));

else  % new run of the program
    if rotate_area==1
        [xylist(:,2:3),rotbest,shiftx,shifty,coordtrans]=rotatearea(xylist(:,2:3));

        fig1=figure('Visible','off');
        plot(xylist(:,2),xylist(:,3),'k.')
        hold on
        set(gca,'Dataaspectratio',[1 1 1])
        axis xy
        xlabel('x [m]')
        ylabel('y [m]')
    end

    % get size of complete gridded area + datatraces matrix + 1000 bytes extra
    ns=profileinfo(1,3);
    channels=profileinfo(1,4);
    areasize=(((max(xylist(:,2))-min(xylist(:,2))+border*2)/dx*(max(xylist(:,3))-min(xylist(:,3))+border*2)/dx*ns)*8+100) + (length(xylist(:,1))*ns*8+100) + 1000; % (each element of the array * bytes in that element (8 for double)) + some overhead (100)
    if memsize/3*2>=areasize
        if force_division==0
            disp('Memory size of computer is larger than required memory! -> Proceed with complete area!')
            num_xrect=1;
            num_yrect=1;
        else
            disp(['Forcing division. Please divide area into rectangles:'])
            num_xrect=str2double(input('How many rectangles in x-direction?  ','s'));
            num_yrect=str2double(input('How many rectangles in y-direction?  ','s'));

            if isempty(num_xrect)
                num_xrect=1;
            end
            if isempty(num_yrect)
                num_yrect=1;
            end
        end
    else
        numreq=2;
        while memsize/3*2<=areasize/numreq
            numreq=numreq+1;
        end
        disp(['Required memory size exceeds computer memory! Please divide area into minimum ',int2str(numreq),' rectangles:'])
        num_xrect=str2double(input('How many rectangles in x-direction?  ','s'));
        num_yrect=str2double(input('How many rectangles in y-direction?  ','s'));
        
        if isempty(num_xrect)
            num_xrect=1;
        end
        if isempty(num_yrect)
            num_yrect=1;
        end  
    end
    
    % Saving configuration file
    disp('Saving configuration file.')
    fid=fopen(fullfile(foldername,'configuration.txt'),'wt');
    fprintf(fid,['Number of channels: ',int2str(profileinfo(1,4)),'\n']);
    fprintf(fid,['Number of samples: ',int2str(ns),'\n']);
    fprintf(fid,['Sampling interval: ',num2str(profileinfo(1,2)),' ns\n']);
    fprintf(fid,['Range: ',num2str((ns-1)*profileinfo(1,2)),' ns\n']);
    fprintf(fid,[int2str(num_xrect),' rectangles in x direction.\n']);
    fprintf(fid,[int2str(num_yrect),' rectangles in y direction.\n']);
    fprintf(fid,['Overlapping border in m: ',num2str(border),'\n']);
    fprintf(fid,['Bin size in m: ',num2str(dx),'\n']);
    if rotate_area==1
        fprintf(fid,['Area rotated by ',int2str(rotbest),' degree.\n']);
        fprintf(fid,['Area shifted by ',num2str(shiftx),' m in x-direction and ',num2str(shifty),' m in y-direction.\n']);
    end
    fclose(fid); 

end
disp('--------------------------------------')


% find min/max of coordinates
minx=floor(min(xylist(:,2)));
maxx=ceil(max(xylist(:,2)));
miny=floor(min(xylist(:,3)));
maxy=ceil(max(xylist(:,3)));


disp('Dividing area into rectangles:')

%%% divide area into rectangles
wid=(maxx-minx)/num_xrect;  % width of rectangles in x-direction
hei=(maxy-miny)/num_yrect;  % height of rectangles in y-direction

disp(['For each rectangle: Area width ',num2str(wid,4),' m, area height ',num2str(hei,4),' m']);


disp('--------------------------------------')
disp('Start reading data in rectangles and binning')

anz=1;
for i=1:num_xrect
    for j=1:num_yrect
        
        tic
        
        % make grid for this rectangle
        xstart=(minx+(i-1)*wid-border)-rem((minx+(i-1)*wid-border)-(minx-border),dx);
        ystart=(miny+(j-1)*hei-border)-rem((miny+(j-1)*hei-border)-(miny-border),dx);
        [x,y]=meshgrid([xstart:dx:minx+i*wid+border],[ystart:dx:miny+j*hei+border]);

        % find traces falling into current rectangle
        ind_in=xylist(:,2)>=minx+(i-1)*wid-border & xylist(:,2)<=minx+i*wid+border & xylist(:,3)>=miny+(j-1)*hei-border & xylist(:,3)<=miny+j*hei+border;

        % interpolate topography onto grid
        z=griddata(xylist(ind_in & xylist(:,5)==channels/2,2),xylist(ind_in & xylist(:,5)==channels/2,3),xylist(ind_in & xylist(:,5)==channels/2,4),x,y);
           
        if any(ind_in==1) % if any data in rectangle
            if (continuedata==1 && continuerectangle<=anz) || continuedata==0
                
                disp(['Reading data in rectangle no. ',int2str(anz),' ...'])
                
                % plot rectangles on coordinate plot
                plot([minx+(i-1)*wid-border minx+i*wid+border minx+i*wid+border minx+(i-1)*wid-border minx+(i-1)*wid-border],[miny+(j-1)*hei-border miny+(j-1)*hei-border miny+j*hei+border miny+j*hei+border miny+(j-1)*hei-border],'Linewidth',2)
                drawnow;
                
                % initialize data grid (3D, 3rd dimension is time) for current rectangle
                data=NaN(length(x(:,1)),length(x(1,:)),ns,'single');
                
                %------------------------------------------------------------------
                %%% Read data files only for current rectangle
                datatraces=zeros(ns,length(ind_in(ind_in==1)),'single');   % initialize matrix for traces
                nn=unique(xylist(ind_in,1)); % profile numbers in rectangle

                if pfw==0
                    for n=1:length(nn) %  loop over profiles in rectangle
                        % load processed data of current profile
                        m=matfile(fullfile(foldername,'profiles2mat',[name,'_',int2str(nn(n)),'_info_proc.mat']));
                        d=matfile(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(nn(n)),'.mat']));
                        
                        % find trace numbers in current profile in rectangle
                        in_rect=xylist(:,6).*ind_in;
                        in_prof=in_rect(xylist(:,1)==nn(n));  % same length as current profile tracenumber
                        
                        % find chunks of neighboring data points
                        chunks=findchunks(in_prof);
                        chunks(:,3)=cumsum(chunks(:,2)-chunks(:,1)+1); % points per chunk
                        
                        % write datatraces and coordinates in matrix
                        data_temp{n}=zeros(ns,length(in_prof(in_prof>0)),'single');
                        x_temp{n}=zeros(1,length(in_prof(in_prof>0)));
                        y_temp{n}=zeros(1,length(in_prof(in_prof>0)));
                        data_temp{n}(:,1:chunks(1,3))=d.traces(:,chunks(1,1):chunks(1,2));
                        x_temp{n}(1:chunks(1,3))=m.info(4,chunks(1,1):chunks(1,2));
                        y_temp{n}(1:chunks(1,3))=m.info(5,chunks(1,1):chunks(1,2));
                        for ii=2:length(chunks(:,1))
                            data_temp{n}(:,chunks(ii-1,3)+1:chunks(ii,3))=d.traces(:,chunks(ii,1):chunks(ii,2));
                            x_temp{n}(chunks(ii-1,3)+1:chunks(ii,3))=m.info(4,chunks(ii,1):chunks(ii,2));
                            y_temp{n}(chunks(ii-1,3)+1:chunks(ii,3))=m.info(5,chunks(ii,1):chunks(ii,2));
                        end
                        % rotate coordinates
                        if rotate_area==1
                            [xy]=apply_rotatearea([x_temp{n}' y_temp{n}'],rotbest,shiftx,shifty);
                            x_temp{n}=xy(:,1)';
                            y_temp{n}=xy(:,2)';
                        end
                    end
                else
                    WaitMessage = parfor_wait(length(nn), 'Waitbar', false,'ReportInterval',1);
                    parfor n=1:length(nn) % parallel loop over profiles in rectangle
                        % load processed data of current profile
                        m=matfile(fullfile(foldername,'profiles2mat',[name,'_',int2str(nn(n)),'_info_proc.mat']));
                        d=matfile(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(nn(n)),'.mat']));
                        
                        % find trace numbers in current profile in rectangle
                        in_rect=xylist(:,6).*ind_in;
                        in_prof=in_rect(xylist(:,1)==nn(n));  % same length as current profile tracenumber
                        
                        % find chunks of neighboring data points
                        chunks=findchunks(in_prof);
                        chunks(:,3)=cumsum(chunks(:,2)-chunks(:,1)+1); % points per chunk
                        
                        % write datatraces and coordinates in matrix
                        data_temp{n}=zeros(ns,length(in_prof(in_prof>0)),'single');
                        x_temp{n}=zeros(1,length(in_prof(in_prof>0)));
                        y_temp{n}=zeros(1,length(in_prof(in_prof>0)));
                        data_temp{n}(:,1:chunks(1,3))=d.traces(:,chunks(1,1):chunks(1,2));
                        x_temp{n}(1:chunks(1,3))=m.info(4,chunks(1,1):chunks(1,2));
                        y_temp{n}(1:chunks(1,3))=m.info(5,chunks(1,1):chunks(1,2));
                        for ii=2:length(chunks(:,1))
                            data_temp{n}(:,chunks(ii-1,3)+1:chunks(ii,3))=d.traces(:,chunks(ii,1):chunks(ii,2));
                            x_temp{n}(chunks(ii-1,3)+1:chunks(ii,3))=m.info(4,chunks(ii,1):chunks(ii,2));
                            y_temp{n}(chunks(ii-1,3)+1:chunks(ii,3))=m.info(5,chunks(ii,1):chunks(ii,2));
                        end
                        % rotate coordinates
                        if rotate_area==1
                            [xy]=apply_rotatearea([x_temp{n}' y_temp{n}'],rotbest,shiftx,shifty);
                            x_temp{n}=xy(:,1)';
                            y_temp{n}=xy(:,2)';
                        end
                        WaitMessage.Send;
                    end
                    WaitMessage.Destroy;
                end

                % put together all data
                datatraces=zeros(ns,length(ind_in(ind_in==1)),'single');
                xx=zeros(1,length(ind_in(ind_in==1)));
                yy=zeros(1,length(ind_in(ind_in==1)));
                ii=1;
                for n=1:length(nn)
                    datatraces(:,ii:ii+length(data_temp{n}(1,:))-1)=data_temp{n};  % -> all traces
                    xx(ii:ii+length(data_temp{n}(1,:))-1)=x_temp{n};
                    yy(ii:ii+length(data_temp{n}(1,:))-1)=y_temp{n};
                    ii=ii+length(data_temp{n}(1,:));
                end
                clear data_temp x_temp y_temp;
                
                disp('Reading of data completed!')
                disp(' ')
                            
                %%%
                % make new folder for each rectangle
                if ~exist(fullfile(foldername,['3D_Grid_R',int2str(anz)]),'dir')
                    mkdir(fullfile(foldername,['3D_Grid_R',int2str(anz)]));
                end
               
                %------------------------------------------------------------------
                disp('Binning of processed data...')
                %%% Binning of processed data
                xrg=x(1,1)-dx/2:dx:x(1,end)+dx/2;
                yrg=y(1,1)-dx/2:dx:y(end,1)+dx/2;

                data=bindata3(datatraces,xx,yy,xrg,yrg); % fast binning

                            
                %%% make mask (1 where data and no overlapping)
                mask=zeros(size(x));
                mask(~isnan(data(:,:,round(length(data(1,1,:))/2))))=1; % set =1 in bins with data
                mask(x>x(1,end)-border & x<=x(1,1)+border & y>y(end,1)-border & y<=y(1,1)+border)=0;    % set =0 in overlapping region
                
                % delete gridding artefacts of z
                z=z.*mask;
                
                disp(' ')
                
                %%% Saving data
                disp('Saving data of current rectangle...')
                
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'data.mat'),'data','-v7.3');
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'x.mat'),'x','-v7.3');
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'y.mat'),'y','-v7.3');
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'z.mat'),'z','-v7.3');
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'t.mat'),'t','-v7.3');
                save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'mask.mat'),'mask','-v7.3');
                if exist('coordtrans','var')
                    save(fullfile(foldername,['3D_Grid_R',int2str(anz)],'coordtrans.mat'),'coordtrans','-v7.3');
                end
               
                % save figure with rectangle
                hf=figure('Name','CurrentRectangle','visible','off');
                hold off
                plot(xylist(:,2),xylist(:,3),'k.')
                hold on
                plot([minx+(i-1)*wid-border minx+i*wid+border minx+i*wid+border minx+(i-1)*wid-border minx+(i-1)*wid-border],[miny+(j-1)*hei-border miny+(j-1)*hei-border miny+j*hei+border miny+j*hei+border miny+(j-1)*hei-border],'Linewidth',2)
                axis xy
                set(gca,'Dataaspectratio',[1 1 1])
                xlabel('x [m]')
                ylabel('y [m]')
                print(fullfile(foldername,['3D_Grid_R',int2str(anz)],'Rectangle_Location.jpg'),'-djpeg');
                close(hf);
                             
                toc
                disp('--------------------------------------')
                
                
                %%% deleting variables for next rectangle
                clear datatraces data x y z mask ind_in;
                
            end
            
            % increase number of rectangle
            anz=anz+1;
        end
    end
end

if num_xrect>1 || num_yrect>1
    fid=fopen(fullfile(foldername,'configuration.txt'),'a+');
    fprintf(fid,['Number of rectangles with data: ',int2str(anz-1),'\n']);
    fclose(fid);
end

% set original path
path(oldpath);

% End of script.

%--------------------------------------------------------------------------

function [xy,rotbest,shiftx,shifty,coordtrans]=rotatearea(xy)
%%% Rotate area for minimum memory
disp('Find optimum rotation angle...')
rot=[-45:5:45];
new=zeros(size(xy));
for r=1:length(rot)
    rmat=[cosd(rot(r)) -sind(rot(r)); sind(rot(r)) cosd(rot(r))]; % rotational matrix
    for rr=1:length(xy(:,1))
        new(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
    end
    area(r)=(max(new(:,1))-min(new(:,1)))*(max(new(:,2))-min(new(:,2)));    % area size of rotated coordinates
end
rotbest=rot(area==min(area)); % best rotation angle => smallest area
disp(['Optimum rotation angle is ',num2str(rotbest),' degree. Area has been rotated. Saving coordtrans.mat for later transformation.'])
rmat=[cosd(rotbest) -sind(rotbest); sind(rotbest) cosd(rotbest)]; % rotational matrix
for rr=1:length(xy(:,1))
    new(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
end
% move origin
shiftx=floor(min(new(:,1)));
shifty=floor(min(new(:,2)));
new(:,1)=new(:,1)-shiftx;
new(:,2)=new(:,2)-shifty;
% save coordinate pairs for later transformation
coordtrans=[new(new(:,1)==min(new(:,1)),:) xy(new(:,1)==min(new(:,1)),:);...
    new(new(:,1)==max(new(:,1)),:) xy(new(:,1)==max(new(:,1)),:);...
    new(new(:,2)==min(new(:,2)),:) xy(new(:,2)==min(new(:,2)),:);...
    new(new(:,2)==max(new(:,2)),:) xy(new(:,2)==max(new(:,2)),:)]; % [local x, local y, global x, global y]
% overwrite coordinates in position
xy=new;
disp(['Area size is now ',int2str(round(range(new(:,1)))),' x ',int2str(round(range(new(:,2)))),' m (x/y).'])
end

function [xy]=apply_rotatearea(xy,rot,shiftx,shifty)
%%% Rotate area with given parameters
rmat=[cosd(rot) -sind(rot); sind(rot) cosd(rot)]; % rotational matrix
for rr=1:length(xy(:,1))
    xy(rr,:)=xy(rr,:)*rmat;   % rotate coordinates
end
% move origin
xy(:,1)=xy(:,1)-shiftx;
xy(:,2)=xy(:,2)-shifty;
end

function chunks=findchunks(in_prof)
% find chunks in in_prof (=blocks of following traces that have to be read
% in)
chunks=[];
flag=0;
for ii=1:length(in_prof)-1
    if (in_prof(ii)==0 && in_prof(ii+1)>0)
        flag=flag+1;
        chunks(flag,1)=ii+1;  % start of interval with data
        if ii+1==length(in_prof) % if only last value ~=0
            chunks(flag,2)=ii+1; % set also this value as ending point
        end
    elseif (ii==1 && in_prof(ii)>0) % beginning of line
        flag=flag+1;
        chunks(flag,1)=ii;  % start of interval with data
    elseif (in_prof(ii)>0 && in_prof(ii+1)==0)
        flag=flag+1;
        chunks(flag-1,2)=ii;    % end of interval with data
    elseif (ii+1==length(in_prof) && in_prof(ii+1)>0) % end of line
        flag=flag+1;
        chunks(flag-1,2)=ii+1;    % end of interval with data
    end
end
chunks(chunks(:,1)==0 & chunks(:,2)==0,:)=[];
end

