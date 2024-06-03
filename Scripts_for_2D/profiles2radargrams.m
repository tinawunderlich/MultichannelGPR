% Script for creating radargrams.mat for larger data sets from processed
% profiles (Mala MIRA or Spidar data)
% -> choose only some profiles and some channels
%
% Dr. Tina Wunderlich, CAU Kiel 2020-2024, tina.wunderlich@ifg.uni-kiel.de
%
% requires files in profiles2mat and profiles2mat/proc


clear all
close all
clc

% Computer system
platform=2; % Linux=1, Mac=2, Windows=3

% Filename (without folder!) (folder is selected later)
name='NIC_Line'; % Name of datafiles without '_number.mat'
profiles=1:16;  % profile numbers

channels=[1:5]; % choose channels (one or several)

% Choose rSlicer folder (MIRA data) or folder with DT1/HD files (Spidar data)!

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
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder or Spidar-data-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder or Spidar-data-folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Select rSlicer-folder or Spidar-data-folder'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder or Spidar-data-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder or Spidar-data-folder');
        end
    else
        foldername=uigetdir([],'Select rSlicer-folder or Spidar-data-folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
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
elseif platform==3 % Windows (noch nicht getestet!)
    [user,sys]=memory;
    memsize=str2double(sys.PhysicalMemory.Available);
end



% file name parts
[pathstr,fname,ext]=fileparts(name);    % divide name into parts

% load profileinfo.mat
load(fullfile(foldername,'profiles2mat','proc','profileinfo.mat'));
ns=profileinfo(1,3);
dt=profileinfo(1,2);
t=[0:dt:dt*(ns-1)];

% load info files and calculate required memory
disp('Determining required memory...')
anz=1;
global_coords=cell(1,length(profiles)*length(channels));
x=cell(1,length(profiles)*length(channels));
pc=zeros(length(profiles)*length(channels),2);
for i=1:length(profiles)
    if exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(profiles(i)),'_info_proc.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(profiles(i)),'_info_proc.mat']));
        for j=1:length(channels)
            numtraces(anz)=length(info(3,info(3,:)==channels(j)));
            num{anz}=info(9,info(3,:)==channels(j));    % column number in file
            global_coords{anz}=[info(4,info(3,:)==channels(j))' info(5,info(3,:)==channels(j))' info(6,info(3,:)==channels(j))'];
            x{anz}=[0; cumsum(sqrt(diff(global_coords{anz}(:,1)).^2+diff(global_coords{anz}(:,2)).^2))];
            pc(anz,:)=[profiles(i) channels(j)]; % profile and channels number
            anz=anz+1;
        end
    elseif ~exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(profiles(i)),'_info_proc.mat']),'file') && exist(fullfile(foldername,'profiles2mat',[name,'_',int2str(profiles(i)),'_info.mat']),'file')
        load(fullfile(foldername,'profiles2mat',[name,'_',int2str(profiles(i)),'_info.mat']));
        for j=1:length(channels)
            numtraces(anz)=length(info(3,info(3,:)==channels(j)));
            num{anz}=info(9,info(3,:)==channels(j));    % column number in file
            global_coords{anz}=[info(4,info(3,:)==channels(j))' info(5,info(3,:)==channels(j))' info(6,info(3,:)==channels(j))'];
            x{anz}=[0; cumsum(sqrt(diff(global_coords{anz}(:,1)).^2+diff(global_coords{anz}(:,2)).^2))];
            pc(anz,:)=[profiles(i) channels(j)]; % profile and channels number
            anz=anz+1;
        end
    else
        numtraces(anz)=0;
        anz=anz+1*length(channels);
    end
end
% delete empty cells
pc(cellfun('isempty',x),:)=[];
numtraces(numtraces==0)=[];
x=x(~cellfun('isempty',x));
num=num(~cellfun('isempty',num));
global_coords=global_coords(~cellfun('isempty',global_coords));

% calculate approx. size of radargrams.mat
% (each element of the array * bytes in that element (8 for double)) + some overhead (100)
datasize=sum(numtraces(:))*ns*8+100;

if datasize>=memsize/3*2 % if radargrams.mat will be larger than 3/2 of memorysize
    disp('Resulting radargrams.mat will probably be larger than memory size. Please reduce number of profiles or channels and start again.')
    return;
else
    disp('Memory size ok. Start reading of data....')
    radargrams=cell(1,numel(numtraces));
    anz=1;
    for i=1:length(pc(:,1)) % list of profiles and channel numbers
        if exist(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(pc(i,1)),'.mat']),'file')
            m=matfile(fullfile(foldername,'profiles2mat','proc',[name,'_',int2str(pc(i,1)),'.mat']));
            % write all profiles in one variable
            radargrams{anz}=m.traces(:,num{anz});
            anz=anz+1;
        end
        if ~mod(i,10)
            disp(['   ',int2str(i),'/',int2str(length(pc(:,1)))])
        end
    end
    
    % delete empty cells
    radargrams(~cellfun('isempty',radargrams));

    disp('Start saving of data...')
    % save all profiles in one variable
    if length(channels)==1
        if ~exist(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)]),'dir')
            mkdir(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)]));
        end
        save(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)],'radargrams.mat'),'radargrams','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)],'t.mat'),'t','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)],'x.mat'),'x','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)],'global_coords.mat'),'global_coords','-v7.3');

        fid=fopen(fullfile(foldername,'profiles2mat','proc',['Channel',int2str(channels)],'radargrams.txt'),'wt');
        fprintf(fid,'Radargrams.mat contains channels\n');
        fprintf(fid,' %d\t',channels);
        fprintf(fid,'\nof profiles\n');
        fprintf(fid,' %d\n',unique(pc(:,1)));
        fclose(fid);
    else
        save(fullfile(foldername,'profiles2mat','proc','radargrams.mat'),'radargrams','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc','t.mat'),'t','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc','x.mat'),'x','-v7.3');
        save(fullfile(foldername,'profiles2mat','proc','global_coords.mat'),'global_coords','-v7.3');

        fid=fopen(fullfile(foldername,'profiles2mat','proc','radargrams.txt'),'wt');
        fprintf(fid,'Radargrams.mat contains channels\n');
        fprintf(fid,' %d\t',channels);
        fprintf(fid,'\nof profiles\n');
        fprintf(fid,' %d\n',unique(pc(:,1)));
        fclose(fid);
    end
end
disp('Done!')