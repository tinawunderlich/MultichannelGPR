clear all
close all
clc

%% Join two datasets with radargrams in mat-format into one dataset
%
% Dr. Tina Wunderlich, CAU Kiel, 2022, tina.wunderlich@ifg.uni-kiel.de
%
% The user has to choose two folders with radargrams in mat-format. The
% joined data is saved in a new folder 'mergedRadargrams' in location of
% first dataset.
% Until now only tested for time domain and with same time vectors...

timedepth=1; % time=1, depth=2



%% ---- DO NOT CHANGE BELOW THIS LINE -------------------------------------

% get folder name
if ispc
    if exist('temp1.temp') % read last opened folder from temp.temp
        fid=fopen('temp1.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            foldername1=uigetdir(fn{1}{1},'Choose folder with dataset 1');
        else
            foldername1=uigetdir([],'Choose folder with dataset 1');
        end
        fid=fopen('temp1.temp','wt');
        fprintf(fid,'%s',foldername1);
        fclose(fid);
    else
        foldername1=uigetdir([],'Choose folder with dataset 1'); % path to radargram-folder
        
        fid=fopen('temp1.temp','wt');
        fprintf(fid,'%s',foldername1);
        fclose(fid);
    end
else
    if exist('.temp1.temp') % read last opened folder from temp.temp
        fid=fopen('.temp1.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername1=uigetdir(fn{1}{1},'Choose folder with dataset 1');
        else
            foldername1=uigetdir([],'Choose folder with dataset 1');
        end
    else
        foldername1=uigetdir([],'Choose folder with dataset 1'); % path to radargram-folder
    end
    
    fid=fopen('.temp1.temp','wt');
    fprintf(fid,'%s',foldername1);
    fclose(fid);
end

if ispc
    if exist('temp2.temp') % read last opened folder from temp.temp
        fid=fopen('temp2.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            foldername2=uigetdir(fn{1}{1},'Choose folder with dataset 2');
        else
            foldername2=uigetdir([],'Choose folder with dataset 2');
        end
        fid=fopen('temp2.temp','wt');
        fprintf(fid,'%s',foldername2);
        fclose(fid);
    else
        foldername2=uigetdir([],'Choose folder with dataset 2'); % path to radargram-folder
        
        fid=fopen('temp2.temp','wt');
        fprintf(fid,'%s',foldername2);
        fclose(fid);
    end
else
    if exist('.temp2.temp') % read last opened folder from temp.temp
        fid=fopen('.temp2.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername2=uigetdir(fn{1}{1},'Choose folder with dataset 2');
        else
            foldername2=uigetdir([],'Choose folder with dataset 2');
        end
    else
        foldername2=uigetdir([],'Choose folder with dataset 1'); % path to radargram-folder
    end
    
    fid=fopen('.temp2.temp','wt');
    fprintf(fid,'%s',foldername2);
    fclose(fid);
end


% dataset 1
r1=load(fullfile(foldername1,'radargrams.mat'));
gc1=load(fullfile(foldername1,'global_coords.mat'));
t1=load(fullfile(foldername1,'t.mat'));
x1=load(fullfile(foldername1,'x.mat'));

dt1=t1.t(2)-t1.t(1);
ns1=length(t1.t);

% dataset 2
r2=load(fullfile(foldername2,'radargrams.mat'));
gc2=load(fullfile(foldername2,'global_coords.mat'));
t2=load(fullfile(foldername2,'t.mat'));
x2=load(fullfile(foldername2,'x.mat'));

dt2=t2.t(2)-t2.t(1);
ns2=length(t2.t);

if timedepth==1
    if dt1==dt2 && ns1==ns2 % same time vector
        if timedepth==1
            disp('Datasets have the same time vector. Merging datasets.')
        else
            disp('Datasets have the same depth vector. Merging datasets.')
        end
        radargrams=r1.radargrams;
        global_coords=gc1.global_coords;
        t=t1.t;
        x=x1.x;
        % add second dataset:
        anz=length(radargrams);
        for i=1:length(r2.radargrams)
            radargrams{anz+i}=r2.radargrams{i};
            global_coords{anz+i}=gc2.global_coords{i};
            x{anz+i}=x2.x{i};
        end
        dt=dt1;
        ns=ns1;
    elseif dt1==dt2 && ns1~=ns2 % same dt, but different sample number
        disp('Datasets have the same dt, but different number of samples. Continue with higehr number of samples and pad other dataset.')
        for i=1:length(r1.radargrams)
            radargrams{i}=zeros(max([ns1,ns2]),length(r1.radargrams{i}(1,:)));
            radargrams{i}(1:ns1,:)=r1.radargrams{i};
            global_coords{i}=gc1.global_coords{i};
            x{i}=x1.x{i};
        end
        anz=length(radargrams);
        for i=1:length(r2.radargrams)
            radargrams{anz+i}=zeros(max([ns1,ns2]),length(r2.radargrams{i}(1,:)));
            radargrams{anz+i}(1:ns2,:)=r2.radargrams{i};
            global_coords{anz+i}=gc2.global_coords{i};
            x{anz+i}=x2.x{i};
        end
        t=0:dt1:(max([ns1,ns2])-1)*dt1;
        dt=dt1;
        ns=length(t);
    elseif dt1~=dt2 % different dt -> interpolate new
        disp(['Datasets have different sample intervals: dt1=',num2str(dt1,2),'ns and dt2=',num2str(dt2,2),'ns'])
        dt=input('Which sample interval do you prefer? dt[ns]= ...'); % ns
        t=0:dt:max([t1.t,t2.t]);
        
        for i=1:length(r1.radargrams)
            radargrams{i}=zeros(length(t),length(r1.radargrams{i}(1,:)));
            for j=1:length(r1.radargrams{i}(1,:))
                radargrams{i}(:,j)=interp1(t1.t,r1.radargrams{i}(:,j),t);
            end
            global_coords{i}=gc1.global_coords{i};
            x{i}=x1.x{i};
        end
        anz=length(radargrams);
        for i=1:length(r2.radargrams)
            radargrams{anz+i}=zeros(length(t),length(r2.radargrams{i}(1,:)));
            for j=1:length(r2.radargrams{i}(1,:))
                radargrams{anz+i}(:,j)=interp1(t2.t,r2.radargrams{i}(:,j),t);
            end
            global_coords{anz+i}=gc2.global_coords{i};
            x{anz+i}=x2.x{i};
        end
        ns=length(t);
    end
else  % depth
    disp(['Datasets have following depth sample intervals: dz1=',num2str(dt1,2),'m and dz2=',num2str(dt2,2),'m'])
    dt=input('Which sample interval do you prefer? dz[m]= ...'); % m
    
    % depth is in t-vector
    t=max([max(t1.t),max(t2.t)]):dt:min([min(t1.t),min(t2.t)]);
    for i=1:length(r1.radargrams)
        radargrams{i}=zeros(length(t),length(r1.radargrams{i}(1,:)));
        for j=1:length(r1.radargrams{i}(1,:))
            radargrams{i}(:,j)=interp1(t1.t,r1.radargrams{i}(:,j),t);
        end
        global_coords{i}=gc1.global_coords{i};
        x{i}=x1.x{i};
    end
    anz=length(radargrams);
    for i=1:length(r2.radargrams)
        radargrams{anz+i}=zeros(length(t),length(r2.radargrams{i}(1,:)));
        for j=1:length(r2.radargrams{i}(1,:))
            radargrams{anz+i}(:,j)=interp1(t2.t,r2.radargrams{i}(:,j),t);
        end
        global_coords{anz+i}=gc2.global_coords{i};
        x{anz+i}=x2.x{i};
    end
    ns=length(t);
end

% save joined data:
if ~exist(fullfile(foldername1,'mergedRadargrams'))
    mkdir(fullfile(foldername1,'mergedRadargrams'));
end
save(fullfile(foldername1,'mergedRadargrams','t.mat'),'t','-v7.3');
save(fullfile(foldername1,'mergedRadargrams','radargrams.mat'),'radargrams','-v7.3');
save(fullfile(foldername1,'mergedRadargrams','x.mat'),'x','-v7.3');
save(fullfile(foldername1,'mergedRadargrams','global_coords.mat'),'global_coords','-v7.3');
% info file:
fid=fopen(fullfile(foldername1,'mergedRadargrams','merge_info.txt'),'wt');
fprintf(fid,'Following datasets have been merged:\n');
fprintf(fid,'%s\n',foldername1);
fprintf(fid,'%s\n',foldername2);
fprintf(fid,'dt of new dataset is %4.2f ns.\nns of new dataset is %d.',[dt ns]);
fclose(fid);

disp(['Datasets have been saved in new folder "mergedRadargrams" in ',foldername1]);