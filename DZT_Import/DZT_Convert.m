%------------------------ DZT_Convert -------------------------------------
%
% Converts DZT-files of all GSSI equipments to segy and mat (compatible with Multichannel-GPR)
%
% Dr. Tina Wunderlich, August 2020, tina.wunderlich@ifg.uni-kiel.de
%
%--------------------------------------------------------------------------


clear all
close all
clc

app='SIR4000';    % Equipment: SIR20 / SIR30 / SIR3000 / SIR4000 / Tablet

dataplot=1; % plot radargram for controlling? 1=yes, 0=no

convert2utm=0; % convert WGS Lat/Long to UTM (=1 if measured with Stonex-GPS)
zone=32; % if convert2utm==1 -> give UTM-zone

offsetGPS_X=0; % Offset between GPS and antenna midpoint crossline (in profile direction GPS left of antenna -> positive)
offsetGPS_Y=0; % Offset between GPS and antenna midpoint in profile direction (if GPS behind antenna midpoint -> positive)

% Options for calculating inline coordinates for each trace:
coords_opt=1;   % =1: trace coordinate is difference to beginning of profile (only use this for straight profiles!)
                % =2: trace coordinates are calculated by taking the cumulative sum of the coordinate differences between subsequent traces (better for curvy profiles, but not useful for strong GPS-antenna movements)

% Attention: still experimental! If you think that the output radargrams
% are wrong, set to 0!
removeOutliers=0; % do you want to remove coordinate outliers? 
    % 0= no, use raw coordinates
    % 1= in middle of profile
    % 2= at the end/beginning of profile

% Export to other formats
export2mat=1; % export to Multichannel-GPR format for radargrams (mat-files)
export2segy=0; % export all radargrams as segy-files
constoff=0; % if=1: a constant coordinate offset will be subtracted and coordinates will be in mm accuracy in segy file (offsets will be saved in Inline3D (x) and Crossline3D (y))


%---------------------------- DO NOT CHANGE FROM HERE ON ----------------------------
% 
% LOAD DATA

% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'Export_Import'),fullfile(curFold,'Subfunctions'));


test = exist('OCTAVE_VERSION', 'builtin'); % check if running with matlab or octave
if test==0
    matlab=1;
else
    matlab=0;
end

% get file names
if ispc
    if exist('dzt.temp') % read last opened folder from temp.temp
        fid=fopen('dzt.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZT-file(s)');
        else
            pfad=uigetdir([],'Choose folder with DZT-file(s)');
        end
        fileattrib('dzt.temp','-h');
        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('dzt.temp','+h');
    else
        pfad=uigetdir([],'Choose folder with DZT-file(s)'); % path to radargram-folder

        fid=fopen('dzt.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('dzt.temp','+h');
    end
else
    if exist('.dzt.temp') % read last opened folder from temp.temp
        fid=fopen('.dzt.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with DZT-file(s)');
        else
            pfad=uigetdir([],'Choose folder with DZT-file(s)');
        end
    else
        pfad=uigetdir([],'Choose folder with DZT-file(s)'); % path to radargram-folder
    end

    fid=fopen('.dzt.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end


% get list of files in this folder
list=dir(fullfile(pfad,'*.DZT')); % GSSI format
% check for names starting with .
ii=1;
while ii<length(list)
    if strcmp(list(ii).name(1),'.')
        list(ii,:)=[];
    else
        ii=ii+1;
    end
end

if convert2utm==0 % no conversion
    zone=0;
end

radargrams=[];
global_coords=[];
x=[];
marker=[];
anz=1;
listrad.name=[];
listrad.chan=[];
for i=1:length(list)
    disp(['File ',int2str(i),' of ',int2str(length(list))])
    
    dztname=list(i).name;
        [~,name{i},ext]=fileparts(dztname);
        

        if strcmp(app,'SIR4000')
            [data,trh,h]=readdzt_4000(fullfile(pfad,[name{i},'.DZT']),dataplot,removeOutliers,zone);
        elseif strcmp(app,'Tablet')
            [data,trh,h]=readdzt_Tablet(fullfile(pfad,[name{i},'.DZT']),dataplot);
        elseif strcmp(app,'SIR30')
            [data,trh,h]=readdzt_30(fullfile(pfad,[name{i},'.DZT']),dataplot);
        else
            [data,trh,h]=readdzt_all(fullfile(pfad,[name{i},'.DZT']),dataplot);
        end
        
        if ~strcmp(app,'SIR4000') && convert2utm==1 && exist([pfad,'/', name{i}, '.DZG'])
            lontemp=num2str(trh.x','%.8f');
            lattemp=num2str(trh.y','%.8f');
            for j=1:length(lattemp(:,1))
                temp=strsplit(lattemp(j,:),'.');
                lat(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
                temp=strsplit(lontemp(j,:),'.');
                lon(j)=str2num(temp{1}(1:end-2))+str2num([temp{1}(end-1:end),'.',temp{2}])/60;
            end
            [xneu,yneu]=wgs2utm(lat,lon,zone,'N');
            trh.x=xneu;
            trh.y=yneu;
        end
        
        % correct GPS-antenna offset:
        if offsetGPS_X~=0 || offsetGPS_Y~=0
            anz2=round(0.5/mean(sqrt(diff(trh.x).^2+diff(trh.y).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
            if anz2/2==round(anz2/2)
                anz2=anz2+1; % make odd
            end
            for ii=1:length(trh.x)-anz2
                dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
            end
            anz1=anz2;
            % calculation for the last few traces
            anz2=anz2-1;
            for ii=length(trh.x)-anz1+1:length(trh.x)-1
                dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2);
                temp=helmert([offsetGPS_X offsetGPS_Y],[0 0; 0 dist],[trh.x(ii) trh.y(ii); trh.x(ii+anz2) trh.y(ii+anz2)]);
                trh.x(ii)=temp(1);
                trh.y(ii)=temp(2);
                anz2=anz2-1;
            end
            % extrapolate for last trace
            trh.x(end)=interp1([1 2],[trh.x(end-2) trh.x(end-1)],3,'linear','extrap');
            trh.y(end)=interp1([1 2],[trh.y(end-2) trh.y(end-1)],3,'linear','extrap');
        end
        
        
        % if readdzt_all was used no field time is created
        if ~isfield(trh,'time')
            trh.time = zeros(size(trh.x));
        end
        
        if h.nchan==2
            if i==1
                radargrams=[{data(:,trh.channum==1)}];
                radargrams=[radargrams; {data(:,trh.channum==2)}];
                global_coords=[{[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                global_coords=[global_coords; {[trh.x(trh.channum==2)' trh.y(trh.channum==2)' trh.z(trh.channum==2)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                xtemp=trh.x(trh.channum==2);
                ytemp=trh.y(trh.channum==2);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                time = [{trh.time(trh.channum==1) }];
                time = [time; {trh.time(trh.channum==2) }];
                if ~all(trh.mark==0) % marker present
                    marker=[{trh.mark(trh.channum==1)}];
                    marker=[marker; {trh.mark(trh.channum==2)}];
                end
                listrad.name=[{name{i}}; {name{i}}];
                listrad.chan=[1; 2];
            else
                radargrams=[radargrams; {data(:,trh.channum==1)}];
                radargrams=[radargrams; {data(:,trh.channum==2)}];
                global_coords=[global_coords; {[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                global_coords=[global_coords; {[trh.x(trh.channum==2)' trh.y(trh.channum==2)' trh.z(trh.channum==2)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];     
                end
                xtemp=trh.x(trh.channum==2);
                ytemp=trh.y(trh.channum==2);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];  
                end
                time = [time; {trh.time(trh.channum==1) }];
                time = [time; {trh.time(trh.channum==2) }];
                if ~all(trh.mark==0) % marker present
                    marker=[marker; {trh.mark(trh.channum==1)}];
                    marker=[marker; {trh.mark(trh.channum==2)}];
                end
                listrad.name=[listrad.name; {name{i}}; {name{i}}];
                listrad.chan=[listrad.chan; 1; 2];
            end
            t=double(0:h.dt:h.dt*(h.ns-1));
            if max(t)<1 % for DF-Antenna is dt in s 
                t=t.*1e9; % convert to ns
            end

            if isfield(h,'dt2')
                t2=double(0:h.dt2:h.dt2*(h.ns-1));
                if max(t2)<1
                    t2=t2.*1e9; % convert to ns for DF antenna
                end
            end
            
            chan(anz)=1;
            chan(anz+1)=2;
            anz=anz+2;
        else % 1 channel
            if i==1
                radargrams=[{data(:,trh.channum==1)}];
                global_coords=[{[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[{[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[{[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];
                end
                
                time = [{trh.time(trh.channum==1) }];

                if ~all(trh.mark==0) % marker present
                    marker=[{trh.mark(trh.channum==1)}];
                end
                
                listrad.name=[{name{i}}];
                listrad.chan=[1];
            else
                radargrams=[radargrams; {data(:,trh.channum==1)}];
                global_coords=[global_coords; {[trh.x(trh.channum==1)' trh.y(trh.channum==1)' trh.z(trh.channum==1)']}];
                xtemp=trh.x(trh.channum==1);
                ytemp=trh.y(trh.channum==1);
                if coords_opt==1 % difference to beginning
                    x=[x; {[sqrt((xtemp-xtemp(1)).^2+(ytemp-ytemp(1)).^2)']}];
                elseif coords_opt==2 % cumulative sum
                    x=[x; {[0; cumsum(sqrt(diff(xtemp).^2+diff(ytemp).^2)')]}];    
                end
                
                time = [time; {trh.time(trh.channum==1) }];
                if ~all(trh.mark==0) % marker present
                    marker=[marker; {trh.mark(trh.channum==1)}];
                end
                
                listrad.name=[listrad.name; {name{i}}];
                listrad.chan=[listrad.chan; 1];
            end
            
            t=double(0:h.dt:h.dt*(h.ns-1));
            if max(t)<1 % for DF-Antenna is dt in s 
                t=t.*1e9; % convert to ns
            end
            chan(anz)=1;
            anz=anz+1;
        end
       % collect headers in cell
       headers{i} = h;

end
%------------------------------ WRITE DATA --------------------------------

if export2mat==1
   
    % delete empty cells
    radargrams(~cellfun('isempty',radargrams));
    
    disp('Start saving of mat-files...')
    if matlab==1
        % save all profiles in one variable
        if isfield(h,'dt2') % split into different folders for different time ranges (DF antenna)
            % temporary variables for splitting:
            rtemp=radargrams;
            xtemp=x;
            ctemp=global_coords;
            mtemp=marker;
            ttemp=t;
            timetemp = time;
            clear global_coords;
            clear x;
            clear radargrams;
            clear marker;
            clear time;
            
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch1'),'dir')
                mkdir(fullfile(pfad,'mat_ch1'));
            end
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch2'),'dir')
                mkdir(fullfile(pfad,'mat_ch2'));
            end
            
            % Channel1 
            anz=1;
            for i=1:length(chan)
                if chan(i)==1
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            save(fullfile(pfad,'mat_ch1','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat_ch1','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat_ch1','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat_ch1','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat_ch1','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat_ch1','time.mat'),'time','-v7.3');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch1','marker.mat'),'marker','-v7.3');
            end
            fid=fopen(fullfile(pfad,'mat_ch1','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==1
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % Channel2
            anz=1;
            for i=1:length(chan)
                if chan(i)==2
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            t=t2;
            save(fullfile(pfad,'mat_ch2','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat_ch2','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat_ch2','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat_ch2','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat_ch2','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat_ch2','time.mat'),'time','-v7.3');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch2','marker.mat'),'marker','-v7.3');
            end
            fid=fopen(fullfile(pfad,'mat_ch2','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==2
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % set variables back to original (for sgy export)
            global_coords=ctemp;
            t=ttemp;
            marker=mtemp;
            radargrams=rtemp;
            x=xtemp;
            
        else  % for one channel/not DF antenna
            % make new folder
            if ~exist(fullfile(pfad,'mat'),'dir')
                mkdir(fullfile(pfad,'mat'));
            end
            save(fullfile(pfad,'mat','radargrams.mat'),'radargrams','-v7.3');
            save(fullfile(pfad,'mat','t.mat'),'t','-v7.3');
            save(fullfile(pfad,'mat','x.mat'),'x','-v7.3');
            save(fullfile(pfad,'mat','global_coords.mat'),'global_coords','-v7.3');
            save(fullfile(pfad,'mat','h.mat'),'headers','-v7.3');
            save(fullfile(pfad,'mat','time.mat'),'time','-v7.3');
            if ~isempty(marker)
                save(fullfile(pfad,'mat','marker.mat'),'marker','-v7.3');
            end
            
            fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                fprintf(fid,'%d\t',i);
                fprintf(fid,'%s\t',listrad.name{i});
                fprintf(fid,'%d\n',listrad.chan(i));
            end
            fclose(fid);
        end
    else  % Octave
        % save all profiles in one variable
        if isfield(h,'dt2') % split into different folders for different time ranges (DF antenna)
            % temporary variables for splitting:
            rtemp=radargrams;
            xtemp=x;
            ctemp=global_coords;
            mtemp=marker;
            ttemp=t;
            timetemp = time;
            clear global_coords;
            clear x;
            clear radargrams;
            clear marker;
            clear time;
            
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch1'),'dir')
                mkdir(fullfile(pfad,'mat_ch1'));
            end
            % make new folder
            if ~exist(fullfile(pfad,'mat_ch2'),'dir')
                mkdir(fullfile(pfad,'mat_ch2'));
            end
            
            % Channel1 
            anz=1;
            for i=1:length(chan)
                if chan(i)==1
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            save(fullfile(pfad,'mat_ch1','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat_ch1','t.mat'),'t');
            save(fullfile(pfad,'mat_ch1','x.mat'),'x');
            save(fullfile(pfad,'mat_ch1','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat_ch1','h.mat'),'headers');
            save(fullfile(pfad,'mat_ch1','time.mat'),'time');
            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch1','marker.mat'),'marker');
            end
            fid=fopen(fullfile(pfad,'mat_ch1','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==1
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % Channel2
            anz=1;
            for i=1:length(chan)
                if chan(i)==2
                    radargrams{anz}=rtemp{i};
                    x{anz}=xtemp{i};
                    global_coords{anz}=ctemp{i};
                    time{anz}=timetemp{i};
                    if ~isempty(mtemp)
                        marker{anz}=mtemp{i};
                    end
                    anz=anz+1;
                end
            end
            t=t2;
            save(fullfile(pfad,'mat_ch2','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat_ch2','t.mat'),'t');
            save(fullfile(pfad,'mat_ch2','x.mat'),'x');
            save(fullfile(pfad,'mat_ch2','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat_ch2','h.mat'),'headers');
            save(fullfile(pfad,'mat_ch2','time.mat'),'time');

            if ~isempty(mtemp)
                save(fullfile(pfad,'mat_ch2','marker.mat'),'marker');
            end
            fid=fopen(fullfile(pfad,'mat_ch2','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                if chan(i)==2
                    fprintf(fid,'%d\t',i);
                    fprintf(fid,'%s\t',listrad.name{i});
                    fprintf(fid,'%d\n',listrad.chan(i));
                end
            end
            fclose(fid);
            
            % set variables back to original (for sgy export)
            global_coords=ctemp;
            t=ttemp;
            marker=mtemp;
            radargrams=rtemp;
            x=xtemp;
            
        else  % for one channel/not DF antenna
            % make new folder
            if ~exist(fullfile(pfad,'mat'),'dir')
                mkdir(fullfile(pfad,'mat'));
            end
            save(fullfile(pfad,'mat','radargrams.mat'),'radargrams');
            save(fullfile(pfad,'mat','t.mat'),'t');
            save(fullfile(pfad,'mat','x.mat'),'x');
            save(fullfile(pfad,'mat','global_coords.mat'),'global_coords');
            save(fullfile(pfad,'mat','h.mat'),'headers');
            save(fullfile(pfad,'mat','time.mat'),'time');
            if ~isempty(marker)
                save(fullfile(pfad,'mat','marker.mat'),'marker');
            end
            
            fid=fopen(fullfile(pfad,'mat','radargrams.txt'),'wt');
            fprintf(fid,'Nr.\tName\tChannel\n');
            for i=1:length(listrad.name)
                fprintf(fid,'%d\t',i);
                fprintf(fid,'%s\t',listrad.name{i});
                fprintf(fid,'%d\n',listrad.chan(i));
            end
            fclose(fid);
        end
    end
   
    disp('   Finished!')
end


if export2segy==1
    disp('Start saving as sgy...')
    % make new folder
    if ~exist(fullfile(pfad,'SEGY'),'dir')
        mkdir(fullfile(pfad,'SEGY'));
    end
    
    for i=1:length(radargrams)
        if length(global_coords{i}(1,:))<3
            global_coords{i}(:,3)=zeros(size(global_coords{i}(:,1))); % if no topography present, set to zero
        end
        if chan(i)==2 && isfield(h,'dt2')
            export2sgy2D(radargrams{i},h.dt2,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);
        else
            export2sgy2D(radargrams{i},h.dt,global_coords{i}(:,1),global_coords{i}(:,2),fullfile(pfad,'SEGY',[listrad.name{i},'_Chan',int2str(listrad.chan(i)),'.sgy']),global_coords{i}(:,3),constoff);
        end
    end
    
    disp('   Finished!')
end

% set original path
path(oldpath);