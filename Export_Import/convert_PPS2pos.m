clear all
close all
clc


%% Script for reading the PPS-related-files outside the rSlicer folder and converting them to *.pos-files (For measurements with GPS mouse/PPS only!)
% (1) Reads the utm, cpdt and cpss-files outside the rSlicer-folder, (2) filters the coordinates for
% removing traces at the same location, (3) determines a time offset by
% comparing PPS signals and applies it to teh coordinates
% and (4) writes *.pos-files with global coordinates (UTM).
% The already available *.pos-files with local coordinates are moved to an
% extra folder named 'local_coords'.
% When asked, select the rSlicer-folder of your data.
%
% Dr. Tina Wunderlich, CAU Kiel, tina.wunderlich@ifg.uni-kiel.de, 2023
%%-------------------------------------------------------------------------


filter_coords=0; % possible values: 0, 1, 2:
% if =0: no filtering
% if =1: coordinates are filtered and traces at
% approx. the same position will be removed (e.g. at the beginning or
% end of the profile): a local difference
% over m coordinates for x and y is calculated separately. If this
% absolute local difference is smaller than nstd*std(diff) it will be removed.
% if =2: coordinates are filtered: local differences over m coordinates are
% calculated for x and y separately and points with differences around zero
% will be deleted (traces at the same place) using nstd as a threshold in
% meters (-/+)
m=9; % number of samples for local difference calculation (must be odd!)
nstd=0.02; % if filter_coords==1: nstd*standard deviation as valid range, if filter_coords==2: +/-nstd [m] will be deleted
show_filter_plot=1; % if =1: show plot for quality control of filtering

coords_smooth=1; % if=1: smooth coordinates along profile (will be applied after filtering)
nsmooth=5; % number of traces for smoothing of coords


%% DO NOT CHANGE FROM HERE ON
if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        foldername=uigetdir([],'Select rSlicer-folder'); % path to radargram-folder

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
            foldername=uigetdir(fn{1}{1},'Select rSlicer-folder');
        else
            foldername=uigetdir([],'Select rSlicer-folder');
        end
    else
        foldername=uigetdir([],'Select rSlicer-folder'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',foldername);
    fclose(fid);
end


% get all available *.pos-files with local coords
list=dir([foldername,'/*.pos']);
weg=cellfun(@(x) startsWith(x,'.'),{list.name}); % these are just temporary files...
list(weg)=[]; % ... delete them from list

if ~isempty(list)
    disp([int2str(length(list)),' *.pos-files found. Moving to folder local_coords.'])
    % make new folder for local files
    if ~exist(fullfile(foldername,'local_coords'),'dir')
        mkdir(fullfile(foldername,'local_coords'));
    end
    for i=1:length(list)
        % get profile number
        test=extractBetween(list(i).name,'_','_');
        number(i,1)=str2num(test{1});

        disp(['   Profile ',int2str(number(i,1))])

        % read utm coords
        temp=extractBefore(list(i).name,'_');
        utm=load(fullfile(foldername,'..',[temp,'_',test{1},'.utm'])); % trace number, N, E, height, UTMZone
        utm(utm(:,1)==-1,:)=[];

        % move original file to new folder
        movefile(fullfile(foldername,list(i).name),fullfile(foldername,'local_coords',list(i).name));

        % determine time offset with PPS and apply it
        utm=PPS_synchro(foldername,list(i).name,utm);

        if isempty(utm)
            disp(['File ',int2str(number(i)),' is too short (test file). Moving into folder "notValid".'])

            if ~exist(fullfile(foldername,'notValid'),'dir')
                mkdir(fullfile(foldername,'notValid'))
            end
            temp=extractBefore(list(i).name,'_G01.pos');
            movefile(fullfile(foldername,'local_coords',list(i).name),fullfile(foldername,'notValid',list(i).name)) % pos-file
            movefile(fullfile(foldername,[temp,'.rad']),fullfile(foldername,'notValid',[temp,'.rad'])) % rad-file
            movefile(fullfile(foldername,[temp,'.rd3']),fullfile(foldername,'notValid',[temp,'.rd3'])) % rd3-file

            continue % skip rest of this loop iteration
        end

        if all(isnan(utm(:,2))) % no correct locations (problems with utm/pps-files)
            disp(['File ',int2str(number(i)),' is corrupt. Moving into folder "notValid".'])

            if ~exist(fullfile(foldername,'notValid'),'dir')
                mkdir(fullfile(foldername,'notValid'))
            end
            temp=extractBefore(list(i).name,'_G01.pos');
            movefile(fullfile(foldername,'local_coords',list(i).name),fullfile(foldername,'notValid',list(i).name)) % pos-file
            movefile(fullfile(foldername,[temp,'.rad']),fullfile(foldername,'notValid',[temp,'.rad'])) % rad-file
            movefile(fullfile(foldername,[temp,'.rd3']),fullfile(foldername,'notValid',[temp,'.rd3'])) % rd3-file

            continue % skip rest of this loop iteration
        end

        if filter_coords>0
            % local differences over m sample
            for j=(m-1)/2+1:length(utm(:,1))-(m-1)/2
                Ndiff(j)=diff(utm([j-(m-1)/2 j+(m-1)/2],2));
                Ediff(j)=diff(utm([j-(m-1)/2 j+(m-1)/2],3));
            end
            Ndiff(j+1:length(utm(:,1)))=0;
            Ediff(j+1:length(utm(:,1)))=0;

            if filter_coords==1
                % mean differences:
                meandiffE=mean(Ediff);
                meandiffN=mean(Ndiff);

                % standard deviation of local diffs
                Nstd=std(Ndiff);
                Estd=std(Ediff);

                % find good points
                inE=(Ediff>=meandiffE-Estd*nstd & Ediff<=meandiffE+Estd*nstd);
                inN=(Ndiff>=meandiffN-Nstd*nstd & Ndiff<=meandiffN+Nstd*nstd);

            elseif filter_coords==2
                % find good points
                inE=~(Ediff>=-nstd & Ediff<=nstd);
                inN=~(Ndiff>=-nstd & Ndiff<=nstd);
            end

            utmfilt=utm(inE & inN,:); % filtered utm list

            if coords_smooth==1 && ~isempty(utmfilt)
                % interpolate coords for all traces
                ind=[utmfilt(1,1):utmfilt(end,1)]';
                x=interp1(utmfilt(:,1),utmfilt(:,2),ind); 
                y=interp1(utmfilt(:,1),utmfilt(:,3),ind);
                z=interp1(utmfilt(:,1),utmfilt(:,4),ind);
                
                % smooth coords
                xs=smooth(x,nsmooth);
                ys=smooth(y,nsmooth);
                zs=smooth(z,nsmooth);

                % extract coords at specific trace numbers again
                utmfiltsmooth=zeros(size(utmfilt));
                for k=1:length(utmfilt(:,1))
                    utmfiltsmooth(k,:)=[utmfilt(k,1) xs(utmfilt(k,1)==ind) ys(utmfilt(k,1)==ind) zs(utmfilt(k,1)==ind)];
                end
            end


            if show_filter_plot==1 && filter_coords==1
                figure
                subplot(2,2,1)
                plot(Ndiff,'*')
                hold on
                yline(meandiffN+Nstd*nstd,'r','LineWidth',2)
                yline(meandiffN-Nstd*nstd,'r','LineWidth',2)
                xlabel('Sample number')
                ylabel('Northing difference [m]')
                set(gca,'Fontsize',20)
                grid on

                subplot(2,2,2)
                plot(Ediff,'*')
                hold on
                yline(meandiffE+Estd*nstd,'r','Linewidth',2)
                yline(meandiffE-Estd*nstd,'r','Linewidth',2)
                xlabel('Sample number')
                ylabel('Easting difference [m]')
                set(gca,'Fontsize',20)
                grid on

                subplot(2,2,3)
                plot(utm(:,3),utm(:,2),'b*','DisplayName','All points')
                hold on
                plot(utm(~inE,3),utm(~inE,2),'kd','DisplayName','filtered out via easting')
                plot(utm(~inN,3),utm(~inN,2),'r+','DisplayName','filtered out via northing')
                set(gca,'FontSize',20)
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                grid on
                legend
                axis equal

                subplot(2,2,4)
                plot(utmfilt(:,3),utmfilt(:,2),'b*','DisplayName','All points filtered')
                hold on
                if exist('utmfiltsmooth','var')
                    plot(utmfiltsmooth(:,3),utmfiltsmooth(:,2),'r+','DisplayName','smoothed')
                end
                set(gca,'FontSize',20)
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                grid on
                legend
                axis equal
            elseif show_filter_plot==1 && filter_coords==2
                figure
                subplot(2,2,1)
                plot(Ndiff,'*')
                hold on
                yline(nstd,'r','LineWidth',2)
                yline(-nstd,'r','LineWidth',2)
                xlabel('Sample number')
                ylabel('Northing difference [m]')
                set(gca,'Fontsize',20)
                grid on

                subplot(2,2,2)
                plot(Ediff,'*')
                hold on
                yline(nstd,'r','Linewidth',2)
                yline(-nstd,'r','Linewidth',2)
                xlabel('Sample number')
                ylabel('Easting difference [m]')
                set(gca,'Fontsize',20)
                grid on

                subplot(2,2,3)
                plot(utm(:,3),utm(:,2),'b*','DisplayName','All points')
                hold on
                plot(utm(~inE,3),utm(~inE,2),'kd','DisplayName','filtered out via easting')
                plot(utm(~inN,3),utm(~inN,2),'r+','DisplayName','filtered out via northing')
                set(gca,'FontSize',20)
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                grid on
                legend
                axis equal

                subplot(2,2,4)
                plot(utmfilt(:,3),utmfilt(:,2),'b*','DisplayName','All points filtered')
                hold on
                if exist('utmfiltsmooth','var')
                    plot(utmfiltsmooth(:,3),utmfiltsmooth(:,2),'r+','DisplayName','smoothed')
                end
                set(gca,'FontSize',20)
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                grid on
                legend
                axis equal
            end

            clear Ndiff;
            clear Ediff;
        else % no filtering
            if coords_smooth==1
                % interpolate coords for all traces
                ind=[utm(1,1):utm(end,1)]';
                x=interp1(utm(:,1),utm(:,2),ind); 
                y=interp1(utm(:,1),utm(:,3),ind);
                z=interp1(utm(:,1),utm(:,4),ind);
                
                % smooth coords
                xs=smooth(x,nsmooth);
                ys=smooth(y,nsmooth);
                zs=smooth(z,nsmooth);

                % extract coords at specific trace numbers again
                utmsmooth=zeros(size(utm));
                for k=1:length(utm(:,1))
                    utmsmooth(k,:)=[utm(k,1) xs(utm(k,1)==ind) ys(utm(k,1)==ind) zs(utm(k,1)==ind)];
                end
            end
        end



        if exist('utmfiltsmooth','var')
            utm=utmfiltsmooth;
        elseif exist('utmsmooth','var')
            utm=utmsmooth;
        else
            utm=utmfilt;
        end


        % write utm coords in pos-file
        fid=fopen(fullfile(foldername,list(i).name),'wt');
        fprintf(fid,'UNITS:m\n');
        fprintf(fid,'%d\t%f\t%f\t%f\n',utm');
        fclose(fid);
    end
else
    disp(['No *.pos-files found.'])
end



%% function for PPS synchronization
function utm=PPS_synchro(foldername,filename,utm)
% split filename:
temp=strsplit(filename,'_');

fid=fopen(fullfile(foldername,'..',[temp{1},'_',temp{2},'.cpss']),'r');

% Header
temp=fread(fid,8,'short');
version=temp(1);
year=temp(2);
month=temp(3);
day=temp(4);
hour=temp(5);
minute=temp(6);
seconds=temp(7);
tie_pos=temp(8);

i=1;
while ~feof(fid) % for every loop 44 bytes forward
    % Trace number
    temp=fread(fid,1,'int'); % 32 bit (4 bytes)
    if isempty(temp)
        break;
    else
        trnum(i,1)=temp;
    end

    % Coordinates (will not be needed here!) -> doesn't matter if wrong
    %     temp=fread(fid,3,'double','b'); % 64 bit, 8 bytes
    %     x(i,1)=temp(1);
    %     y(i,1)=temp(2);
    %     z(i,1)=temp(3);
    fseek(fid,24,'cof'); % jump over coordinates

    % Times
    temp=fread(fid,8,'ushort'); % 16 bit (2 bytes)
    pos_hour(i,1)=temp(1); % time from GPS
    pos_min(i,1)=temp(2);
    pos_sec(i,1)=temp(3);
    pos_ms(i,1)=temp(4);
    pc_hour(i,1)=temp(5); % time from PC
    pc_min(i,1)=temp(6);
    pc_sec(i,1)=temp(7);
    gps_status(i,1)=temp(8); % status GPS

    i=i+1;
end
fclose(fid);

% log: tracenumber PosTime[ms] PCTime[ms] pos_ms[ms]
log_cpss=[trnum pos_hour*60*60*1e3+pos_min*60*1e3+pos_sec*1e3+pos_ms ...
    pc_hour*60*60*1e3+pc_min*60*1e3+pc_sec*1e3 pos_ms];
log_cpss(log_cpss(:,1)<=0,:)=[]; % delete lines with tracenumber <=0
if log_cpss(end-1,2)==log_cpss(end,2)
    log_cpss(end,:)=[];
end
if log_cpss(end-1,1)==log_cpss(end,1)
    log_cpss(end,:)=[]; % delete if the same trace number
end


%% CPDT:
% split filename:
temp=strsplit(filename,'_');

fid=fopen(fullfile(foldername,'..',[temp{1},'_',temp{2},'.cpdt']),'r','l');

% Header:
temp=fread(fid,8,'short'); % short (2 byte)
version=temp(1);
year=temp(2);
month=temp(3);
day=temp(4);
hour=temp(5);
minute=temp(6);
seconds=temp(7);
wheel_num=temp(8);

%ftell(fid) % 16 byte

% Data records:
buffer = fread(fid, inf, 'int8=>int8'); % read in all bytes
fclose(fid);

a=1;
for i=1:32:length(buffer)-32
    cu_time(a,1)=swapbytes(typecast(buffer(i+16:i+17),'int16')); % in 1/1024 s !! time from GPR unit
    wpos(a,1)=swapbytes(typecast(buffer(i+4:i+5),'int16')); % wheel
    trn(a,1)=swapbytes(typecast(buffer(i:i+1),'int16')); % trace number
    pc_time_hr(a,1)=typecast(buffer(i+21:i+24),'int32'); % PC time hours
    pc_time_min(a,1)=typecast(buffer(i+25:i+28),'int32');
    pc_time_sec(a,1)=swapbytes(typecast(buffer(i+29),'int8'));
    a=a+1;
end

PCTime_ms=double(pc_time_hr)*60*60*1e3+double(pc_time_min)*60*1e3+double(pc_time_sec)*1e3;

% log: tracenumber cutime[ms] PCtime[ms] MWheel (hier auch cu_time
% umrechnen von 1/1024s in ms)
log_cpdt=[double(trn) round(double(cu_time)/1024*1e3) PCTime_ms double(wpos)];
log_cpdt(log_cpdt(:,1)<=0,:)=[];
% make new cu_time[ms]
log_cpdt(:,2)=[0; 20; cumsum(diff(log_cpdt(2:end,2)))+20];


%% PPS signal filtering:
% find places where abs(diff(wpos))==1 -> PPS signal comes here
wo=find(abs(diff(log_cpdt(:,4)))==1)+1;

% calculate trace and time differences between PPS signals:
tr_diff=diff(wo);
cu_time_diff=diff(log_cpdt(wo,2));
pps=[wo [0;tr_diff] [0;cu_time_diff]]; % tracenum tracediff timediff

% delete first 5 seconds
pps(log_cpdt(wo,2)<=5e3,:)=[];

% delete pps signals with very small number of trace differences
pps(pps(:,2)<=10,:)=[];

% calculate trace and time differences for remaining PPS signals
pps=[pps(:,1) [0;diff(pps(:,1))] [0;diff(log_cpdt(pps(:,1),2))]];

% if pss was missing inbetween -> tr_diff large -> divide
% if median(pps(:,2))>40 && median(pps(:,2))<60
%     thresh=55; % should be around 50 for 50 Hz (time dist 0.02s) GPR data (=50 traces per second)
% elseif median(pps(:,2))>70 && median(pps(:,2))<130
%     thresh=110; % should be around 100 for 100 Hz (time dist 0.01s) GPR data (=100 traces per second)
% else
%     thresh=median(pps(:,2));
% end
thresh=55;
for i=1:length(pps(:,1))
    temp=pps(i,2);
    a=2;
    while temp>=thresh
        temp=pps(i,2)/a;
        a=a+1;
    end
    factor(i,1)=a-1; % division factor
    pps(i,2)=temp;
end
% also divide time differences by this factor
pps(:,3)=pps(:,3)./factor;

% time diff should be around 1 s -> delete other points
f=pps(:,3)<=900 | pps(:,3)>=1100; % bad points (time diff too small/large) % original: 990/1010
pps(f,:)=[];

% only take the lines where utm and pps have the same trace number
b=1;
weg=zeros(length(log_cpss(:,1)),1);
for i=1:length(log_cpss(:,1))
    findline=find(utm(:,1)==log_cpss(i,1),1,'first'); % line in log_cpss with trace number as in utm
    if isempty(findline) % no match found, because too many locations in utm were filtered out
        weg(i)=1; % this line will be deleted in the log_cpss, because no match found in utm-file (due to filtering)
    else % line was found with match
        utmnew(b,:)=utm(findline,:); % get this line from utm file
        b=b+1;
    end
end
log_cpss(logical(weg),:)=[]; % delete these lines with no match in utm file
utm=utmnew;
%clear utmnew;

wo=pps(:,1); % trace number were PPS signal is ok and utm coordinates are present

%% connect everything:
% go through locations of PPS signals (from cpdt-file) and find next
% tracenumber with 0 ms (postime) from cpss-file:
for i=1:length(wo)
    % find similar trace number for start of search
    strn=find(wo(i)<=log_cpss(:,1),1,'first'); % line number in log_cpss for similar trace number
    % starting from strn-line find next line with pos_ms==0 in cpss-file
    temp1=log_cpss(strn-1+find(log_cpss(strn:end,4)==0,1,'first'),1); % this is the trace number of the corresponding trace
    temp2=log_cpss(strn-1+find(log_cpss(strn:end,4)>=98 & log_cpss(strn:end,4)<=102,1,'first'),1); % if there is no 0 ms in this cycle, take the one around 100 ms
    if ~isempty(temp1) || ~isempty(temp2)
        temp=min([temp1 temp2]); % take the one which comes first (either 0 ms or 101 ms of the next cycle)
        wtrace(i,1)=temp;
        if wtrace(i)>length(log_cpdt(:,1))
            wtrace(i)=length(log_cpdt(:,1));
        end

        tr_offset(i,1)=wtrace(i)-wo(i); % trace offset for this PPS signal
        time_offset(i,1)=log_cpdt(wo(i),2)-log_cpdt(wtrace(i),2); % difference in CU_time
    end
end

if ~exist('tr_offset','var')
    % problem with files, probably only test file...
    utm=[];
else
    % sometimes there is no next 0 ms trace (but only 100 ms after it goes up
    % again)... -> delete these offsets:
    del=tr_offset>=thresh;
    time_offset(del)=[];

    % median offset:
    moffset=median(time_offset); % in ms


    % interpolate coordinates:
    trnum=[1:log_cpss(end,1)]';
    postime_int=round(interp1(log_cpss(:,1),log_cpss(:,2),trnum)); % interpolate pos-time for every  trace
    shift_postime=postime_int-moffset; 
    % sometimes log_cpss(:,2) is not monotonically increasing, but has the
    % same values in between -> make unique
    [~,a,~]=unique(log_cpss(:,2));
    %new_trn=interp1(log_cpss(a,2),log_cpss(a,1),shift_postime);  % new trace numbers
    x=interp1(log_cpss(a,2),utm(a,2),shift_postime,'linear','extrap'); % new x-coords
    y=interp1(log_cpss(a,2),utm(a,3),shift_postime,'linear','extrap'); % new y-coords
    z=interp1(log_cpss(a,2),utm(a,4),shift_postime,'linear','extrap'); % new z-coords

    utm=[utm(:,1) x(utm(:,1)) y(utm(:,1)) z(utm(:,1))];
end

end