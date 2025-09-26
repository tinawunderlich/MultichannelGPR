function []=processingTestGUI(folder,name,profilelist,ch_list,equipmentflag,utmzone)

% GUI for plotting of profiles and interactive testing of processing steps
% Mala MIRA and Spidar data and Impulse Radar data
%
% Dr. Tina Wunderlich, CAU Kiel 2020-2025, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% folder: path and name of rSlicer-folder
% name: name of data
% profilelist: number of profiles

if equipmentflag==1
    S.equipmentflag=1; % Mala Mira
elseif equipmentflag==2
    S.equipmentflag=2; % Spidar
    S.ch_list=ch_list;

    % get gps-channel:
    temp=dir(fullfile(folder,'/*.GPS'));
    for i=1:length(temp)
        if ~startsWith(temp(i).name,'.')
            temp1=strsplit(temp(i).name,'_'); % split in two parts (e.g. NIC01 and Line001.GPS)
            S.gps_channel=str2double(temp1{1}(end-1:end)); % channelNumber for gps data
            break;
        end
    end
elseif equipmentflag==3 % Impulse Radar
    S.equipmentflag=3;
    S.ch_list=ch_list;
end

S.folder=folder;
S.name=name;
S.profilelist=profilelist;
S.utmzone=utmzone;

% load first profile
if S.equipmentflag==1
    % Mala
    [S.traces,S.dt,S.ns,S.x,S.y,S.z]=readmala(S.folder,S.name,S.profilelist(1));
    S.ch_list=1:length(S.traces);
elseif S.equipmentflag==2
    %Spidar
    [temptraces,S.dt,S.ns,tempx,tempy,tempz]=readspidar(S.folder,S.name,S.profilelist(1),ch_list,zeros(length(ch_list),1),zeros(length(ch_list),1),S.gps_channel,S.utmzone,0);
    trPerChan=length(temptraces(1,:))/length(ch_list);
    for i=1:length(ch_list)
        S.traces{i}=temptraces(:,(i-1)*trPerChan+1:i*trPerChan);
        S.x{i}=tempx(:,(i-1)*trPerChan+1:i*trPerChan)';
        S.y{i}=tempy(:,(i-1)*trPerChan+1:i*trPerChan)';
        S.z{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
    end
elseif S.equipmentflag==3
    % Impulse Radar
    [temptraces,S.dt,S.ns,tempx,tempy,tempz,numchannels]=readImpulseRadar(S.folder,S.name,S.profilelist(1),0,0,S.utmzone,1,5);
    trPerChan=size(temptraces,2)/numchannels;
    for i=1:numchannels
        S.traces{i}=temptraces(:,(i-1)*trPerChan+1:i*trPerChan);
        S.x{i}=tempx(:,(i-1)*trPerChan+1:i*trPerChan)';
        S.y{i}=tempy(:,(i-1)*trPerChan+1:i*trPerChan)';
        S.z{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
    end
end

S.channellist=[{'All'}];
S.raw=[];
S.info=[]; % channel number corresponding to traces in raw/data
S.zz=[];
for i=1:length(S.traces)
    S.xprof{i}=[0; cumsum(sqrt(diff(S.x{i}).^2+diff(S.y{i}).^2))]; % local coords for this channel in m
    weg=find(isnan(S.xprof{i})); % find traces with NaN coordinates (e.g. if standing at the end of the profile)
    S.traces{i}(:,weg)=[];
    S.x{i}(weg)=[];
    S.y{i}(weg)=[];
    S.z{i}(weg)=[];
    S.zz=[S.zz; S.z{i}];
    S.xprof{i}(weg)=[];
    S.raw=[S.raw S.traces{i}];
    if S.equipmentflag==1 || S.equipmentflag==3 % Mala or impulseRadar
        S.channellist=[S.channellist; {int2str(i)}];
        S.info=[S.info zeros(1,length(S.xprof{i}))+i]; % channelnumber
    else % spidar
        S.channellist=[S.channellist; {int2str(ch_list(i))}];
        S.info=[S.info zeros(1,length(S.xprof{i}))+ch_list(i)]; % channelnumber
    end
end
S.info_proc=S.info; % channelnumber of processed data

dx=mean(diff(S.xprof{1})); % mean trace spacing
S.xloc=0:dx:(length(S.info)-1)*dx; % local coords for all channels behind each other

S.t=0:S.dt:S.dt*(S.ns-1); % time vector


% traces are original raw data in cells, raw ist raw data of all channels, proc
% is processed data of all channels, choose channels via info (channel
% number for each trace), xloc is x-axis for all channels, xprof is x-axsi
% for each channel in cells

% get screensize
S.siz = get( 0, 'Screensize' );
S.fh.Position=S.siz;  % default: screensize

% Add the UI components
S=addcomponents(S);

% Make figure visible after adding components
S.fh.Visible='on';

    function S=addcomponents(S)
        
        % create figure handle:
        S.fh = figure('units','pixels',...
            'position',S.siz,...
            'menubar','none',...
            'name','Processing test',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui);
        
        % make initial plot of radargram
        S.ax=axes('unit','pix',...
            'position',[250 300 S.siz(3)-300 S.siz(4)-320]);
        S.rad=imagesc(S.xloc,S.t,S.raw);
        hold on
        for i=1:length(S.ch_list)-1
            plot([max(S.xloc(S.info==S.ch_list(i))) max(S.xloc(S.info==S.ch_list(i)))],[0 max(S.t)],'r')
        end
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('t [ns]')
        set(S.ax,'DataAspectratio',[1 1 1])
        grid on
        axis tight
        
        %% UI controls:
        % UI control for channel number
        S.chan=uicontrol(S.fh,'style','listbox','unit','pix','position',[20 310 55 100],'String',S.channellist,'Value',1,'callback',@channum_call);
        S.chantext=uicontrol(S.fh,'style','text','unit','pix','position',[20 420 55 15],'String','Channel','HorizontalAlignment','left');
        
        % UI control for profile number
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[80 310 55 100],'String',S.profilelist,'Value',1,'callback',@profnum_call);
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[80 420 55 15],'String','Profile','HorizontalAlignment','left');
        
        % UI control aspectratio
        S.asp_String=[{'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'}];
        S.aspval=1/2;
        S.asp = uicontrol(S.fh,'style','listbox','unit','pix','position',[140 310 55 100],'String',S.asp_String,'Value',7,'Callback',@asp_call);
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[140 420 85 15],'String','Aspect ratio','HorizontalAlignment','left');
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol(S.fh,'Style','radiobutton','String','Auto','Position',[20 590-80 100 15],'Value',1,'Callback',@colorscale_call); % this button is on
        S.rb2=uicontrol(S.fh,'Style','radiobutton','String','1 %','Position',[20 570-80 100 15],'Value',0,'Callback',@colorscale_call);
        S.rb3=uicontrol(S.fh,'Style','radiobutton','String','3 %','Position',[20 550-80 100 15],'Value',0,'Callback',@colorscale_call);
        S.coltext=uicontrol(S.fh,'style','text','unit','pix','position',[20 610-80 100 15],'String','Color scale','HorizontalAlignment','left');
        S.flag_cs=1;
        
        % UI control - radiobuttons for data (raw/proc)
        S.d1=uicontrol(S.fh,'Style','radiobutton','String','raw','Position',[120 590-80 100 15],'Value',1,'Callback',@data_call); % this button is on
        S.d2=uicontrol(S.fh,'Style','radiobutton','String','processed','Position',[120 570-80 100 15],'Value',0,'Enable','off','Callback',@data_call);
        S.datatext=uicontrol(S.fh,'style','text','unit','pix','position',[120 610-80 100 15],'String','Data','HorizontalAlignment','left');
        S.flag_data=1;
        
        % axes for waitbar:
        S.wbax=axes('unit','pix','position',[20 270 180 10]);
        S.wb=patch([0 100 100 0 0],[0 0 1 1 0],[1 1 1]);
        set(S.wbax,'XTick',[]','XTicklabel',[],'YTick',[],'YTickLabel',[])
        
        % UI control for processing steps list
        S.proclist=reorderableListbox(S.fh,'Style','listbox','unit','pix','position',[20 110 180 120],'String',[]);
        S.proclisttext=uicontrol(S.fh,'style','text','unit','pix','position',[20 240 180 15],'String','Processing steps','HorizontalAlignment','left');
        
        % UI control for Apply/delete/save
        S.apply=uicontrol(S.fh,'Style','pushbutton','unit','pix','position',[20 80 180 20],'Enable','off','String','Apply','Callback',@apply_call);
        S.save=uicontrol(S.fh,'Style','pushbutton','unit','pix','position',[20 50 180 20],'Enable','off','String','Save settings','Callback',@save_call);
        S.delete=uicontrol(S.fh,'Style','pushbutton','unit','pix','position',[20 20 180 20],'Enable','off','String','Delete list','Callback',@del_call);
        
        
        %%% Processing steps
        % 1st column
        off1=70;
        % DCremoval
        S.DCrem=uicontrol(S.fh,'Style','checkbox','String','Remove amplitude offset','FontWeight','bold','Position',[300-off1 300 200 15],'Value',0,'Callback',@DCrem_call);
        S.tstart=uicontrol(S.fh,'Style','edit','String','0','Value',0,'Enable','off','Position',[320-off1 270 35 20]);
        S.tend=uicontrol(S.fh,'Style','edit','String','100','Value',100,'Enable','off','Position',[320-off1 240 35 20]);
        S.tstarttext=uicontrol(S.fh,'Style','text','String','tstart [ns]','Position',[360-off1 270 140 20],'HorizontalAlignment','left');
        S.tendtext=uicontrol(S.fh,'Style','text','String','tend [ns]','Position',[360-off1 240 140 20],'HorizontalAlignment','left');
        
        % Trace interpolation
        S.badtrace=uicontrol(S.fh,'Style','checkbox','String','Bad trace removal','FontWeight','bold','Position',[300-off1 210 200 15],'Value',0,'Callback',@badtraceremoval_call);
        S.minfactor=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[320-off1 180 35 20]);
        S.maxfactor=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[320-off1 150 35 20]);
        S.minfactortext=uicontrol(S.fh,'Style','text','String','minfactor','Position',[360-off1 180 140 20],'HorizontalAlignment','left');
        S.maxfactortext=uicontrol(S.fh,'Style','text','String','maxfactor','Position',[360-off1 150 140 20],'HorizontalAlignment','left');
        
        % Spherical divergence
        S.spherDiv=uicontrol(S.fh,'Style','checkbox','String','Spherical divergence','FontWeight','bold','Position',[300-off1 120 200 15],'Value',0,'Callback',@spherdiv_call);
        
        % attenuation correction
        S.attenuation=uicontrol(S.fh,'Style','checkbox','String','Attenuation correction','FontWeight','bold','Position',[300-off1 90 200 15],'Value',0,'Callback',@attenuation_call);
        S.eps=uicontrol(S.fh,'Style','edit','String','9','Enable','off','Position',[320-off1 60 35 20]);
        S.sigma=uicontrol(S.fh,'Style','edit','String','0.02','Enable','off','Position',[320-off1 30 35 20]);
        S.epstext=uicontrol(S.fh,'Style','text','String','rel. permittivity','Position',[360-off1 60 140 20],'HorizontalAlignment','left');
        S.sigmatext=uicontrol(S.fh,'Style','text','String','conductivity [S/m]','Position',[360-off1 30 140 20],'HorizontalAlignment','left');
        
        % 2nd column
        off2=30;
        % t0 correction
        S.t0corr=uicontrol(S.fh,'Style','checkbox','String','t0 correction','FontWeight','bold','Position',[500-off1-off2 300 200 15],'Value',0,'Callback',@t0_call);
        S.t0list=uicontrol(S.fh,'Style','listbox','unit','pix','Enable','off','position',[520-off1-off2 240 170 55],'Callback',@t0list_call,'String',[{'Channelshift for all profiles'}; {'Channelshift individually'}; {'Amplitude threshold'}; {'t0 correction with reference trace'}]);
        S.thresh=uicontrol(S.fh,'Style','edit','String','-1000','Enable','off','Position',[520-off1-off2 210 35 20]);
        S.refprof=uicontrol(S.fh,'Style','edit','String',num2str(S.profilelist(1)),'Enable','off','Position',[520-off1-off2 180 35 20]);
        S.refchan=uicontrol(S.fh,'Style','edit','String','1','Enable','off','Position',[520-off1-off2 150 35 20]);
        S.reftrace=uicontrol(S.fh,'Style','edit','String','1','Enable','off','Position',[520-off1-off2 120 35 20]);
        S.t0=uicontrol(S.fh,'Style','edit','String','0.5','Enable','off','Position',[520-off1-off2 90 35 20]);
        S.maxshift=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[520-off1-off2 60 35 20]);
        S.shiftsamples=uicontrol(S.fh,'Style','edit','String','','Enable','off','Position',[520-off1-off2 30 70 20]);
        S.threshtext=uicontrol(S.fh,'Style','text','String','threshold','Position',[560-off1-off2 210 140 20],'HorizontalAlignment','left');
        S.refproftext=uicontrol(S.fh,'Style','text','String','reference profile','Position',[560-off1-off2 180 140 20],'HorizontalAlignment','left');
        S.refchantext=uicontrol(S.fh,'Style','text','String','reference channel','Position',[560-off1-off2 150 140 20],'HorizontalAlignment','left');
        S.reftracetext=uicontrol(S.fh,'Style','text','String','reference trace','Position',[560-off1-off2 120 140 20],'HorizontalAlignment','left');
        S.t0text=uicontrol(S.fh,'Style','text','String','t0 [ns]','Position',[560-off1-off2 90 140 20],'HorizontalAlignment','left');
        S.maxshifttext=uicontrol(S.fh,'Style','text','String','maxshift','Position',[560-off1-off2 60 140 20],'HorizontalAlignment','left');
        S.shitfsamplestext=uicontrol(S.fh,'Style','text','String','shiftsamples','Position',[595-off1-off2 30 100 20],'HorizontalAlignment','left');
        
        
        % 3rd column
        off3=-20;
        % crossprofileshift
        S.cps=uicontrol(S.fh,'Style','checkbox','String','CrossProfileShift','FontWeight','bold','Position',[700+off3-off1 300 200 15],'Value',0,'Callback',@crossps_call);
        
        % bandpass
        S.bandpass=uicontrol(S.fh,'Style','checkbox','String','Bandpass','FontWeight','bold','Position',[700+off3-off1 270 200 15],'Value',0,'Callback',@bandpass_call);
        S.fstart=uicontrol(S.fh,'Style','edit','String','100','Enable','off','Position',[720+off3-off1 240 35 20]);
        S.fend=uicontrol(S.fh,'Style','edit','String','600','Enable','off','Position',[720+off3-off1 210 35 20]);
        S.fstarttext=uicontrol(S.fh,'Style','text','String','fstart [MHz]','Position',[760+off3-off1 240 140 20],'HorizontalAlignment','left');
        S.fendtext=uicontrol(S.fh,'Style','text','String','fend [MHz]','Position',[760+off3-off1 210 140 20],'HorizontalAlignment','left');
        
        % cut TWT
        S.cutTWT=uicontrol(S.fh,'Style','checkbox','String','Cut TWT','FontWeight','bold','Position',[700+off3-off1 180 200 15],'Value',0,'Callback',@cutt_call);
        S.tmax=uicontrol(S.fh,'Style','edit','String','80','Enable','off','Position',[720+off3-off1 150 35 20]);
        S.tmaxtext=uicontrol(S.fh,'Style','text','String','tmax [ns]','Position',[760+off3-off1 150 140 20],'HorizontalAlignment','left');
        
        % normalization
        S.norm=uicontrol(S.fh,'Style','checkbox','String','Normalization','FontWeight','bold','Position',[700+off3-off1 120 200 15],'Value',0,'Callback',@norm_call);
        S.qclip=uicontrol(S.fh,'Style','edit','String','0.98','Enable','off','Position',[720+off3-off1 90 35 20]);
        S.qcliptext=uicontrol(S.fh,'Style','text','String','qclip','Position',[760+off3-off1 90 140 20],'HorizontalAlignment','left');
        
        % k highpass
        S.khigh=uicontrol(S.fh,'Style','checkbox','String','k-highpass','FontWeight','bold','Position',[700+off3-off1 60 200 15],'Value',0,'Callback',@khigh_call);
        S.kcutoff=uicontrol(S.fh,'Style','edit','String','0.1','Enable','off','Position',[720+off3-off1 30 35 20]);
        S.kcutofftext=uicontrol(S.fh,'Style','text','String','kcutoff [1/m]','Position',[760+off3-off1 30 140 20],'HorizontalAlignment','left');
        
        % 4th column:
        off4=70;
        % gain
        S.gain=uicontrol(S.fh,'Style','checkbox','String','Gain','FontWeight','bold','Position',[900-off1-off4 300 200 15],'Value',0,'Callback',@gain_call);
        S.g1=uicontrol(S.fh,'Style','edit','String','-20','Enable','off','Position',[920-off1-off4 270 35 20]);
        S.g2=uicontrol(S.fh,'Style','edit','String','0','Enable','off','Position',[920-off1-off4 240 35 20]);
        S.g3=uicontrol(S.fh,'Style','edit','String','10','Enable','off','Position',[920-off1-off4 210 35 20]);
        S.g4=uicontrol(S.fh,'Style','edit','String','20','Enable','off','Position',[920-off1-off4 180 35 20]);
        S.g5=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[920-off1-off4 150 35 20]);
        S.g1text=uicontrol(S.fh,'Style','text','String','g1 [dB]','Position',[960-off1-off4 270 140 20],'HorizontalAlignment','left');
        S.g2text=uicontrol(S.fh,'Style','text','String','g2 [dB]','Position',[960-off1-off4 240 140 20],'HorizontalAlignment','left');
        S.g3text=uicontrol(S.fh,'Style','text','String','g3 [dB]','Position',[960-off1-off4 210 140 20],'HorizontalAlignment','left');
        S.g4text=uicontrol(S.fh,'Style','text','String','g4 [dB]','Position',[960-off1-off4 180 140 20],'HorizontalAlignment','left');
        S.g5text=uicontrol(S.fh,'Style','text','String','g5 [dB]','Position',[960-off1-off4 150 140 20],'HorizontalAlignment','left');
        
        % remove mean trace
        S.mean=uicontrol(S.fh,'Style','checkbox','String','Remove mean trace','FontWeight','bold','Position',[900-off1-off4 120 200 15],'Value',0,'Callback',@mean_call);
        S.numtrace=uicontrol(S.fh,'Style','edit','String','0','Enable','off','Position',[920-off1-off4 90 35 20]);
        S.numtracetext=uicontrol(S.fh,'Style','text','String','numtraces','Position',[960-off1-off4 90 140 20],'HorizontalAlignment','left');
        
        % remove median trace
        S.median=uicontrol(S.fh,'Style','checkbox','String','Remove median trace','FontWeight','bold','Position',[900-off1-off4 60 200 15],'Value',0,'Callback',@median_call);
        S.numtrace2=uicontrol(S.fh,'Style','edit','String','0','Enable','off','Position',[920-off1-off4 30 35 20]);
        S.numtrace2text=uicontrol(S.fh,'Style','text','String','numtraces','Position',[960-off1-off4 30 140 20],'HorizontalAlignment','left');
        
        % 5th column:
        off5=110;
        % isochrone migration 2d (reduced to 1d-v-function for easier input)
        S.migration=uicontrol(S.fh,'Style','checkbox','String','Isochrone migration v(z)','FontWeight','bold','Position',[1100-off1-off5 300 200 15],'Value',0,'Callback',@migration_call);
        S.v=uicontrol(S.fh,'Style','edit','String','0.1','Enable','off','Position',[1120-off1-off5 270 35 20]);
        S.aperture=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[1120-off1-off5 240 35 20]);
        S.vtext=uicontrol(S.fh,'Style','text','String','v [m/ns]','Position',[1160-off1-off5 270 140 20],'HorizontalAlignment','left');
        S.aperturetext=uicontrol(S.fh,'Style','text','String','aperture [°]','Position',[1160-off1-off5 240 140 20],'HorizontalAlignment','left');
        
        % Topomigration/korrektur (reduced to constant v for easier input)
        S.topo=uicontrol(S.fh,'Style','checkbox','String','Topographic migration (const. v)','FontWeight','bold','Position',[1100-off1-off5 210 200 15],'Value',0,'Callback',@topo_call);
        S.v2=uicontrol(S.fh,'Style','edit','String','0.1','Enable','off','Position',[1120-off1-off5 180 35 20]);
        S.flagtopo=uicontrol(S.fh,'Style','edit','String','1','Enable','off','Position',[1120-off1-off5 150 35 20]);
        S.aperture2=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[1120-off1-off5 120 35 20]);
        S.v2text=uicontrol(S.fh,'Style','text','String','v [m/ns]','Position',[1160-off1-off5 180 140 20],'HorizontalAlignment','left');
        S.flagtopotext=uicontrol(S.fh,'Style','text','String','correction(1), migration(2)','Position',[1160-off1-off5 150 140 20],'HorizontalAlignment','left');
        S.aperture2text=uicontrol(S.fh,'Style','text','String','aperture [°]','Position',[1160-off1-off5 120 140 20],'HorizontalAlignment','left');

        % median filter along x
        S.medfilt_x=uicontrol(S.fh,'Style','checkbox','String','X-wise median filter','FontWeight','bold','Position',[1100-off1-off5 90 200 15],'Value',0,'Callback',@medfilt_x_call);
        S.numsamp_x=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[1120-off1-off5 60 35 20]);
        S.tstart_x=uicontrol(S.fh,'Style','edit','String','0','Enable','off','Position',[1120-off1-off5 30 35 20]);
        S.numsamptext_x=uicontrol(S.fh,'Style','text','String','number of samples','Position',[1160-off1-off5 60 140 20],'HorizontalAlignment','left');
        S.tstarttext_x=uicontrol(S.fh,'Style','text','String','tstart [ns]','Position',[1160-off1-off5 30 140 20],'HorizontalAlignment','left');
        

        
        % 6th column:
        off6=120;
        % spectralWhitening
        S.sw=uicontrol(S.fh,'Style','checkbox','String','Spectral Whitening','FontWeight','bold','Position',[1300-off1-off6 300 200 15],'Value',0,'Callback',@sw_call);
        S.fmin_sw=uicontrol(S.fh,'Style','edit','String','100','Enable','off','Position',[1320-off1-off6 270 35 20]);
        S.fmin_sw_text=uicontrol(S.fh,'Style','text','String','fmin [MHz]','Position',[1360-off1-off6 270 140 20],'HorizontalAlignment','left');
        S.fmax_sw=uicontrol(S.fh,'Style','edit','String','600','Enable','off','Position',[1320-off1-off6 240 35 20]);
        S.fmax_sw_text=uicontrol(S.fh,'Style','text','String','fmax [MHz]','Position',[1360-off1-off6 240 140 20],'HorizontalAlignment','left');
        S.alpha=uicontrol(S.fh,'Style','edit','String','0.01','Enable','off','Position',[1320-off1-off6 210 35 20]);
        S.alpha_text=uicontrol(S.fh,'Style','text','String','alpha','Position',[1360-off1-off6 210 140 20],'HorizontalAlignment','left');

        off6=-80;
        % Const Trace Distance
        S.constDist=uicontrol(S.fh,'Style','checkbox','String','Constant trace distance','FontWeight','bold','Position',[1100-off1-off6 120 200 15],'Value',0,'Callback',@constDist_call);
        S.dist=uicontrol(S.fh,'Style','edit','String','0.02','Value',0,'Enable','off','Position',[1120-off1-off6 90 35 20]);
        S.disttext=uicontrol(S.fh,'Style','text','String','dx [m]','Position',[1160-off1-off6 90 140 20],'HorizontalAlignment','left');

        % Trace interpolation
        S.traceinterp=uicontrol(S.fh,'Style','checkbox','String','Trace interpolation','FontWeight','bold','Position',[1100-off1-off6 180 200 15],'Value',0,'Callback',@traceinterp_call);
        S.gap=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[1120-off1-off6 150 35 20]);
        S.gaptext=uicontrol(S.fh,'Style','text','String','number of traces','Position',[1160-off1-off6 150 140 20],'HorizontalAlignment','left');

        % median filter on traces
        S.medfilt=uicontrol(S.fh,'Style','checkbox','String','Trace-wise median filter','FontWeight','bold','Position',[1100-off1-off6 60 200 15],'Value',0,'Callback',@medfilt_call);
        S.numsamp=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[1120-off1-off6 30 35 20]);
        S.numsamptext=uicontrol(S.fh,'Style','text','String','number of samples','Position',[1160-off1-off6 30 140 20],'HorizontalAlignment','left');
    end


% save handles
guidata(S.fh,S);


%%%---------- Callback functions ----------------

    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        hei_fig=S.fh.Position(4); % figure height
        
        % change size of axes
        S.ax.Position = [300 360 wid_fig-320 hei_fig-380];
    end


    function [] = asp_call(varargin)
        % Callback for aspectratio
        S=guidata(gcbf);
        S.aspval=S.asp_String(S.asp.Value);
        S.aspval=str2num(S.aspval{1});
        set(S.ax,'DataAspectratio',[S.aspval 1 1]);
        drawnow;
        guidata(gcbf,S); % Update
    end

    function [] = colorscale_call(varargin)
        % Callback for colorscale limits
        S=guidata(gcbf);
        if isfield(S,'rad') % only for radargram, not for amplitude spectrum
            temp=get(S.rad,'CData');
            coldata=sort(unique(temp(~isnan(temp))));
            S.cmin1=coldata(round(length(coldata)/100*1));
            S.cmax1=coldata(end-round(length(coldata)/100*1));
            S.cmin3=coldata(round(length(coldata)/100*3));
            S.cmax3=coldata(end-round(length(coldata)/100*3));
            val1=S.rb1.Value;
            val2=S.rb2.Value;
            val3=S.rb3.Value;
            if val1==1 && S.flag_cs~=1
                S.rb2.Value=0;
                S.rb3.Value=0;
                set(S.ax,'ClimMode','auto');
                drawnow;
                S.flag_cs=1;
            elseif val2==1 && S.flag_cs~=2
                S.rb1.Value=0;
                S.rb3.Value=0;
                set(S.ax,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
                drawnow;
                S.flag_cs=2;
            elseif val3==1 && S.flag_cs~=3
                S.rb1.Value=0;
                S.rb2.Value=0;
                set(S.ax,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
                drawnow;
                S.flag_cs=3;
            end
        end
        guidata(gcbf,S); % Update (for flag_cs)
    end

    function [] = data_call(varargin)
        % Callback for raw/proc data
        S=guidata(gcbf);
        temp=S.channellist(S.chan.Value);
        num=str2num(temp{1}); % channel number
        if S.d1.Value==1 && S.flag_data~=1
            if isempty(num) % all channels
                set(S.rad,'CData',S.raw,'YData',S.t);
            else
                set(S.rad,'CData',S.raw(:,S.info==num),'YData',S.t);
            end
            S.flag_data=1;
            S.d2.Value=0;
        elseif S.d2.Value==1 && S.flag_data~=2
            if isempty(num) % all channels
                set(S.rad,'CData',S.proc);
            else
                set(S.rad,'CData',S.proc(:,S.info_proc==num));
            end
            S.flag_data=2;
            S.d1.Value=0;
        end
        drawnow;
        guidata(gcbf,S); % Update (for flag_data)
        plot_call();
    end

    function [] = save_call(varargin)
        % Callback for raw/proc data
        S=guidata(gcbf);
        steps=S.proclist.String; % list of processing steps
        order=1:length(steps);
        fid=fopen(fullfile(S.folder,'settings.txt'),'wt');
        if any(ismember(steps,'Remove amplitude offset'))
            fprintf(fid,['do_amplitudeOffset ',int2str(order(ismember(steps,'Remove amplitude offset'))),'\ntstart[ns] ',S.tstart.String,'\ntend[ns] ',S.tend.String,'\n\n']);
        else
            fprintf(fid,['do_amplitudeOffset 0\ntstart[ns] 0\ntend[ns] 100\n\n']);
        end
        if any(ismember(steps,'Bad trace removal'))
            fprintf(fid,['do_badTraceRemoval ',int2str(order(ismember(steps,'Bad trace removal'))),'\nminfactor ',S.minfactor.String,'\nmaxfactor ',S.maxfactor.String,'\n\n']);
        else
            fprintf(fid,['do_badTraceRemoval 0\nminfactor 2\nmaxfactor 2\n\n']);
        end
        if any(ismember(steps,'Constant trace distance'))
            fprintf(fid,['do_constTraceDist ',int2str(order(ismember(steps,'Constant trace distance'))),'\ndx[m] ',S.dist.String,'\n\n']);
        else
            fprintf(fid,['do_constTraceDist 0\ndx[m] 0.02\n\n']);
        end
        if any(ismember(steps,'Trace interpolation'))
            fprintf(fid,['do_traceInterpolation ',int2str(order(ismember(steps,'Trace interpolation'))),'\ngap ',S.gap.String,'\n\n']);
        else
            fprintf(fid,['do_traceInterpolation 0\ngap 10\n\n']);
        end
        if any(ismember(steps,'Trace-wise median filter'))
            fprintf(fid,['do_medfilt ',int2str(order(ismember(steps,'Trace-wise median filter'))),'\nnumsamp ',S.numsamp.String,'\n\n']);
        else
            fprintf(fid,['do_medfilt 0\nnumsamp 3\n\n']);
        end
        if any(ismember(steps,'X-wise median filter'))
            fprintf(fid,['do_medfilt_x ',int2str(order(ismember(steps,'X-wise median filter'))),'\nnumsamp_x ',S.numsamp_x.String,'\ntstart_x ',S.tstart_x.String,'\n\n']);
        else
            fprintf(fid,['do_medfilt_x 0\nnumsamp_x 3\ntstart_x 0\n\n']);
        end
        if any(ismember(steps,'Spherical divergence'))
            fprintf(fid,['do_sphericalDivergence ',int2str(order(ismember(steps,'Spherical divergence'))),'\n\n']);
        else
            fprintf(fid,['do_sphericalDivergence 0\n\n']);
        end
        if any(ismember(steps,'Attenuation correction'))
            fprintf(fid,['do_attenuationCorrection ',int2str(order(ismember(steps,'Attenuation correction'))),'\nsigma[S/m] ',S.sigma.String,'\neps ',S.eps.String,'\n\n']);
        else
            fprintf(fid,['do_attenuationCorrection 0\nsigma[S/m] 0.02\neps 9\n\n']);
        end
        
        if any(ismember(steps,'Spectral Whitening'))
            fprintf(fid,['do_spectralWhitening ',int2str(order(ismember(steps,'Spectral Whitening'))),'\nfmin_sw[MHz] ',S.fmin_sw.String,'\nfmax_sw[MHz] ',S.fmax_sw.String,'\nalpha ',S.alpha.String,'\n\n']);
        else
            fprintf(fid,['do_spectralWhitening 0\nfmin_sw[MHz] 100\nfmax_sw[MHz] 600\nalpha 0.01\n\n']);
        end
        
        if any(ismember(steps,'Bandpass'))
            fprintf(fid,['do_bandpass ',int2str(order(ismember(steps,'Bandpass'))),'\nfstart[MHz] ',S.fstart.String,'\nfend[MHz] ',S.fend.String,'\n\n']);
        else
            fprintf(fid,['do_bandpass 0\nfstart[MHz] 100\nfend[MHz] 600\n\n']);
        end
        if any(ismember(steps,'Isochrone migration'))
            fprintf(fid,['do_migration2d_vz ',int2str(order(ismember(steps,'Isochrone migration'))),'\nv-file[m/ns] ',S.v.String,'\ntv-file[ns] ',S.tv.String,'\naperture_m[degree] ',S.aperture.String,'\n\n']);
        else
            fprintf(fid,['do_migration2d_vz 0\nv-file[m/ns] 0.1\ntv-file[ns] \naperture_m[degree] 30\n\n']);
        end
        if any(ismember(steps,'Topomigration/correction'))
            fprintf(fid,['do_topomig2d ',int2str(order(ismember(steps,'Topomigration/correction'))),'\nv[m/ns] ',S.v2.String,'\nflag ',S.flagtopo.String,'\naperture_t[degree] ',S.aperture2.String,'\n\n']);
        else
            fprintf(fid,['do_topomig2d 0\nv[m/ns] 0.1\nflag 1\naperture_t[degree] 30\n\n']);
        end
        if any(ismember(steps,'Cut TWT'))
            fprintf(fid,['do_cutTWT ',int2str(order(ismember(steps,'Cut TWT'))),'\ntmax[ns] ',S.tmax.String,'\n\n']);
        else
            fprintf(fid,['do_cutTWT 0\ntmax[ns] 80\n\n']);
        end
        if any(ismember(steps,'k-highpass'))
            fprintf(fid,['do_khighpass ',int2str(order(ismember(steps,'k-highpass'))),'\nkcutoff[1/m] ',S.kcutoff.String,'\n\n']);
        else
            fprintf(fid,['do_khighpass 0\nkcutoff[1/m] 0.1\n\n']);
        end
        if any(ismember(steps,'Gain'))
            fprintf(fid,['do_applygain ',int2str(order(ismember(steps,'Gain'))),'\ng1 ',S.g1.String,'\ng2 ',S.g2.String,'\ng3 ',S.g3.String,'\ng4 ',S.g4.String,'\ng5 ',S.g5.String,'\n\n']);
        else
            fprintf(fid,['do_applygain 0\ng1 -20\ng2 0\ng3 10\ng4 20\ng5 30\n\n']);
        end
        if any(ismember(steps,'Normalization'))
            fprintf(fid,['do_normalization ',int2str(order(ismember(steps,'Normalization'))),'\nqclip ',S.qclip.String,'\n\n']);
        else
            fprintf(fid,['do_normalization 0\nqclip 0.98\n\n']);
        end
        if any(ismember(steps,'Remove mean trace'))
            fprintf(fid,['do_removeMeanMedianTrace ',int2str(order(ismember(steps,'Remove mean trace'))),'\nmeanmedian ',int2str(1),'\nnumtraces ',S.numtrace.String,'\n\n']);
        elseif any(ismember(steps,'Remove median trace'))
            fprintf(fid,['do_removeMeanMedianTrace ',int2str(order(ismember(steps,'Remove median trace'))),'\nmeanmedian ',int2str(2),'\nnumtraces ',S.numtrace2.String,'\n\n']);
        else
            fprintf(fid,['do_removeMeanMedianTrace 0\nmeanmedian 1\nnumtraces 0\n\n']);
        end
        if any(ismember(steps,'t0 correction'))
            if S.t0list.Value==3
                fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
                fprintf(fid,['do_t0CorrectionThreshold ',int2str(order(ismember(steps,'t0 correction'))),'\nthreshold ',S.thresh.String,'\n\n']);
                fprintf(fid,['do_t0CorrectionReferencetrace 0\nprofilenumber 0\nchannelnumber 1\ntracenumber 1\nt0[ns] 0.1\n\n']);
                if S.cps.Value==1
                    fprintf(fid,['do_crossprofileshift ',int2str(order(ismember(steps,'CrossProfileShift (only in settings.txt)'))),'\nmaxshift_cps ',S.maxshift.String,'\n\n']);
                else
                    fprintf(fid,['do_crossprofileshift 0\nmaxshift_cps 30\n\n']);
                end
                fprintf(fid,['do_channelshift 0\nOneOrAll 1\nmaxshift 30\nrefprofile 0\nshiftsamples\n']);
            elseif S.t0list.Value==4
                fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
                fprintf(fid,['do_t0CorrectionThreshold 0\nthreshold -1000\n\n']);
                fprintf(fid,['do_t0CorrectionReferencetrace ',int2str(order(ismember(steps,'t0 correction'))),'\nprofilenumber ',S.refprof.String,'\nchannelnumber ',S.refchan.String,'\ntracenumber ',S.reftrace.String,'\nt0[ns] ',S.t0.String,'\n\n']);
                if S.cps.Value==1
                    fprintf(fid,['do_crossprofileshift ',int2str(order(ismember(steps,'CrossProfileShift (only in settings.txt)'))),'\nmaxshift_cps ',S.maxshift.String,'\n\n']);
                else
                    fprintf(fid,['do_crossprofileshift 0\nmaxshift_cps 30\n\n']);
                end
                fprintf(fid,['do_channelshift 0\nOneOrAll 1\nmaxshift 30\nrefprofile 0\nshiftsamples\n']);
            elseif S.t0list.Value==1
                fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
                fprintf(fid,['do_t0CorrectionThreshold 0\nthreshold -1000\n\n']);
                fprintf(fid,['do_t0CorrectionReferencetrace 0\nprofilenumber 0\nchannelnumber 1\ntracenumber 1\nt0[ns] 0.1\n\n']);
                if S.cps.Value==1
                    fprintf(fid,['do_crossprofileshift ',int2str(order(ismember(steps,'CrossProfileShift (only in settings.txt)'))),'\nmaxshift_cps ',S.maxshift.String,'\n\n']);
                else
                    fprintf(fid,['do_crossprofileshift 0\nmaxshift_cps 30\n\n']);
                end
                fprintf(fid,['do_channelshift ',int2str(order(ismember(steps,'t0 correction'))),'\nOneOrAll ',int2str(1),'\nmaxshift ',S.maxshift.String,'\nrefprofile ',S.refprof.String,'\nshiftsamples ',S.shiftsamples.String,'\n']);
            elseif S.t0list.Value==2
                fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
                fprintf(fid,['do_t0CorrectionThreshold 0\nthreshold -1000\n\n']);
                fprintf(fid,['do_t0CorrectionReferencetrace 0\nprofilenumber 0\nchannelnumber 1\ntracenumber 1\nt0[ns] 0.1\n\n']);
                if S.cps.Value==1
                    fprintf(fid,['do_crossprofileshift ',int2str(order(ismember(steps,'CrossProfileShift (only in settings.txt)'))),'\nmaxshift_cps ',S.maxshift.String,'\n\n']);
                else
                    fprintf(fid,['do_crossprofileshift 0\nmaxshift_cps 30\n\n']);
                end
                fprintf(fid,['do_channelshift ',int2str(order(ismember(steps,'t0 correction'))),'\nOneOrAll ',int2str(2),'\nmaxshift ',S.maxshift.String,'\nrefprofile ',S.refprof.String,'\nshiftsamples ',S.shiftsamples.String,'\n']);
            end
        else
            fprintf(fid,'do_t0shift 0\nt0s[ns] 5.0\n\n');
            fprintf(fid,['do_t0CorrectionThreshold 0\nthreshold -1000\n\n']);
            fprintf(fid,['do_t0CorrectionReferencetrace 0\nprofilenumber 0\nchannelnumber 1\ntracenumber 1\nt0[ns] 0.1\n\n']);
            if S.cps.Value==1
                fprintf(fid,['do_crossprofileshift ',int2str(order(ismember(steps,'CrossProfileShift (only in settings.txt)'))),'\nmaxshift_cps ',S.maxshift.String,'\n\n']);
            else
                fprintf(fid,['do_crossprofileshift 0\nmaxshift_cps 30\n\n']);
            end
            fprintf(fid,['do_channelshift 0\nOneOrAll 1\nmaxshift 30\nrefprofile 0\nshiftsamples\n']);
        end
        fclose(fid);
        guidata(gcbf,S); % Update (for flag_data)
    end

%% Processing steps calls

    function [] = medfilt_call(varargin) % trace median filter
        S=guidata(gcbf);
        if S.medfilt.Value==1
            S.proclist.String=[S.proclist.String; {'Trace-wise median filter'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
            S.numsamp.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Trace-wise median filter'))=[];
            S.numsamp.Enable='off';
        end
        guidata(gcbf,S); % Update
    end

    function [] = medfilt_x_call(varargin) % x median filter
        S=guidata(gcbf);
        if S.medfilt_x.Value==1
            S.proclist.String=[S.proclist.String; {'X-wise median filter'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
            S.numsamp_x.Enable='on';
            S.tstart_x.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'X-wise median filter'))=[];
            S.numsamp_x.Enable='off';
            S.tstart_x.Enable='off';
        end
        guidata(gcbf,S); % Update
    end

    function [] = constDist_call(varargin)
        S=guidata(gcbf);
        if S.constDist.Value==1
            S.proclist.String=[S.proclist.String; {'Constant trace distance'}];
            S.dist.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Constant trace distance'))=[];
            S.dist.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = traceinterp_call(varargin)
        S=guidata(gcbf);
        if S.traceinterp.Value==1
            S.proclist.String=[S.proclist.String; {'Trace interpolation'}];
            S.gap.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Trace interpolation'))=[];
            S.gap.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = sw_call(varargin) % spectral whitening
        S=guidata(gcbf);
        if S.sw.Value==1
            S.proclist.String=[S.proclist.String; {'Spectral Whitening'}];
            S.fmin_sw.Enable='on';
            S.fmax_sw.Enable='on';
            S.alpha.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Spectral Whitening'))=[];
            S.fmin_sw.Enable='off';
            S.fmax_sw.Enable='off';
            S.alpha.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = migration_call(varargin)
        S=guidata(gcbf);
        if S.migration.Value==1
            S.proclist.String=[S.proclist.String; {'Isochrone migration'}];
            S.v.Enable='on';
            S.aperture.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Isochrone migration'))=[];
            S.v.Enable='off';
            S.aperture.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = topo_call(varargin)
        S=guidata(gcbf);
        if S.topo.Value==1
            S.proclist.String=[S.proclist.String; {'Topomigration/correction'}];
            S.v2.Enable='on';
            S.flagtopo.Enable='on';
            S.aperture2.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Topomigration/correction'))=[];
            S.flagtopo.Enable='off';
            S.v2.Enable='off';
            S.aperture2.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = DCrem_call(varargin)
        S=guidata(gcbf);
        if S.DCrem.Value==1
            S.proclist.String=[S.proclist.String; {'Remove amplitude offset'}];
            S.tstart.Enable='on';
            S.tend.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Remove amplitude offset'))=[];
            S.tstart.Enable='off';
            S.tend.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = badtraceremoval_call(varargin)
        S=guidata(gcbf);
        if S.badtrace.Value==1
            S.proclist.String=[S.proclist.String; {'Bad trace removal'}];
            S.minfactor.Enable='on';
            S.maxfactor.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Bad trace removal'))=[];
            S.minfactor.Enable='off';
            S.maxfactor.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = spherdiv_call(varargin)
        S=guidata(gcbf);
        if S.spherDiv.Value==1
            S.proclist.String=[S.proclist.String; {'Spherical divergence'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Spherical divergence'))=[];
        end
        guidata(gcbf,S); % Update
    end
    function [] = attenuation_call(varargin)
        S=guidata(gcbf);
        if S.attenuation.Value==1
            S.proclist.String=[S.proclist.String; {'Attenuation correction'}];
            S.sigma.Enable='on';
            S.eps.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Attenuation correction'))=[];
            S.sigma.Enable='off';
            S.eps.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = t0_call(varargin)
        S=guidata(gcbf);
        if S.t0corr.Value==1
            S.proclist.String=[S.proclist.String; {'t0 correction'}];
            S.t0list.Enable='on';
            S.thresh.Enable='on';
            S.refprof.Enable='on';
            S.refchan.Enable='on';
            S.reftrace.Enable='on';
            S.t0.Enable='on';
            S.maxshift.Enable='on';
            S.shiftsamples.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
            t0list_call();
        else
            S.proclist.String(ismember(S.proclist.String,'t0 correction'))=[];
            S.t0list.Enable='off';
            S.thresh.Enable='off';
            S.refprof.Enable='off';
            S.refchan.Enable='off';
            S.reftrace.Enable='off';
            S.t0.Enable='off';
            S.maxshift.Enable='off';
            S.shiftsamples.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = t0list_call(varargin)
        % Callback for t0list
        S=guidata(gcbf);
        if strcmp(S.t0list.String{S.t0list.Value},'Channelshift for all profiles')
            S.thresh.Enable='off';
            S.refprof.Enable='on';
            S.refchan.Enable='off';
            S.reftrace.Enable='off';
            S.t0.Enable='off';
            S.maxshift.Enable='on';
            S.shiftsamples.Enable='on';
        elseif strcmp(S.t0list.String{S.t0list.Value},'Channelshift individually')
            S.thresh.Enable='off';
            S.refprof.Enable='off';
            S.refchan.Enable='off';
            S.reftrace.Enable='off';
            S.t0.Enable='off';
            S.maxshift.Enable='on';
            S.shiftsamples.Enable='off';
        elseif strcmp(S.t0list.String{S.t0list.Value},'Amplitude threshold')
            S.thresh.Enable='on';
            S.refprof.Enable='off';
            S.refchan.Enable='off';
            S.reftrace.Enable='off';
            S.t0.Enable='off';
            S.maxshift.Enable='off';
            S.shiftsamples.Enable='off';
        elseif strcmp(S.t0list.String{S.t0list.Value},'t0 correction with reference trace')
            S.thresh.Enable='off';
            S.refprof.Enable='on';
            S.refchan.Enable='on';
            S.reftrace.Enable='on';
            S.t0.Enable='on';
            S.maxshift.Enable='off';
            S.shiftsamples.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = crossps_call(varargin)
        S=guidata(gcbf);
        if S.cps.Value==1
            S.proclist.String=[S.proclist.String; {'CrossProfileShift (only in settings.txt)'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'CrossProfileShift (only in settings.txt)'))=[];
        end
        guidata(gcbf,S); % Update
    end
    function [] = bandpass_call(varargin)
        S=guidata(gcbf);
        if S.bandpass.Value==1
            S.proclist.String=[S.proclist.String; {'Bandpass'}];
            S.fstart.Enable='on';
            S.fend.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Bandpass'))=[];
            S.fstart.Enable='off';
            S.fend.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = cutt_call(varargin)
        S=guidata(gcbf);
        if S.cutTWT.Value==1
            S.proclist.String=[S.proclist.String; {'Cut TWT'}];
            S.tmax.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Cut TWT'))=[];
            S.tmax.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = norm_call(varargin)
        S=guidata(gcbf);
        if S.norm.Value==1
            S.proclist.String=[S.proclist.String; {'Normalization'}];
            S.qclip.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Normalization'))=[];
            S.qclip.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = khigh_call(varargin)
        S=guidata(gcbf);
        if S.khigh.Value==1
            S.proclist.String=[S.proclist.String; {'k-highpass'}];
            S.kcutoff.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'k-highpass'))=[];
            S.kcutoff.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = gain_call(varargin)
        S=guidata(gcbf);
        if S.gain.Value==1
            S.proclist.String=[S.proclist.String; {'Gain'}];
            S.g1.Enable='on';
            S.g2.Enable='on';
            S.g3.Enable='on';
            S.g4.Enable='on';
            S.g5.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Gain'))=[];
            S.g1.Enable='off';
            S.g2.Enable='off';
            S.g3.Enable='off';
            S.g4.Enable='off';
            S.g5.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = mean_call(varargin)
        S=guidata(gcbf);
        if S.mean.Value==1
            S.proclist.String=[S.proclist.String; {'Remove mean trace'}];
            S.numtrace.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Remove mean trace'))=[];
            S.numtrace.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = median_call(varargin)
        S=guidata(gcbf);
        if S.median.Value==1
            S.proclist.String=[S.proclist.String; {'Remove median trace'}];
            S.numtrace2.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Remove median trace'))=[];
            S.numtrace2.Enable='off';
        end
        guidata(gcbf,S); % Update
    end


%% choose new profile or channel
    function [] = profnum_call(varargin)
        % callback for new profile
        S=guidata(gcbf);
        S.profnum=S.profilelist(S.prof.Value);
        
        % load new data
        if S.equipmentflag==1
            % Mala
            [S.traces,S.dt,S.ns,S.x,S.y,S.z]=readmala(S.folder,S.name,S.profnum);
        elseif S.equipmentflag==2
            % SPidar
            [temptraces,S.dt,S.ns,tempx,tempy,tempz]=readspidar(S.folder,S.name,S.profnum,S.ch_list,zeros(length(S.ch_list),1),zeros(length(S.ch_list),1),S.gps_channel,S.utmzone,0);
            trPerChan=length(temptraces(1,:))/length(S.ch_list);
            for i=1:length(S.ch_list)
                S.traces{i}=temptraces(:,(i-1)*trPerChan+1:i*trPerChan);
                S.x{i}=tempx(:,(i-1)*trPerChan+1:i*trPerChan)';
                S.y{i}=tempy(:,(i-1)*trPerChan+1:i*trPerChan)';
                S.z{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
            end
        elseif S.equipmentflag==3
            % Impulse Radar
            [temptraces,S.dt,S.ns,tempx,tempy,tempz,numchannels]=readImpulseRadar(S.folder,S.name,S.profnum,0,0,S.utmzone,1,5);
            trPerChan=size(temptraces,2)/numchannels;
            for i=1:numchannels
                S.traces{i}=temptraces(:,(i-1)*trPerChan+1:i*trPerChan);
                S.x{i}=tempx(:,(i-1)*trPerChan+1:i*trPerChan)';
                S.y{i}=tempy(:,(i-1)*trPerChan+1:i*trPerChan)';
                S.z{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
            end
        end
        S.raw=[];
        S.info=[];
        S.zz=[];
        for i=1:length(S.traces)
            S.raw=[S.raw S.traces{i}];
            S.xprof{i}=[0; cumsum(sqrt(diff(S.x{i}).^2+diff(S.y{i}).^2))];
            if S.equipmentflag==1 || S.equipmentflag==3
                S.info=[S.info zeros(1,length(S.xprof{i}))+i];
            elseif S.equipmentflag==2
                S.info=[S.info zeros(1,length(S.xprof{i}))+S.ch_list(i)];
            end
            S.zz=[S.zz; S.z{i}];
        end
        S.info_proc=S.info;
        S.xloc=[S.xprof{1}];
        for i=2:length(S.traces)
            S.xloc=[S.xloc; S.xprof{i}+max(S.xloc)];
        end
        S.t=0:S.dt:S.dt*(S.ns-1);
        guidata(gcbf,S); % Update
        if S.d2.Value==1 % processed data
            apply_call(); % process data
        else    % raw data
            plot_call();
        end
    end

    function [] = channum_call(varargin)
        % callback for new channel
        S=guidata(gcbf);
        guidata(gcbf,S); % Update
        plot_call();
    end

%% Plot call
    function [] = plot_call(varargin)
        % callback for plotting
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Processing test'), 'pointer', 'watch');
        temp=S.channellist(S.chan.Value);
        num=str2num(temp{1});
        child=get(S.ax,'Children');
        if length(child)>1
            for i=1:length(S.xprof)-1
                delete(child(i)); % delete lines between channels
            end
        end
        if S.d1.Value==1 % raw data
            if isempty(num)  % all channels
                set(S.rad,'CData',S.raw,'XData',S.xloc,'YData',S.t);
                for i=1:length(S.ch_list)-1
                    plot(S.ax,[max(S.xloc(S.info==S.ch_list(i))) max(S.xloc(S.info==S.ch_list(i)))],[0 max(S.t)],'r')
                end
            else
                set(S.rad,'CData',S.raw(:,S.info==num),'XData',S.xloc(S.info==num)-min(S.xloc(S.info==num)),'YData',S.t);
            end
            ylabel(S.ax,'t [ns]')
            set(S.ax,'YDir','reverse')
        else % proc data
            if isempty(num) % all channels
                if isempty(S.zmig)
                    set(S.rad,'CData',S.proc,'XData',S.xproc,'YData',S.tproc);
                    for i=1:length(S.xprof)-1
                        plot(S.ax,[max(S.xproc(S.info_proc==S.ch_list(i))) max(S.xproc(S.info_proc==S.ch_list(i)))],[0 max(S.tproc)],'r')
                    end
                    ylabel(S.ax,'t [ns]')
                    set(S.ax,'YDir','reverse')
                else
                    set(S.rad,'CData',S.proc,'XData',S.xproc,'YData',S.zmig);
                    for i=1:length(S.xprof)-1
                        plot(S.ax,[max(S.xproc(S.info_proc==S.ch_list(i))) max(S.xproc(S.info_proc==S.ch_list(i)))],[min(S.zmig) max(S.zmig)],'r')
                    end
                    ylabel(S.ax,'z [m]')
                    set(S.ax,'YDir','normal')
                end
            else % only one channel
                if isempty(S.zmig)
                    set(S.rad,'CData',S.proc(:,S.info_proc==num),'XData',S.xproc(S.info_proc==num)-min(S.xproc(S.info_proc==num)),'YData',S.tproc);
                    ylabel(S.ax,'t [ns]')
                    set(S.ax,'YDir','reverse')
                else
                    set(S.rad,'CData',S.proc(:,S.info_proc==num),'XData',S.xproc(S.info_proc==num)-min(S.xproc(S.info_proc==num)),'YData',S.zmig);
                    ylabel(S.ax,'z [m]')
                    set(S.ax,'YDir','normal')
                end
            end
        end
        if isfield(S,'rad')
            temp=get(S.rad,'CData');
            coldata=sort(unique(temp(~isnan(temp))));
            S.cmin1=coldata(round(length(coldata)/100*1));
            S.cmax1=coldata(end-round(length(coldata)/100*1));
            S.cmin3=coldata(round(length(coldata)/100*3));
            S.cmax3=coldata(end-round(length(coldata)/100*3));
            val1=S.rb1.Value;
            val2=S.rb2.Value;
            val3=S.rb3.Value;
            if val1==1
                S.rb2.Value=0;
                S.rb3.Value=0;
                set(S.ax,'ClimMode','auto');
                drawnow;
                S.flag_cs=1;
            elseif val2==1
                S.rb1.Value=0;
                S.rb3.Value=0;
                set(S.ax,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
                drawnow;
                S.flag_cs=2;
            elseif val3==1
                S.rb1.Value=0;
                S.rb2.Value=0;
                set(S.ax,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
                drawnow;
                S.flag_cs=3;
            end
        end
        asp_call();
        guidata(gcbf,S); % Update
        set(findobj('Type','Figure','Name','Processing test'), 'pointer', 'arrow');
    end

%% Apply processing steps
    function [] = apply_call(varargin)
        % callback for application of processing steps
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Processing test'), 'pointer', 'watch');
        steps=S.proclist.String; % list of processing steps
        order=1:length(steps);
        % Parameter
        params=struct('tstart',str2num(S.tstart.String),'tend',str2num(S.tend.String),...
            'tmax',str2num(S.tmax.String),'dist',str2num(S.dist.String),'gap',str2num(S.gap.String),'numsamp',str2num(S.numsamp.String),'numtraces_mean',str2num(S.numtrace.String),...
            'numtraces_median',str2num(S.numtrace2.String),'minfactor',str2num(S.minfactor.String),...
            'tstart_x',str2num(S.tstart_x.String),'numsamp_x',str2num(S.numsamp_x.String),...
            'maxfactor',str2num(S.maxfactor.String),'g1',str2num(S.g1.String),...
            'g2',str2num(S.g2.String),'g3',str2num(S.g3.String),'g4',str2num(S.g4.String),...
            'g5',str2num(S.g5.String),'threshold',str2num(S.thresh.String),...
            't0',str2num(S.t0.String),'dt',S.dt,...
            'maxshift',str2num(S.maxshift.String),'shiftsamples',str2num(S.shiftsamples.String),...
            'sigma',str2num(S.sigma.String),'eps',str2num(S.eps.String),...
            'fstart',str2num(S.fstart.String),'fend',str2num(S.fend.String),...
            'qclip',str2num(S.qclip.String),'dx',mean(diff(S.xloc)),'kcutoff',str2num(S.kcutoff.String),...
            'method',S.t0list.Value,'v',str2num(S.v.String),'v2',str2num(S.v2.String),'aperture',str2num(S.aperture.String),...
            'aperture2',str2num(S.aperture2.String),'flagtopo',str2num(S.flagtopo.String),...
            'fmin_sw',str2num(S.fmin_sw.String),'fmax_sw',str2num(S.fmax_sw.String),'alpha',str2num(S.alpha.String));
        
        if S.t0corr.Value==1 && S.t0list.Value~=2 && S.t0list.Value~=3 %only for t0 with reference trace (4) und channelshift all profiles (1)
            % prepare reference trace/profile
            % load reference trace
            if S.equipmentflag==1
                % Mala
                [temptraces,~,~,~,~,ztemp]=readmala(S.folder,S.name,str2num(S.refprof.String));
            elseif S.equipmentflag==2
                % Spidar
                [temtraces,S.dt,S.ns,tempx,tempy,tempz]=readspidar(S.folder,S.name,str2num(S.refprof.String),S.ch_list,zeros(length(S.ch_list),1),zeros(length(S.ch_list),1),S.gps_channel,S.utmzone,0);
                trPerChan=length(temtraces(1,:))/length(S.ch_list);
                if exist('temptraces','var')
                    clear temptraces;
                end
                for i=1:length(S.ch_list)
                    temptraces{i}=temtraces(:,(i-1)*trPerChan+1:i*trPerChan);
                    ztemp{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
                end
            elseif S.equipmentflag==3
                % Impulse Radar
                [temtraces,S.dt,S.ns,tempx,tempy,tempz,numchannels]=readImpulseRadar(S.folder,S.name,str2num(S.refprof.String),0,0,S.utmzone,1,5);
                trPerChan=size(temtraces,2)/numchannels;
                if exist('temptraces','var')
                    clear temptraces;
                end
                for i=1:numchannels
                    temptraces{i}=temtraces(:,(i-1)*trPerChan+1:i*trPerChan);
                    ztemp{i}=tempz(:,(i-1)*trPerChan+1:i*trPerChan)';
                end
            end
            temp=[];
            infotemp=[];
            zztemp=[];
            for i=1:length(temptraces)
                temp=[temp temptraces{i}];
                infotemp=[infotemp zeros(1,length(temptraces{i}(1,:)))+S.ch_list(i)];
                zztemp=[ztemp ztemp{i}];
            end
            % process reference profile until t0-correction
            [temp]=processing(steps(1:find(ismember(steps,'t0 correction'))-1),order(1:find(ismember(steps,'t0 correction'))-1),temp,infotemp,S.xloc,S.zz,S.t,params,[]);
            % get reference trace
            referencechannel=temp(:,str2num(S.refchan.String)==infotemp);
            referencetrace=referencechannel(:,str2num(S.reftrace.String));
            params.reftrace=referencetrace; %add reference trace to struct
            if params.method==1 && isempty(params.shiftsamples)
                % make channelshift for reference profile to determine shiftsamples
                for i=1:length(temptraces)
                    temp2{i}=temp(:,S.ch_list(i)==infotemp);
                end
                [~,shiftsamples]=channelshift(temp2,params.maxshift);
                params.shiftsamples=shiftsamples;
            end
        end
        % do processing:
        [S.proc,S.ns,S.tproc,S.zmig,S.xproc,S.info_proc]=processing(steps,order,S.raw,S.info,S.xloc,S.zz,S.t,params,S);
        % set to processed data
        S.d1.Value=0;
        S.d2.Value=1;
        S.save.Enable='on';
        S.d2.Enable='on';
        guidata(gcbf,S); % Update
        % plotting
        plot_call();
        guidata(gcbf,S); % Update
        set(findobj('Type','Figure','Name','Processing test'), 'pointer', 'arrow');
    end

%% delete list of processing steps
    function [] = del_call(varargin)
        % callback for deletion of processing steps
        S=guidata(gcbf);
        S.proclist.String=[];
        S.d1.Value=1;
        S.d2.Value=0;
        S.save.Enable='off';
        S.apply.Enable='off';
        S.d2.Enable='off';
        % turn checkboxes off
        S.DCrem.Value=0;
        S.medfilt_x.Value=0;
        S.traceinterp.Value=0;
        S.badtrace.Value=0;
        S.constDist.Value=0;
        S.spherDiv.Value=0;
        S.attenuation.Value=0;
        S.t0corr.Value=0;
        S.bandpass.Value=0;
        S.cutTWT.Value=0;
        S.norm.Value=0;
        S.khigh.Value=0;
        S.medfilt.Value=0;
        S.gain.Value=0;
        S.mean.Value=0;
        S.median.Value=0;
        S.migration.Value=0;
        S.topo.Value=0;
        S.cps.Value=0;
        S.sw.Value=0;
        % enable parameters off
        S.fmin_sw.Enable='off';
        S.fmax_sw.Enable='off';
        S.numsamp_x.Enable='off';
        S.tstart_x.Enable='off';
        S.alpha.Enable='off';
        S.v.Enable='off';
        S.gap.Enable='off';
        S.dist.Enable='off';
        S.numsamp.Enable='off';
        S.aperture.Enable='off';
        S.flagtopo.Enable='off';
        S.v2.Enable='off';
        S.aperture2.Enable='off';
        S.tstart.Enable='off';
        S.tend.Enable='off';
        S.minfactor.Enable='off';
        S.maxfactor.Enable='off';
        S.sigma.Enable='off';
        S.eps.Enable='off';
        S.t0list.Enable='off';
        S.thresh.Enable='off';
        S.refprof.Enable='off';
        S.refchan.Enable='off';
        S.reftrace.Enable='off';
        S.t0.Enable='off';
        S.maxshift.Enable='off';
        S.shiftsamples.Enable='off';
        S.fstart.Enable='off';
        S.fend.Enable='off';
        S.tmax.Enable='off';
        S.qclip.Enable='off';
        S.kcutoff.Enable='off';
        S.g1.Enable='off';
        S.g2.Enable='off';
        S.g3.Enable='off';
        S.g4.Enable='off';
        S.g5.Enable='off';
        S.numtrace.Enable='off';
        S.numtrace2.Enable='off';
        guidata(gcbf,S); % Update
        plot_call();
    end

end



%% subfunction:
function [datatraces,ns,tproc,zmig,xproc,info]=processing(steps,order,datatraces,info,x,z,t,params,S)
% (info=channel number for each trace)

% initialize values
ns=length(t);
tproc=t;
zmig=[];
xproc=x;

if ~isempty(S) % S is only given for real processing, for preparation give []
    % get parts in waitbar
    numsteps=100/length(order(order>0));
end

for k=1:length(order(order>0))  % for all processing steps in right order
    
    if ~isempty(S)
        % plot red box in waitbar
        S.wb1=patch(S.wbax,[numsteps*(k-1) numsteps*k numsteps*k numsteps*(k-1) numsteps*(k-1)],[0 0 1 1 0],[1 0 0]);
        set(S.wbax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[])
        drawnow;
    end

    if strcmp(steps{order==k},'Constant trace distance')
        if S.equipmentflag==1 || S.equipmentflag==3 % Mala or impulse Radar
            channels=max(info); % number of channels

            temp=cell(channels,1);
            xx=cell(channels,1);
            for ch=1:channels
                % create local global coords
                xyz=zeros(length(x(info==ch)),3);
                xyz(:,1)=x(info==ch);
                % make const trace dist for each channel
                minx(ch)=min(x(info==ch)); % for reducing x for each channel
                [temp{ch},xx{ch},~]=constTraceDist(datatraces(:,info==ch),params.dist,x(info==ch)-minx(ch),xyz);
            end
            numtrch=min(cellfun(@(x) length(x),xx)); % minimum number of traces per channel

            % replace new data (has different size than before!)
            info=zeros(1,channels*numtrch);
            datatraces=zeros(length(t),channels*numtrch);
            xproc=zeros(1,channels*numtrch);
            for ii=1:channels
                info((ii-1)*numtrch+1:ii*numtrch)=ii;
                xproc((ii-1)*numtrch+1:ii*numtrch)=xx{ii}(1:numtrch)+minx(ii); % add offset again for having coordinates along all channels
                datatraces(:,(ii-1)*numtrch+1:ii*numtrch)=temp{ii}(:,1:numtrch);
            end
            clear temp;
            clear xx;
        else % Spidar
            channels=length(unique(info)); % number of channels

            temp=cell(channels,1);
            xx=cell(channels,1);
            for ch=1:channels
                % create local global coords
                xyz=zeros(length(x(info==S.ch_list(ch))),3);
                xyz(:,1)=x(info==S.ch_list(ch));
                % make const trace dist for each channel
                minx(S.ch_list(ch))=min(x(info==S.ch_list(ch))); % for reducing x for each channel
                [temp{ch},xx{ch},~]=constTraceDist(datatraces(:,info==S.ch_list(ch)),params.dist,x(info==S.ch_list(ch))-minx(S.ch_list(ch)),xyz);
            end
            numtrch=min(cellfun(@(x) length(x),xx)); % minimum number of traces per channel

            % replace new data (has different size than before!)
            info=zeros(1,channels*numtrch);
            datatraces=zeros(length(t),channels*numtrch);
            xproc=zeros(1,channels*numtrch);
            for ii=1:channels
                info((ii-1)*numtrch+ii-(ii-1):ii*numtrch)=S.ch_list(ii);
                xproc((ii-1)*numtrch+ii-(ii-1):ii*numtrch)=xx{ii}(1:numtrch)+minx(S.ch_list(ii)); % add offset again for having coordinates along all channels
                datatraces(:,(ii-1)*numtrch+ii-(ii-1):ii*numtrch)=temp{ii}(:,1:numtrch);
            end
            clear temp;
            clear xx;
        end
    end

    if strcmp(steps{order==k},'Trace interpolation')
        for ch=1:length(unique(info))
            datatraces(:,info==S.ch_list(ch))=interpolation(datatraces(:,info==S.ch_list(ch)),params.gap);
        end
    end

    if strcmp(steps{order==k},'Trace-wise median filter')
        [datatraces]=medfilt(datatraces,params.numsamp);
    end

    if strcmp(steps{order==k},'X-wise median filter')
        [datatraces]=medfilt_x(datatraces,tproc,params.numsamp_x,params.tstart_x);
    end

    if strcmp(steps{order==k},'Spectral Whitening')
        for i=1:length(unique(info)) % for each channel
            datatemp(:,info==S.ch_list(i))=spectralWhitening(datatraces(:,info==S.ch_list(i)),params.dt,params.fmin_sw,params.fmax_sw,params.alpha);
        end
        datatraces=datatemp;
    end
    
    if strcmp(steps{order==k},'Topomigration/correction')
        % check for constant trace spacing:
        xtemp=x(info==S.ch_list(1));
        dx=xtemp(2)-xtemp(1); %[m] Trace distance
        if round(dx*1000)~=round(mean(diff(xtemp))*1000)
            mode.Interpreter='tex';
            mode.WindowStyle='non-modal';
            msgbox('\fontsize{15}Constant trace spacing is neccessary before topomigration! Please add "constant trace distance" (and optionally "trace interpolation") before!','Error','warn',mode);
        else
            for i=1:length(unique(info)) % for each channel
                [datatemp(:,info==S.ch_list(i)),zmig]=topomig2d_varV(datatraces(:,info==S.ch_list(i)),x(info==S.ch_list(i)),tproc,z(info==S.ch_list(i)),params.v2,params.aperture2,params.flagtopo,0,[],[]);
            end
            datatraces=datatemp;
            ns=length(datatraces(:,1));
        end
    end
    
    if strcmp(steps{order==k},'Isochrone migration')
        % read v
        v=params.v;
        for i=1:length(unique(info)) % for each channel
            % make column vectors
            tp=tproc;
            if length(tp(:,1))<length(tp(1,:))
                tp=tp';
            end
            % interpolate vgrid
            vgrid=zeros(size(datatraces(:,info==S.ch_list(i))))+v;
            % migration:
            [datatemp(:,info==S.ch_list(i)),zmig]=isochrone_mig_2d_varV(datatraces(:,info==S.ch_list(i)),x(info==S.ch_list(i)),tp,vgrid,params.aperture,0);
        end
        datatraces=datatemp;
        ns=length(datatraces(:,1));
    end
    
    if strcmp(steps{order==k},'Remove amplitude offset')
        datatraces=DCremoval(datatraces,tproc,params.tstart,params.tend);
    end
    
    if strcmp(steps{order==k},'Cut TWT')
        [datatraces,tproc,ns]=cutTWT(datatraces,tproc,params.tmax);
    end
    
    if strcmp(steps{order==k},'Remove mean trace')
        mm='mean';
        for ch=1:length(unique(info))
            datatraces(:,info==S.ch_list(ch))=removeHorizontalLines(datatraces(:,info==S.ch_list(ch)),mm,params.numtraces_mean);
        end
    end
    
    if strcmp(steps{order==k},'Remove median trace')
        mm='median';
        for ch=1:max(info)
            datatraces(:,info==ch)=removeHorizontalLines(datatraces(:,info==ch),mm,params.numtraces_median);
        end
    end
    
    if strcmp(steps{order==k},'do_badtraceremoval')
        datatraces=traceInterpolation(datatraces,params.minfactor,params.maxfactor);
    end
    
    if strcmp(steps{order==k},'Gain')
        datatraces=applygain(datatraces,[params.g1 params.g2 params.g3 params.g4 params.g5]);
    end
    
    if strcmp(steps{order==k},'t0 correction')
        if params.method==3
            % threshold
            datatraces=t0corr_thresh(datatraces,params.threshold);
        elseif params.method==4
            % t0 reference trace
            datatraces=t0correction(datatraces,params.reftrace,params.t0,params.dt);
        elseif params.method<3
            % Channelshift same for all
            % sort traces into channels
            tempch=cell(length(unique(info)),1);
            if S.equipmentflag==1 || S.equipmentflag==3
                for ch=1:length(tempch)
                    tempch{ch}=datatraces(:,info==ch);
                end
            else
                for ch=1:length(tempch)
                    tempch{ch}=datatraces(:,info==S.ch_list(ch));
                end
            end
            if params.method==1 % apply shiftsamples from reference profile or use given shiftsamples
                tempch=channelshift(tempch,params.maxshift,params.shiftsamples);
            elseif params.method==2 % apply shift individually for each profile
                tempch=channelshift(tempch,params.maxshift);
            end
            % sort traces back
            if S.equipmentflag==1 || S.equipmentflag==3
                for ch=1:length(tempch)
                    datatraces(:,info==ch)=tempch{ch};
                end
            else
                for ch=1:length(tempch)
                    datatraces(:,info==S.ch_list(ch))=tempch{ch};
                end
            end
        end
    end
    
    if strcmp(steps{order==k},'Spherical divergence')
        datatraces=sphericalDivergence(datatraces,tproc);
    end
    
    if strcmp(steps{order==k},'Attenuation correction')
        datatraces=attenuationcorrection(datatraces,tproc,params.sigma,params.eps);
    end
    
    if strcmp(steps{order==k},'Bandpass')
        datatraces=bandpass_gpr(datatraces,params.dt,params.fstart,params.fend);
    end
    
    if strcmp(steps{order==k},'Normalization')
        datatraces=normalize2d(datatraces,params.qclip);
    end
    
    if strcmp(steps{order==k},'k-highpass')
        datatraces=k_highpass(datatraces,params.dx,params.kcutoff);
    end
end

if ~isempty(S)
    % set white box in waitbar again
    set(S.wb,'XData',[0 100 100 0 0],'CData',[1 1 1]);
    set(S.wbax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[])
    % delete red boxes
    childwb=get(S.wbax,'Children');
    if length(childwb)>1
        for i=1:length(childwb)-1
            delete(childwb(i)); % delete red boxes
        end
    end
    drawnow;
end
end