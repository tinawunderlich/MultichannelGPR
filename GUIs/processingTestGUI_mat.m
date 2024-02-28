function []=processingTestGUI_mat(radargrams,t,x,global_coords,folder)

% GUI for plotting of profiles and interactive testing of processing steps
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% radargrams, t, x, global_coords, folder


S.profilelist=1:length(radargrams);
S.folder=folder;

% load first profile
S.traces=radargrams;
S.raw=radargrams{1};
S.dt=abs(t(2)-t(1));
S.ns=length(t);
S.xyz=global_coords{1};
S.gc=global_coords;
S.x=x;
S.xprof=x{1};

S.channellist=[{'All'}];

dx=mean(diff(S.xprof)); % mean trace spacing

S.t=t; % time vector

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
        S.rad=imagesc(S.xprof,S.t,S.raw);
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('t [ns]')
        set(S.ax,'DataAspectratio',[1 1 1])
        grid on
        axis tight
        
        %% UI controls:
        % UI control for channel number
        S.chan=uicontrol(S.fh,'style','listbox','unit','pix','position',[20 310 55 100],'String','1','Value',1,'callback',@channum_call);
        S.chantext=uicontrol(S.fh,'style','text','unit','pix','position',[20 420 55 15],'String','Channel','HorizontalAlignment','left');
        
        % UI control for profile number
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[80 310 55 100],'String',S.profilelist,'Value',1,'callback',@profnum_call);
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[80 420 55 15],'String','Profile','HorizontalAlignment','left');
        
        % UI control aspectratio
        S.asp_String=[{'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'} {'1/15'} {'1/20'}];
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
        
        % Const Trace Distance
        S.constDist=uicontrol(S.fh,'Style','checkbox','String','Constant trace distance','FontWeight','bold','Position',[300-off1 270 200 15],'Value',0,'Callback',@constDist_call);
        S.dist=uicontrol(S.fh,'Style','edit','String','0','Value',0,'Enable','off','Position',[320-off1 240 35 20]);
        S.disttext=uicontrol(S.fh,'Style','text','String','dx [m]','Position',[360-off1 240 140 20],'HorizontalAlignment','left');
        
        % Trace interpolation
        S.traceinterp=uicontrol(S.fh,'Style','checkbox','String','Trace interpolation','FontWeight','bold','Position',[300-off1 210 200 15],'Value',0,'Callback',@traceinterp_call);
        S.gap=uicontrol(S.fh,'Style','edit','String','3','Enable','off','Position',[320-off1 180 35 20]);
        S.gaptext=uicontrol(S.fh,'Style','text','String','number of traces','Position',[360-off1 180 140 20],'HorizontalAlignment','left');
        
        % Spherical divergence
        S.spherDiv=uicontrol(S.fh,'Style','checkbox','String','Spherical divergence','FontWeight','bold','Position',[300-off1 150 200 15],'Value',0,'Callback',@spherdiv_call);
        
        % attenuation correction
        S.attenuation=uicontrol(S.fh,'Style','checkbox','String','Attenuation correction','FontWeight','bold','Position',[300-off1 120 200 15],'Value',0,'Callback',@attenuation_call);
        S.eps=uicontrol(S.fh,'Style','edit','String','9','Enable','off','Position',[320-off1 90 35 20]);
        S.sigma=uicontrol(S.fh,'Style','edit','String','0.02','Enable','off','Position',[320-off1 60 35 20]);
        S.epstext=uicontrol(S.fh,'Style','text','String','rel. permittivity','Position',[360-off1 90 140 20],'HorizontalAlignment','left');
        S.sigmatext=uicontrol(S.fh,'Style','text','String','conductivity [S/m]','Position',[360-off1 60 140 20],'HorizontalAlignment','left');
        
        % helmert
        S.helmert=uicontrol(S.fh,'Style','checkbox','String','Coordinate transformation (no effect here)','FontWeight','bold','Position',[300-off1 30 200 15],'Value',0,'Callback',@helmert_call);


        % 2nd column
        off2=30;
        % t0 correction
        S.t0corr=uicontrol(S.fh,'Style','checkbox','String','t0 correction','FontWeight','bold','Position',[500-off1-off2 300 200 15],'Value',0,'Callback',@t0_call);
        S.t0list=uicontrol(S.fh,'Style','listbox','unit','pix','Enable','off','position',[520-off1-off2 240 170 55],'Callback',@t0list_call,'String',[{'t0 shift'}; {'Amplitude threshold'}; {'t0 correction with reference trace'}]);
        S.thresh=uicontrol(S.fh,'Style','edit','String','-1000','Enable','off','Position',[520-off1-off2 210 35 20]);
        S.t0s=uicontrol(S.fh,'Style','edit','String','5','Enable','off','Position',[520-off1-off2 180 35 20]);
        S.t0=uicontrol(S.fh,'Style','edit','String','0.5','Enable','off','Position',[520-off1-off2 150 35 20]);
        S.threshtext=uicontrol(S.fh,'Style','text','String','threshold','Position',[560-off1-off2 210 140 20],'HorizontalAlignment','left');
        S.t0stext=uicontrol(S.fh,'Style','text','String','t0s [ns]','Position',[560-off1-off2 180 140 20],'HorizontalAlignment','left');
        S.t0text=uicontrol(S.fh,'Style','text','String','t0 [ns]','Position',[560-off1-off2 150 140 20],'HorizontalAlignment','left');

        % reduce Number of samples
        S.redNumSamples=uicontrol(S.fh,'Style','checkbox','String','Reduce number of samples','FontWeight','bold','Position',[500-off1-off2 120 200 15],'Value',0,'Callback',@rednum_call);
        S.rednum=uicontrol(S.fh,'Style','edit','String','2','Enable','off','Position',[520-off1-off2 90 35 20]);
        S.rednumtext=uicontrol(S.fh,'Style','text','String','take every ? sample','Position',[560-off1-off2 90 140 20],'HorizontalAlignment','left');

        % turn profiles
        S.turn=uicontrol(S.fh,'Style','checkbox','String','Turn profiles','FontWeight','bold','Position',[500-off1-off2 60 200 15],'Value',0,'Callback',@turn_call);

        % exchange x y
        S.ex_xy=uicontrol(S.fh,'Style','checkbox','String','Exchange x and y (no effect here)','FontWeight','bold','Position',[500-off1-off2 30 200 15],'Value',0,'Callback',@ex_xy_call);


        % 3rd column
        off3=-20;
        % amplitude spectrum
        S.ampspec=uicontrol(S.fh,'Style','checkbox','String','Amplitude spectrum','FontWeight','bold','Position',[700+off3-off1 300 200 15],'Value',0,'Callback',@ampspec_call);
        
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
        S.v=uicontrol(S.fh,'Style','edit','String','','Enable','off','Position',[1120-off1-off5 270 35 20],'Callback',@readv_call);
        S.tv=uicontrol(S.fh,'Style','edit','String','','Enable','off','Position',[1120-off1-off5 240 35 20],'Callback',@readtv_call);
        S.aperture=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[1120-off1-off5 210 35 20]);
        S.vtext=uicontrol(S.fh,'Style','text','String','v [m/ns]','Position',[1160-off1-off5 270 140 20],'HorizontalAlignment','left');
        S.tvtext=uicontrol(S.fh,'Style','text','String','corresp. t [ns]','Position',[1160-off1-off5 240 140 20],'HorizontalAlignment','left');
        S.aperturetext=uicontrol(S.fh,'Style','text','String','aperture [°]','Position',[1160-off1-off5 210 140 20],'HorizontalAlignment','left');
        
        % Topomigration/korrektur (reduced to constant v for easier input)
        S.topo=uicontrol(S.fh,'Style','checkbox','String','Topographic migration (const. v)','FontWeight','bold','Position',[1100-off1-off5 180 200 15],'Value',0,'Callback',@topo_call);
        S.v2=uicontrol(S.fh,'Style','edit','String','0.1','Enable','off','Position',[1120-off1-off5 150 35 20]);
        S.flagtopo=uicontrol(S.fh,'Style','edit','String','1','Enable','off','Position',[1120-off1-off5 120 35 20]);
        S.aperture2=uicontrol(S.fh,'Style','edit','String','30','Enable','off','Position',[1120-off1-off5 90 35 20]);
        S.v2text=uicontrol(S.fh,'Style','text','String','v [m/ns]','Position',[1160-off1-off5 150 140 20],'HorizontalAlignment','left');
        S.flagtopotext=uicontrol(S.fh,'Style','text','String','correction(1), migration(2)','Position',[1160-off1-off5 120 140 20],'HorizontalAlignment','left');
        S.aperture2text=uicontrol(S.fh,'Style','text','String','aperture [°]','Position',[1160-off1-off5 90 140 20],'HorizontalAlignment','left');
        
        % 6th column:
        off6=130;
        % spectralWhitening
%         S.sw=uicontrol(S.fh,'Style','checkbox','String','Spectral Whitening','FontWeight','bold','Position',[1300-off1-off6 300 200 15],'Value',0,'Callback',@sw_call);
%         S.fmin_sw=uicontrol(S.fh,'Style','edit','String','100','Enable','off','Position',[1320-off1-off6 270 35 20]);
%         S.fmin_sw_text=uicontrol(S.fh,'Style','text','String','fmin [MHz]','Position',[1360-off1-off6 270 140 20],'HorizontalAlignment','left');
%         S.fmax_sw=uicontrol(S.fh,'Style','edit','String','600','Enable','off','Position',[1320-off1-off6 240 35 20]);
%         S.fmax_sw_text=uicontrol(S.fh,'Style','text','String','fmax [MHz]','Position',[1360-off1-off6 240 140 20],'HorizontalAlignment','left');
%         S.alpha=uicontrol(S.fh,'Style','edit','String','0.01','Enable','off','Position',[1320-off1-off6 210 35 20]);
%         S.alpha_text=uicontrol(S.fh,'Style','text','String','alpha','Position',[1360-off1-off6 210 140 20],'HorizontalAlignment','left');
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

    function [] = readv_call(varargin)
        % Callback for path and name of v-file (migration mig)
        S=guidata(gcbf);
        [S.vfile,S.vpath]=uigetfile('*.mat','Select v-vector file',S.folder);
        S.v.String=fullfile(S.vpath,S.vfile);
        guidata(gcbf,S); % Update
    end

    function [] = readtv_call(varargin)
        % Callback for path and name of tv-file (migration mig)
        S=guidata(gcbf);
        [S.tvfile,S.tvpath]=uigetfile('*.mat','Select corresponding t-vector file',S.folder);
        S.tv.String=fullfile(S.tvpath,S.tvfile);
        guidata(gcbf,S); % Update
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
        guidata(gcbf,S); % Update (for flag_cs)
    end

    function [] = data_call(varargin)
        % Callback for raw/proc data
        S=guidata(gcbf);
        temp=S.channellist(S.chan.Value);
        num=str2num(temp{1}); % channel number
        if S.d1.Value==1 && S.flag_data~=1
            set(S.rad,'CData',S.raw,'YData',S.t);
            S.flag_data=1;
            S.d2.Value=0;
        elseif S.d2.Value==1 && S.flag_data~=2
            set(S.rad,'CData',S.proc);
            S.flag_data=2;
            S.d1.Value=0;
        end
        drawnow;
        guidata(gcbf,S); % Update (for flag_data)
        plot_call();
    end

    function [] = save_call(varargin)
        % Callback for saving of settings.txt
        S=guidata(gcbf);
        steps=S.proclist.String; % list of processing steps
        order=1:length(steps);
        fid=fopen(fullfile(S.folder,'settings.txt'),'wt');
        if any(ismember(steps,'Remove amplitude offset'))
            fprintf(fid,['do_DCremoval ',int2str(order(ismember(steps,'Remove amplitude offset'))),'\n\n']);
        else
            fprintf(fid,['do_DCremoval 0\n\n']);
        end

        if any(ismember(steps,'Amplitude spectrum'))
            fprintf(fid,['do_makeAmpSpec ',int2str(order(ismember(steps,'Amplitude spectrum'))),'\n\n']);
        else
            fprintf(fid,['do_makeAmpSpec 0\n\n']);
        end

        if any(ismember(steps,'Turn profiles'))
            fprintf(fid,['do_turnProfiles ',int2str(order(ismember(steps,'Turn profiles'))),'\n\n']);
        else
            fprintf(fid,['do_turnProfiles 0\n\n']);
        end

        if any(ismember(steps,'Exchange x and y'))
            fprintf(fid,['do_exchange_x_y ',int2str(order(ismember(steps,'Exchange x and y'))),'\n\n']);
        else
            fprintf(fid,['do_exchange_x_y 0\n\n']);
        end

        if any(ismember(steps,'Coordinate transformation'))
            fprintf(fid,['do_helmertTransformation ',int2str(order(ismember(steps,'Coordinate transformation'))),'\ncoordsfile coords.txt\n\n']);
        else
            fprintf(fid,['do_helmertTransformation 0\ncoordsfile coords.txt\n\n']);
        end

        if any(ismember(steps,'Constant trace distance'))
            fprintf(fid,['do_constTraceDist ',int2str(order(ismember(steps,'Constant trace distance'))),'\ndist ',S.dist.String,'\n\n']);
        else
            fprintf(fid,['do_constTraceDist 0\ndist 0.02\n\n']);
        end

        if any(ismember(steps,'Reduce number of samples'))
            fprintf(fid,['do_reduceNumberOfSamples ',int2str(order(ismember(steps,'Reduce number of samples'))),'\nn ',S.rednum.String,'\n\n']);
        else
            fprintf(fid,['do_reduceNumberOfSamples 0\nn 2\n\n']);
        end

        if any(ismember(steps,'Trace interpolation'))
            fprintf(fid,['do_interpolation ',int2str(order(ismember(steps,'Trace interpolation'))),'\ngap ',S.gap.String,'\n\n']);
        else
            fprintf(fid,['do_interpolation 0\ngap 3\n\n']);
        end
        if any(ismember(steps,'Spherical divergence'))
            fprintf(fid,['do_sphericalDivergence ',int2str(order(ismember(steps,'Spherical divergence'))),'\n\n']);
        else
            fprintf(fid,['do_sphericalDivergence 0\n\n']);
        end
        if any(ismember(steps,'Attenuation correction'))
            fprintf(fid,['do_attenuationCorrection ',int2str(order(ismember(steps,'Attenuation correction'))),'\nsigma ',S.sigma.String,'\neps ',S.eps.String,'\n\n']);
        else
            fprintf(fid,['do_attenuationCorrection 0\nsigma 0.02\neps 9\n\n']);
        end
               
        if any(ismember(steps,'Bandpass'))
            fprintf(fid,['do_bandpass ',int2str(order(ismember(steps,'Bandpass'))),'\nfstart ',S.fstart.String,'\nfend ',S.fend.String,'\n\n']);
        else
            fprintf(fid,['do_bandpass 0\nfstart 100\nfend 600\n\n']);
        end
        if any(ismember(steps,'Isochrone migration'))
            fprintf(fid,['do_migration ',int2str(order(ismember(steps,'Isochrone migration'))),'\nvfile_mig ',S.v.String,'\naperture_mig ',S.aperture.String,'\nverbose_mig 0\n\n']);
        else
            fprintf(fid,['do_migration 0\nvfile_mig \naperture_mig 30\nverbose_mig 0\n\n']);
        end
        if any(ismember(steps,'Topomigration/correction'))
            fprintf(fid,['do_topomigration ',int2str(order(ismember(steps,'Topomigration/correction'))),'\ntopofile topo.mat\nvfile_topomig ',S.v2.String,'\naperture_topomig ',S.aperture2.String,'\nflag ',S.flagtopo.String,'\nverbose_topomig 1\nzmin \nzmax \n\n']);
        else
            fprintf(fid,['do_topomigration 0\ntopofile topo.mat\nvfile_topomig vgrid.mat\naperture_topomig 30\nflag 1\nverbose_topomig 1\nzmin \nzmax \n\n']);
        end
        if any(ismember(steps,'Cut TWT'))
            fprintf(fid,['do_cutTWT ',int2str(order(ismember(steps,'Cut TWT'))),'\ncutT ',S.tmax.String,'\n\n']);
        else
            fprintf(fid,['do_cutTWT 0\ncutT 80\n\n']);
        end
        if any(ismember(steps,'k-highpass'))
            fprintf(fid,['do_kHighpass ',int2str(order(ismember(steps,'k-highpass'))),'\nkcutoff ',S.kcutoff.String,'\n\n']);
        else
            fprintf(fid,['do_kHighpass 0\nkcutoff 0.1\n\n']);
        end

        if any(ismember(steps,'Normalization'))
            fprintf(fid,['do_normalization ',int2str(order(ismember(steps,'Normalization'))),'\nqclip ',S.qclip.String,'\n\n']);
        else
            fprintf(fid,['do_normalization 0\nqclip 0.98\n\n']);
        end
        if any(ismember(steps,'Remove mean trace'))
            fprintf(fid,['do_removeMeanTrace ',int2str(order(ismember(steps,'Remove mean trace'))),'\nmeanMedian ',int2str(1),'\nnumberOfTraces ',S.numtrace.String,'\n\n']);
        elseif any(ismember(steps,'Remove median trace'))
            fprintf(fid,['do_removeMeanTrace ',int2str(order(ismember(steps,'Remove median trace'))),'\nmeanMedian ',int2str(2),'\nnumberOfTraces ',S.numtrace2.String,'\n\n']);
        else
            fprintf(fid,['do_removeMeanTrace 0\nmeanMedian 1\nnumberOfTraces 0\n\n']);
        end
        if any(ismember(steps,'t0 correction'))
            if S.t0list.Value==2 % Threshold
                fprintf(fid,['do_t0Threshold ',int2str(order(ismember(steps,'t0 correction'))),'\nthreshold ',S.thresh.String,'\n\n']);
                fprintf(fid,['do_t0correction 0\nt0 5\n\n']);
                fprintf(fid,['do_t0shift 0\nt0s 5\n\n']);
            elseif S.t0list.Value==1 % t0 shift
                fprintf(fid,['do_t0Threshold 0\nthreshold 1000\n\n']);
                fprintf(fid,['do_t0correction 0\nt0 5\n\n']);
                fprintf(fid,['do_t0shift ',int2str(order(ismember(steps,'t0 correction'))),'\nt0s ',num2str(S.t0s.String),'\n\n']);
            elseif S.t0list.Value==3 % t0 correction (correlation with first trace)
                fprintf(fid,['do_t0Threshold 0\nthreshold 1000\n\n']);
                fprintf(fid,['do_t0correction ',int2str(order(ismember(steps,'t0 correction'))),'\nt0 ',num2str(S.t0.String),'\n\n']);
                fprintf(fid,['do_t0shift 0\nt0s 5\n\n']);
            end
        else % if no t0 correction
            fprintf(fid,['do_t0Threshold 0\nthreshold 1000\n\n']);
            fprintf(fid,['do_t0correction 0\nt0 5\n\n']);
            fprintf(fid,['do_t0shift 0\nt0s 5\n\n']);
        end
        if any(ismember(steps,'Gain'))
            fprintf(fid,['do_applyGain ',int2str(order(ismember(steps,'Gain'))),'\ng ',S.g1.String,'\t',S.g2.String,'\t',S.g3.String,'\t',S.g4.String,'\t',S.g5.String,'\n\n']);
        else
            fprintf(fid,['do_applyGain 0\ng -20\t0\t10\t20\t30\n\n']);
        end
        fclose(fid);
        guidata(gcbf,S); % Update (for flag_data)
    end

%% Processing steps calls
    function [] = turn_call(varargin) % turn profiles
        S=guidata(gcbf);
        if S.turn.Value==1
            S.proclist.String=[S.proclist.String; {'Turn profiles'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Turn profiles'))=[];
        end
        guidata(gcbf,S); % Update
    end

    function [] = ex_xy_call(varargin) % exchange x y
        S=guidata(gcbf);
        if S.ex_xy.Value==1
            S.proclist.String=[S.proclist.String; {'Exchange x and y'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Exchange x and y'))=[];
        end
        guidata(gcbf,S); % Update
    end

    function [] = ampspec_call(varargin) % amplitude spectrum
        S=guidata(gcbf);
        if S.ampspec.Value==1
            S.proclist.String=[S.proclist.String; {'Amplitude spectrum'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Amplitude spectrum'))=[];
        end
        guidata(gcbf,S); % Update
    end

    function [] = helmert_call(varargin) % helmert
        S=guidata(gcbf);
        if S.helmert.Value==1
            S.proclist.String=[S.proclist.String; {'Coordinate transformation'}];
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Coordinate transformation'))=[];
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

    function [] = rednum_call(varargin) % reduce number of samples
        S=guidata(gcbf);
        if S.redNumSamples.Value==1
            S.proclist.String=[S.proclist.String; {'Reduce number of samples'}];
            S.rednum.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Reduce number of samples'))=[];
            S.rednum.Enable='off';
        end
        guidata(gcbf,S); % Update
    end

    function [] = migration_call(varargin)
        S=guidata(gcbf);
        if S.migration.Value==1
            S.proclist.String=[S.proclist.String; {'Isochrone migration'}];
            S.v.Enable='on';
            S.tv.Enable='on';
            S.aperture.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Isochrone migration'))=[];
            S.v.Enable='off';
            S.tv.Enable='off';
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
            S.apply.Enable='on';
            S.delete.Enable='on';
        else
            S.proclist.String(ismember(S.proclist.String,'Remove amplitude offset'))=[];
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
        if S.t0corr.Value==1 % t0 shift
            S.proclist.String=[S.proclist.String; {'t0 correction'}];
            S.t0list.Enable='on';
            S.thresh.Enable='off';
            S.t0.Enable='off';
            S.t0s.Enable='on';
            S.apply.Enable='on';
            S.delete.Enable='on';
            t0list_call();
        elseif S.t0corr.Value==2 % threshold
            S.proclist.String=[S.proclist.String; {'t0 correction'}];
            S.t0list.Enable='on';
            S.thresh.Enable='on';
            S.t0.Enable='off';
            S.t0s.Enable='off';
            S.apply.Enable='on';
            S.delete.Enable='on';
            t0list_call();
        elseif S.t0corr.Value==3 % correlation with reference trace
            S.proclist.String=[S.proclist.String; {'t0 correction'}];
            S.t0list.Enable='on';
            S.thresh.Enable='off';
            S.t0.Enable='on';
            S.t0s.Enable='off';
            S.apply.Enable='on';
            S.delete.Enable='on';
            t0list_call();
        else
            S.proclist.String(ismember(S.proclist.String,'t0 correction'))=[];
            S.t0list.Enable='off';
            S.thresh.Enable='off';
            S.t0.Enable='off';
            S.t0s.Enable='off';
        end
        guidata(gcbf,S); % Update
    end
    function [] = t0list_call(varargin)
        % Callback for t0list
        S=guidata(gcbf);
        if strcmp(S.t0list.String{S.t0list.Value},'Amplitude threshold')
            S.thresh.Enable='on';
            S.t0.Enable='off';
            S.t0s.Enable='off';
        elseif strcmp(S.t0list.String{S.t0list.Value},'t0 shift')
            S.thresh.Enable='off';
            S.t0.Enable='off';
            S.t0s.Enable='on';
        elseif strcmp(S.t0list.String{S.t0list.Value},'t0 correction with reference trace')
            S.thresh.Enable='off';
            S.t0.Enable='on';
            S.t0s.Enable='off';
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
        S.raw=S.traces{S.profnum};
        S.xprof=S.x{S.profnum};
        S.xyz=S.gc{S.profnum};
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

        child=get(S.ax,'Children');
        if length(child)>1
            for i=1:length(S.xprof)-1
                delete(child(i)); % delete lines between channels
            end
        end
        if S.d1.Value==1 % raw data
            set(S.rad,'CData',S.raw,'XData',S.xprof,'YData',S.t);
            ylabel(S.ax,'t [ns]')
            set(S.ax,'YDir','reverse')
        else % proc data
            if isempty(S.zmig) && isempty(S.f) % normal time data
                set(S.rad,'CData',S.proc,'XData',S.xprof,'YData',S.tproc);
                %S.rad=imagesc(S.ax,S.xprof,S.tproc,S.proc);
                ylabel(S.ax,'t [ns]')
                set(S.ax,'YDir','reverse')
            elseif ~isempty(S.zmig) && isempty(S.f) % z-mig data
                set(S.rad,'CData',S.proc,'XData',S.xprof,'YData',S.zmig);
                %S.rad=imagesc(S.ax,S.xprof,S.zmig,S.proc);
                ylabel(S.ax,'z [m]')
                set(S.ax,'YDir','normal')
            else % amplitude spectrum
                plot(S.ax,S.f,S.ampl);
                ylabel(S.ax,'Amplitudes')
                xlabel(S.ax,'f [MHz]')
                set(S.ax,'YDir','normal')
            end
        end
        drawnow;
        guidata(gcbf,S); % Update
        colorscale_call();
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
        params=struct('tmax',str2num(S.tmax.String),'numtraces_mean',str2num(S.numtrace.String),...
            'numtraces_median',str2num(S.numtrace2.String),'gap',str2num(S.gap.String),...
            'g1',str2num(S.g1.String),...
            'g2',str2num(S.g2.String),'g3',str2num(S.g3.String),'g4',str2num(S.g4.String),...
            'g5',str2num(S.g5.String),'threshold',str2num(S.thresh.String),...
            't0',str2num(S.t0.String),'t0s',str2num(S.t0s.String),'dt',S.dt,...
            'sigma',str2num(S.sigma.String),'eps',str2num(S.eps.String),...
            'fstart',str2num(S.fstart.String),'fend',str2num(S.fend.String),...
            'qclip',str2num(S.qclip.String),'dx',mean(diff(S.xprof)),'kcutoff',str2num(S.kcutoff.String),...
            'method',S.t0list.Value,'v',S.v.String,'v2',str2num(S.v2.String),'aperture',str2num(S.aperture.String),...
            'aperture2',str2num(S.aperture2.String),'flagtopo',str2num(S.flagtopo.String),'tv',S.tv.String,...
            'dist',str2num(S.dist.String),'rednum',str2num(S.rednum.String));
        
        % do processing:
        [S.proc,S.ns,S.tproc,S.zmig,S.f,S.ampl]=processing(steps,order,S.raw,S.xprof,S.xyz,S.t,params,S);
        % set to processed data
        S.d1.Value=0;
        S.d2.Value=1;
        S.save.Enable='on';
        S.d2.Enable='on';
        guidata(gcbf,S); % Update
        % plotting
        plot_call();
        guidata(gcbf,S); % Update
        colorscale_call();
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
        S.traceinterp.Value=0;
        S.spherDiv.Value=0;
        S.attenuation.Value=0;
        S.t0corr.Value=0;
        S.bandpass.Value=0;
        S.cutTWT.Value=0;
        S.norm.Value=0;
        S.khigh.Value=0;
        S.gain.Value=0;
        S.mean.Value=0;
        S.median.Value=0;
        S.migration.Value=0;
        S.topo.Value=0;
        S.constDist.Value=0;
        S.redNumSamples.Value=0;
        S.ampspec.Value=0;
        S.turn.Value=0;
        S.ex_xy.Value=0;
        S.helmert.Value=0;
        % enable parameters off
        S.fmin_sw.Enable='off';
        S.fmax_sw.Enable='off';
        S.alpha.Enable='off';
        S.v.Enable='off';
        S.tv.Enable='off';
        S.aperture.Enable='off';
        S.flagtopo.Enable='off';
        S.v2.Enable='off';
        S.aperture2.Enable='off';
        S.gap.Enable='off';
        S.sigma.Enable='off';
        S.eps.Enable='off';
        S.t0list.Enable='off';
        S.thresh.Enable='off';
        S.t0.Enable='off';
        S.t0s.Enable='off';
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
        S.dist.Enable='off';
        S.rednum.Enable='off';
   
        guidata(gcbf,S); % Update
        plot_call();
    end

end



%% subfunction:
function [datatraces,ns,tproc,zmig,f,ampl]=processing(steps,order,datatraces,x,xyz,t,params,S)

ns=length(t);
tproc=t;
zmig=[];
f=[];
ampl=[];

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
    
    if strcmp(steps{order==k},'Turn profiles')
        [datatraces,x,xyz]=turnProfiles(datatraces,x,xyz);
    end

    if strcmp(steps{order==k},'Constant trace distance')
        [datatraces,x,xyz]=constTraceDist(datatraces,params.dist,x,xyz);
    end

    if strcmp(steps{order==k},'Reduce number of samples')
        [datatraces,tproc]=reduceNumberOfSamples(datatraces,t,params.rednum);
    end

    if strcmp(steps{order==k},'Amplitude spectrum')
        [f,ampl]=makeAmpspec(t,datatraces);
    end
    
    if strcmp(steps{order==k},'Topomigration/correction')
        [datatemp,zmig]=topomig2d_varV(datatraces,x,tproc,xyz(:,3),params.v2,params.aperture2,params.flagtopo,0,[],[]);
        datatraces=datatemp;
        ns=length(datatraces(:,1));
    end
    
    if strcmp(steps{order==k},'Isochrone migration')
        % read v(t)
        temp=load(params.v);
        temp2=fieldnames(temp);
        v=getfield(temp,temp2{1});
        temp=load(params.tv);
        temp2=fieldnames(temp);
        tv=getfield(temp,temp2{1});

        % make column vectors
        if length(v(:,1))<length(v(1,:))
            v=v';
        end
        if length(tv(:,1))<length(tv(1,:))
            tv=tv';
        end
        tp=tproc;
        if length(tp(:,1))<length(tp(1,:))
            tp=tp';
        end
        % interpolate vgrid
        vgrid=repmat(interp1(tv,v,tp),[1 length(datatraces(:,1))]);
        % migration:
        [datatemp,zmig]=isochrone_mig_2d_varV(datatraces,x,tp,vgrid,params.aperture,0);
        datatraces=datatemp;
        ns=length(datatraces(:,1));
    end
    
    if strcmp(steps{order==k},'Remove amplitude offset')
        datatraces=DCremoval(datatraces,tproc,0,max(t));
    end
    
    if strcmp(steps{order==k},'Cut TWT')
        [datatraces,tproc,ns]=cutTWT(datatraces,t,params.tmax);
    end
    
    if strcmp(steps{order==k},'Remove mean trace')
        mm='mean';
        datatraces=removeHorizontalLines(datatraces,mm,params.numtraces_mean);
    end
    
    if strcmp(steps{order==k},'Remove median trace')
        mm='median';
        datatraces=removeHorizontalLines(datatraces,mm,params.numtraces_median);
    end
    
    if strcmp(steps{order==k},'do_traceInterpolation')
        datatraces=traceInterpolation(datatraces,params.minfactor,params.maxfactor);
    end
    
    if strcmp(steps{order==k},'Gain')
        datatraces=applygain(datatraces,[params.g1 params.g2 params.g3 params.g4 params.g5]);
    end
    
    if strcmp(steps{order==k},'t0 correction')
        if params.method==2
            %% threshold
            datatraces=t0corr_thresh(datatraces,params.threshold);
        elseif params.method==1
            %% t0 shift
            datatraces=t0shift(datatraces,params.t0s,params.dt);
        elseif params.method==3
            %% t0 correction (correlation)
            datatraces=t0correction(datatraces,datatraces(:,1),params.t0,params.dt)
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