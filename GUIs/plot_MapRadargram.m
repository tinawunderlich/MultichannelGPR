function []=plot_MapRadargramm(radar,t,coords,profileslist,timedepth)

% GUI for plotting of map and radargrams with connected cursor
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% (created by make_Radargrams.m, Bin2Radargram.m or own creation)
% radar: cells with data of radargrams 
% t: time vector of radargrams [ns]
% coords: cells with starting and ending coordinates of profile (or
% complete profile lists for curved profiles)=global_coords
% profilelist: number of profiles
% timedepth: flag: time=1, depth=2


S.timedepth=timedepth;

% set radargramms (=data)
S.data=radar;
S.t=t;
S.maxElevation=S.t(1);
S.dt=abs(S.t(2)-S.t(1));
for i=1:length(coords)
    S.xprof{i}=coords{i}(:,1);  % x-coordinate along profile
    S.yprof{i}=coords{i}(:,2);  % y-coordinate along profile
    S.loc_coord{i}=[0; cumsum(sqrt(diff(S.xprof{i}).^2+diff(S.yprof{i}).^2))];  % local coordinate along profile
    % delete possibly same coordinates:
    [S.loc_coord{i},a]=unique(S.loc_coord{i});
    S.xprof{i}=S.xprof{i}(a);
    S.yprof{i}=S.yprof{i}(a);
    S.data{i}=S.data{i}(:,a);
end
S.profileslist=profileslist;
S.rdgnum='1'; % radargram number
S.ar=1/10;


% colorscale limits radargrams
for i=1:length(S.data)
    coldata=sort(unique(S.data{i}(~isnan(S.data{i}(:)))));
    S.cmin1R(i)=coldata(round(length(coldata)/100*1));
    S.cmax1R(i)=coldata(end-round(length(coldata)/100*1));
    S.cmin3R(i)=coldata(round(length(coldata)/100*3));
    S.cmax3R(i)=coldata(end-round(length(coldata)/100*3));
end

% get screensize
S.siz = get( 0, 'Screensize' );
S.fh.Position=S.siz;  % default: screensize

% Add the UI components
S=addcomponents(S);

% Make figure visible after adding components
S.fh.Visible='on';

    function S=addcomponents(S)
        
        % create figure handle: (Renderer is important for scrollbar)
        S.fh = figure('units','pixels',...
            'position',S.siz,...
            'menubar','none',...
            'name','Map and Radargrams',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui,'Renderer','painters','KeyPressFcn',@keypress);
        
        S.fh.MenuBar = 'figure';
        set (S.fh, 'WindowButtonMotionFcn', @mouseMove); % for connecting of Tsl and radargrams
        
        %%% RADARGRAM
        % scrollbar for radargram:
        S.panelR=uipanel('Parent',S.fh);
        S.panelRa=uipanel('Parent',S.panelR);
        if range(S.yprof{1})>range(S.xprof{1}) % Map height>width
            % layout1
            set(S.panelR,'Position',[180+(S.siz(3)-200)/3 40 S.siz(3)-200-(S.siz(3)-200)/3 S.siz(4)/3],'Units','pixels'); % panel for cutting frame with slider
            set(S.panelRa,'Position',[180+(S.siz(3)-200)/3 40 S.siz(3)-200-(S.siz(3)-200)/3 S.siz(4)/3],'Units','pixels'); % panel inside cutting frame with image
        else
            % layout2
            set(S.panelR,'Position',[150 40 S.siz(3)-170 (S.siz(4)-90)/3],'Units','pixels'); % panel for cutting frame with slider
            set(S.panelRa,'Position',[150 40 S.siz(3)-170 (S.siz(4)-90)/3],'Units','pixels'); % panel inside cutting frame with image
        end
        S.sliderR=uicontrol('Style','Slider','Parent',S.fh,'Units','pixels','Position',[150 160 S.siz(3)-100 0.05],'Value',0.5,'Callback',@sliderR_callback);

        S.ax_data=axes('units','normalized',...
                'OuterPosition',[0.02 0.02 0.98 0.98]);
        
        % make initial plot RADARGRAM
        S.h_rdg=imagesc(S.ax_data,S.loc_coord{1},S.t,S.data{1});
        colormap(flipud(gray));
        xlabel('ProfileX [m]')
        set(S.ax_data,'Xlim',[min(S.loc_coord{1}) max(S.loc_coord{1})],'CLim',[S.cmin3R(1) S.cmax3R(1)])
        if S.timedepth==1
            ylabel('t [ns]')
            axis ij
        elseif S.timedepth==2
            ylabel('z [m]')
            axis xy
        end
        if S.timedepth==1
            set(S.ax_data,'DataAspectRatio',[1 1/S.ar 1]);
        else
            set(S.ax_data,'DataAspectRatio',[1 1 1]);
        end
        grid on
        axis tight


        %%% Map
        % scrollbars for Map:
        S.panelT=uipanel('Parent',S.fh);
        S.panelTa=uipanel('Parent',S.panelT);
        if range(S.yprof{1})>range(S.xprof{1}) % Tsl height>width
            % layout 1
            set(S.panelT,'Position',[150 60 (S.siz(3)-200)/3 S.siz(4)-60],'Units','pixels'); % panel for cutting frame with slider
            set(S.panelTa,'Position',[150 60 (S.siz(3)-200)/3 S.siz(4)-60],'Units','pixels'); % panel inside cutting frame with image
            S.layout=1;
        else
            % layout 2
            set(S.panelT,'Position',[150 S.siz(4)-(S.siz(4)-90)/3*2 S.siz(3)-170 (S.siz(4)-90)/3*2],'Units','pixels'); % panel for cutting frame with slider
            set(S.panelTa,'Position',[150 S.siz(4)-(S.siz(4)-90)/3*2 S.siz(3)-170 (S.siz(4)-90)/3*2],'Units','pixels'); % panel inside cutting frame with image
            S.layout=2;
        end
        % vertical slider:
        S.sliderTV=uicontrol('Style','Slider','Parent',S.fh,'Units','pixels','Position',[150 160 0.05 S.siz(4)-100],'Value',0.5,'Callback',@sliderT_callback);
        % horizontal slider
        S.sliderTH=uicontrol('Style','Slider','Parent',S.fh,'Units','pixels','Position',[150 160 S.siz(3)-100 0.05],'Value',0.5,'Callback',@sliderT_callback);
        
        S.ax_tsl = axes('units','normalized',...
            'OuterPosition',[0.02 0.02 0.98 0.98],'Parent',S.panelTa);
        % make initial plot map
        S.ax_tsl.NextPlot='add';
        for i=1:length(S.loc_coord)
            S.h_tsl=plot(S.ax_tsl,S.xprof{i},S.yprof{i},'b','Linewidth',2);
        end
        S.line=plot(S.ax_tsl,S.xprof{1},S.yprof{1},'r','Linewidth',2);
        xlabel('x [m]')
        ylabel('y [m]')
        axis xy
        set(S.ax_tsl,'DataAspectratio',[1 1 1])
        grid on
        axis tight
        

        
        % UI control for aspect ratio
        S.aspect=uicontrol(S.fh,'Style','popupmenu','unit','pix','position',[5 195 110 25],'Value',12,'String',[{'AspectRatio'} {'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'} {'1/20'}],'callback',{@aspect_call});
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[10 225 100 15],'String','Aspect ratio','HorizontalAlignment','left','FontWeight','bold');
        if S.timedepth==2
            S.aspect.Value=7;
        end

        % UI control for tsl zoom
        S.zoom=uicontrol(S.fh,'Style','popupmenu','unit','pix','position',[5 500 110 25],'Value',5,'String',[{'10 %'} {'20 %'} {'50 %'} {'75 %'} {'100 %'} {'150 %'} {'175 %'} {'200 %'} {'300 %'}],'callback',{@zoom_call});
        

        % UI control for profile number 
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[10 80 100 85],'String',S.profileslist,'Value',1,'callback',{@profnum_call});
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[10 170 100 15],'String','Choose profile','HorizontalAlignment','left','FontWeight','bold');

        
        % UI control - radiobuttons for layout
        S.lay1=uicontrol('Style','radiobutton','String','Layout 1','Position',[10 35 100 15],'Callback',@layout_call); 
        S.lay2=uicontrol('Style','radiobutton','String','Layout 2','Position',[10 15 100 15],'Callback',@layout_call);
        if S.layout==1
            S.lay1.Value=1;
        else
            S.lay2.Value=1;
        end
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[10 365 100 15],'Callback',@colorscale_call);
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[10 345 100 15],'Callback',@colorscale_call);
        S.rb3=uicontrol('Style','radiobutton','String','3 %','Position',[10 325 100 15],'Callback',@colorscale_call,'Value',1);  % this button is on
        S.flag_cs=3;
        
        % UI control for wiggleplot
        S.cb_wiggle=uicontrol('Style','checkbox','String','Wiggle plot','Position',[10 275 100 30],'Value',0,'Callback',@wiggle_call);
        if S.timedepth==2
            S.cb_wiggle.Enable='off';
        end
        S.wigscal=uicontrol('Style','Edit','String','1','Enable','off','Position',[10 255 45 25],'Callback',@wiggle_call);
        S.wigscaltext=uicontrol('Style','text','String','Scale','Position',[52 255 50 25]);
        
        % UI control for depth scale switch
        if S.timedepth==2
            S.dswitch=uicontrol('Style','checkbox','String','Absolute depth','Position',[560 10 180 30],'Value',1,'Callback',@dswitch_call);
        end

        % UI control for info text
        S.infotext=uicontrol('Style','text','String','For scrolling: move mouse inside plot and use arrow keys','Position',[800 10 200 30]);
    end


% save handles
guidata(S.fh,S);


%%%---------- Callback functions ----------------


    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        hei_fig=S.fh.Position(4); % figure height
        
        % change size of axes
        set(S.ax_data,'OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panelRa);
        set(S.ax_tsl,'OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panelTa);

        if S.layout==1 % Tsl height>width: Layout1
            % change size of panels
            S.panelR.Position=[180+(wid_fig-200)/3 60 wid_fig-200-(wid_fig-200)/3 hei_fig/3];
            S.panelT.Position=[150 60 (wid_fig-200)/3 hei_fig-60];
            % change slider:
            S.sliderR.Position=[180+(wid_fig-200)/3 60 wid_fig-200-(wid_fig-200)/3 0.05];
            S.sliderTH.Position=[150 60 (wid_fig-200)/3 0.05];
            S.sliderTV.Position=[150 60 0.05 hei_fig-60];
        else  % Layout2
            % change size of panels
            S.panelR.Position=[150 60 wid_fig-170 (hei_fig-90)/3];
            S.panelT.Position=[150 hei_fig-(hei_fig-90)/3*2 wid_fig-170 (hei_fig-90)/3*2];
            % change slider:
            S.sliderR.Position=[150 60 wid_fig-170 0.05];
            S.sliderTH.Position=[150 hei_fig-(hei_fig-90)/3*2 wid_fig-170 0.05];
            S.sliderTV.Position=[150 hei_fig-(hei_fig-90)/3*2 0.05 (hei_fig-90)/3*2];
        end

        % slider update
        % radargram
        val1 = get(S.sliderR,'Value');
        set(S.panelRa,'Position',[-val1 0 2 1],'Units','normalized');
        % tsl
        temp=extractBefore(S.zoom.String{S.zoom.Value},' %'); % e.g. '100 %'
        S.factor=1/100*str2num(temp); % scaling factor for panelTa
        val2 = get(S.sliderTH,'Value');
        val3 = get(S.sliderTV,'Value');
        set(S.panelTa,'Position',[-val2 -val3 2 2]*S.factor,'Units','normalized');
    end


    function keypress(hObject,event)
        pos=event.Source.CurrentPoint; % current point in figure
        pos_rad_panel=S.panelR.Position; % position of radargram-panel
        pos_tsl_panel=S.panelT.Position; % position of tsl-panel
       
        % find out if mouse is inside a plot
        if pos(1)>=pos_rad_panel(1) && pos(1)<=pos_rad_panel(1)+pos_rad_panel(3) && pos(2)>=pos_rad_panel(2) && pos(2)<=pos_rad_panel(2)+pos_rad_panel(4)
            % Mouse in radargram-panel
            val = get(S.sliderR,'Value');
            if strcmp(event.Key,'rightarrow')
                set(S.panelRa,'Position',[-(val+0.05) 0 2 1],'Units','normalized')
                S.sliderR.Value=val+0.05;
            elseif strcmp(event.Key,'leftarrow')
                set(S.panelRa,'Position',[-(val-0.05) 0 2 1],'Units','normalized')
                S.sliderR.Value=val-0.05;
            end
        elseif pos(1)>=pos_tsl_panel(1) && pos(1)<=pos_tsl_panel(1)+pos_tsl_panel(3) && pos(2)>=pos_tsl_panel(2) && pos(2)<=pos_tsl_panel(2)+pos_tsl_panel(4)
            % Mouse in Tsl-panel
            val2 = get(S.sliderTH,'Value');
            val3 = get(S.sliderTV,'Value');
            temp=extractBefore(S.zoom.String{S.zoom.Value},' %'); % e.g. '100 %'
            S.factor=1/100*str2num(temp); % scaling factor for panelTa
            
            if strcmp(event.Key,'rightarrow') % horizontal slider:
                set(S.panelTa,'Position',[-(val2+0.05)*S.factor -val3*S.factor 2 2]*S.factor,'Units','normalized')
                S.sliderTH.Value=val2+0.05;
            elseif strcmp(event.Key,'leftarrow')
                set(S.panelTa,'Position',[-(val2-0.05)*S.factor -val3*S.factor 2 2]*S.factor,'Units','normalized')
                S.sliderTH.Value=val2-0.05;   
            elseif strcmp(event.Key,'uparrow')  % vertical slider:
                set(S.panelTa,'Position',[-val2*S.factor -(val3+0.05)*S.factor 2 2]*S.factor,'Units','normalized')
                S.sliderTV.Value=val3+0.05;
            elseif strcmp(event.Key,'downarrow')
                set(S.panelTa,'Position',[-val2*S.factor -(val3-0.05)*S.factor 2 2]*S.factor,'Units','normalized')
                S.sliderTV.Value=val3-0.05;
            end
        end
    end

    function sliderR_callback(varargin) % for scrollbar radargram
        S=guidata(gcbf);
        val = get(S.sliderR,'Value');
        set(S.panelRa,'Position',[-val 0 2 1],'Units','normalized')
        guidata(gcbf,S); % Update
    end

    function sliderT_callback(varargin) % for scrollbars tsl 
        S=guidata(gcbf);
        % zoom scaling factor:
        temp=extractBefore(S.zoom.String{S.zoom.Value},' %'); % e.g. '100 %'
        S.factor=1/100*str2num(temp); % scaling factor for panelTa
        % slider values:
        val2 = get(S.sliderTH,'Value');
        val3 = get(S.sliderTV,'Value');
        set(S.panelTa,'Position',[-val2*S.factor -val3*S.factor 2 2]*S.factor,'Units','normalized');
        guidata(gcbf,S); % Update
    end


    function mouseMove (object, eventdata)
        % get position inside radargram
        C = get (S.ax_data, 'CurrentPoint');
        S=guidata(object);
        
        % if mouse inside radargram:
        if C(1,1)>=S.ax_data.XLim(1) && C(1,1)<=S.ax_data.XLim(2) && C(1,2)>=S.ax_data.YLim(1) && C(1,2)<=S.ax_data.YLim(2)
            set(gcf,'pointer','crosshair');

            % find global coords of point:
            S.xy=[interp1(S.loc_coord{str2num(S.rdgnum)},S.xprof{str2num(S.rdgnum)},C(1,1)) interp1(S.loc_coord{str2num(S.rdgnum)},S.yprof{str2num(S.rdgnum)},C(1,1))];

            % update cross in map-plot:
            obj=get(S.ax_tsl,'Children');
            nobj=length(obj);
            for o=1:nobj-length(S.data)-1
                delete(obj(o)); % delete old lines
            end
            dd=(S.ax_tsl.XLim(2)-S.ax_tsl.XLim(1))/100;  % half width of cross
            plot(S.ax_tsl,[S.xy(1)-dd S.xy(1)+dd],[S.xy(2) S.xy(2)],'k','Linewidth',2);
            plot(S.ax_tsl,[S.xy(1) S.xy(1)],[S.xy(2)-dd S.xy(2)+dd],'k','Linewidth',2);

            % update title with current coordinates:
            if S.timedepth==1
                title(S.ax_data, ['(ProfileX,T) = (', num2str(C(1,1),'%8.2f'), 'm, ',num2str(C(1,2),'%8.2f'), 'ns)','   (X,Y) = (', num2str(S.xy(1),'%8.2f'), 'm, ',num2str(S.xy(2),'%8.2f'), 'm)']); 
            else
                title(S.ax_data, ['(ProfileX,Z) = (', num2str(C(1,1),'%8.2f'), 'm, ',num2str(C(1,2),'%8.2f'), 'm)','   (X,Y) = (', num2str(S.xy(1),'%8.2f'), 'm, ',num2str(S.xy(2),'%8.2f'), 'm)']); 
            end
        else
            set(gcf,'pointer','arrow');
        end
        guidata(gcbf,S); % Update
    end



    function [] = profnum_call(varargin)
        % new profile from list
        S=guidata(gcbf);
        S.rdgnum=S.profileslist(S.prof.Value,:);
        guidata(gcbf,S); % Update
        % plot new red line in tsl-plot:
        newline(varargin);

        % delete all old plotting things
        obj=get(S.ax_data,'Children');
        nobj=length(obj);
        for o=1:nobj-length(S.data)
            delete(obj(o)); % delete old lines
        end
        if S.cb_wiggle.Value==0
            %set(S.h_rdg,'XData',S.loc_coord{str2num(S.rdgnum)},'CData',S.data{str2num(S.rdgnum)});
            S.h_rdg=imagesc(S.ax_data,S.loc_coord{str2num(S.rdgnum)},S.t,S.data{str2num(S.rdgnum)});
            colormap(flipud(gray));
        else
            % Wiggle plot
            axes(S.ax_data); % activate current axes
            wigglesc(S.data{str2num(S.rdgnum)},S.t,S.loc_coord{str2num(S.rdgnum)},str2num(S.wigscal.String));
        end
        set(S.ax_data,'Xlim',[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})]);
        axes(S.ax_data);
        if S.timedepth==1
            ylabel('t [ns]')
        else
            ylabel('z [m]')
        end
        xlabel('ProfileX [m]')
        if S.timedepth==2 && S.dswitch.Value==1
            axis xy
        elseif S.timedepth==2 && S.dswitch.Value==0
            axis ij
        end

        % colorscale
        nR=str2num(S.rdgnum); % radargram-number
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1
            set(S.ax_data,'ClimMode','auto');
            drawnow;
        elseif val2==1
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin1R(nR) S.cmax1R(nR)]);
            drawnow;
        elseif val3==1
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin3R(nR) S.cmax3R(nR)]);
            drawnow;
        end

        guidata(gcbf,S); % Update
        aspect_call();
    end

    function [] = dswitch_call(varargin)
        % for depth-switch
        S=guidata(gcbf);
        if S.dswitch.Value==1
            % absolute depth
            S.t=S.maxElevation-S.t;
            axes(S.ax_data);
            axis xy
        else
            % relative depth from top of radargram
            S.t=abs(S.t-S.maxElevation);
            axes(S.ax_data);
            axis ij
        end
        guidata(gcbf,S); % Update
        profnum_call();
    end


    function [] = layout_call(varargin)
        % Callback for layout
        S=guidata(gcbf);
        if strcmp(varargin{2}.Source.String,'Layout 1')
            % Tsl and radargram next to each other
            S.layout=1;
            S.lay2.Value=0;
        else
            % TSL and radargram over each other
            S.layout=2;
            S.lay1.Value=0;
        end
        resizeui();
        guidata(gcbf,S); % Update
    end

    function [] = wiggle_call(varargin)
        % Callback for wiggle plot of radargram
        S=guidata(gcbf);
        if S.cb_wiggle.Value==1
            S.wigscal.Enable='on';
        else
            S.wigscal.Enable='off';
        end
        guidata(gcbf,S); % Update
        profnum_call();
    end

    function [] = colorscale_call(varargin)
        % Callback for colorscale limits
        S=guidata(gcbf);
        S.rdgnum=S.profileslist(S.prof.Value,:);
        nR=str2num(S.rdgnum); % radargram-number
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1 && S.flag_cs~=1
            S.rb2.Value=0;
            S.rb3.Value=0;
            set(S.ax_data,'ClimMode','auto');
            drawnow;
            S.flag_cs=1;
        elseif val2==1 && S.flag_cs~=2
            S.rb1.Value=0;
            S.rb3.Value=0;
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin1R(nR) S.cmax1R(nR)]);
            drawnow;
            S.flag_cs=2;
        elseif val3==1 && S.flag_cs~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin3R(nR) S.cmax3R(nR)]);
            drawnow;
            S.flag_cs=3;
        end
        guidata(gcbf,S); % Update (for flag_cs)
    end

    function [] = aspect_call(varargin)
        % if aspect ratio is changed
        S=guidata(gcbf);
        if ~strcmp(S.aspect.String{S.aspect.Value},'AspectRatio')
            S.ar=1/str2num(S.aspect.String{S.aspect.Value});
            set(S.ax_data,'DataAspectRatio',[1 S.ar 1]);   
        end
        guidata(gcbf,S); % Update
    end

    function [] = zoom_call(varargin)
        % if zoom of tsl plot is changed
        S=guidata(gcbf);
        temp=extractBefore(S.zoom.String{S.zoom.Value},' %'); % e.g. '100 %'
        S.factor=1/100*str2num(temp); % scaling factor for panelTa
        % update panel Ta:
        val2 = get(S.sliderTH,'Value');
        val3 = get(S.sliderTV,'Value');
        set(S.panelTa,'Position',[-val2*S.factor -val3*S.factor 2 2]*S.factor,'Units','normalized');
        guidata(gcbf,S); % Update
    end


    function [] = newline(varargin)
        % plot new line of current radargram in tsl
        S=guidata(gcbf);
        obj=get(S.ax_tsl,'Children');
        nobj=length(obj);
        for o=1:nobj-length(S.data)
            delete(obj(o)); % delete old lines
        end
        % plot new line
        S.line=plot(S.ax_tsl,S.xprof{str2num(S.rdgnum)},S.yprof{str2num(S.rdgnum)},'r','Linewidth',2);
        drawnow;
        guidata(gcbf,S); % Update
    end
end