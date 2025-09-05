function []=PondView_GUI(xgrid,ygrid,list,mask_interp,t,tz_flag,maxElevation,tr,radar,coords,pfad,coordtrans)

% function [] = PondView_GUI(xgrid,ygrid,list,mask_interp,t,tz_flag,maxElevation,tr,radar,coords,pfad,coordtrans)
%
% Plot thickslices in PondView (Grasmuek & Viggiano 2018)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% xgrid, ygrid: Coordinate grids in m
% list: list of slice-files
% mask_interp: mask with 1 where data is
% t: vector with of time/depth samples
% tz_flag: =2 if depthslices instead of timeslices (tz_flag=1)
% maxElevation: if tz_flag=2: Elevation (a.s.l) at depth=0m (maxElevation=[] for tz_flag=1)
% tr: time/depth vector for radargrams
% radar: radargrams in cells
% coords: gloab_coords in cells for radargrams
% pfad: path to slices folder
%
% requires folder Subfunctions

S.pfad=pfad;

% set data of slices
S.x=xgrid; % in global coordinates
S.y=ygrid;
S.maxElevation=maxElevation;
S.timedepth=tz_flag; % time/depth flag
if S.timedepth==1
    S.unit='ns';
else
    S.unit='m';
end
[S.r,S.c]=size(S.x);

S.t = t; % time vector for slices
S.dx=xgrid(1,2)-xgrid(1,1);
S.dy=ygrid(2,1)-ygrid(1,1);
S.mask=mask_interp;

% initialize values for slices
S.tt=0; % time/depth of starttime of slice (in ns or m)
if S.timedepth==1
    S.dz=10; % thickness of slice in ns
else
    S.dz=0.2; % thickness in m
end
S.thresh=0.2; % Amplitude threshold as min

s=S.tt:S.dz/2:S.t(end)-S.dz/2;
e=S.dz:S.dz/2:S.t(end);
S.t_tsl(:,1)=s(1:min([length(s),length(e)])'); % start of thickslice
S.t_tsl(:,2)=e(1:min([length(s),length(e)])'); % end of thickslice

S.tslnum=1; % current slice number

% load slices:
disp('Loading slices...')
S.temp=load(fullfile(list(1).folder,list(1).name));
[r,c]=size(S.temp.slice); % size of slice 1
S.n=length(list);
S.slices=zeros(r,c,S.n);
for i=1:S.n
    temp=extractBetween(list(i).name,'_','.');
    num=str2double(temp{1}); % number of slice from name
    sl=load(fullfile(list(i).folder,list(i).name));
    S.slices(:,:,num)=sl.slice; % set this slice at correct position
end
rmfield(S,'temp');



% set radargramms (=data)
S.data=radar;
S.ar=1/10; % aspect ratio factor
S.tr=tr; % time vector for radargrams
if S.timedepth==2
    S.tr=abs(S.tr-S.tr(1)); % depth is positive down, starting at 0m -> corresponds to dsl
end

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
S.profileslist=int2str([1:length(coords)]'); % List of profile names

S.rdgnum='1'; % current radargram number

S.factor=1; % zoom factor for tsl plot

S.coordtrans=coordtrans;

% apply coordtrans to radargram (global->lokal)
for i=1:length(S.xprof)
    temp=helmert([S.xprof{i} S.yprof{i}],S.coordtrans(:,3:4),S.coordtrans(:,1:2));
    S.xprof{i}=temp(:,1);
    S.yprof{i}=temp(:,2);
end


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
            'name','PondView',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui,'Renderer','painters','KeyPressFcn',@keypress);
        
        S.fh.MenuBar = 'figure';
        set (S.fh, 'WindowButtonMotionFcn', @mouseMove); % for connecting of Tsl and radargrams
        
        
        %%% Thickslices
        % make first thick slice
        S.ind=find(S.t>=S.tt & S.t<=S.tt+S.dz); % indices of samples in this thickslice
        S.nn=length(S.ind)-mod(length(S.ind),2); % even number of indices
        S.ind=S.ind(1:S.nn); % even number of indices
        S.maxampl=max(abs(S.slices(:,:,S.ind(1):S.ind(end)))); % max abs(amplitude) in this thickslice
        S.ampl=abs(S.slices(:,:,S.ind(1):S.ind(end))./max(S.maxampl,[],2)); % normalize amplitudes in this thickslice
        % apply amplitude threshold:
        S.ampl(S.ampl<S.thresh)=0;

        % create colorscale
        S.nbw=256; % number of black-white values
        S.cmap=[[linspace(1,0,S.nn/2)'; zeros(S.nn/2,1)] [ones(S.nn/2,1); linspace(1,0,S.nn/2)'] [zeros(S.nn/2,1); linspace(0,1,S.nn/2)']];
        S.scale=linspace(0,1,S.nbw); % scaling factor for black-white
        cmap_all=[];
        for i=1:S.nbw
            cmap_all=[cmap_all; S.cmap.*S.scale(i)];
        end
        S.imcb=reshape(linspace(0,1,length(cmap_all(:,1))),[S.nn S.nbw]); % image of colorscale

        % fuse different sample slices with transparency, see here for
        % formula: https://graphicdesign.stackexchange.com/questions/93450/overlap-rgb-values-to-produce-a-single-color
        S.mask2=repmat(~S.mask,[1 1 3]); % RGB picture from mask (=background)
        S.imm=(ones(size(S.mask2))-repmat(S.ampl(:,:,S.nn),[1 1 3])).*S.mask2;
        S.imm=S.imm + repmat(S.ampl(:,:,S.nn),[1 1 3]).*repmat(permute(S.cmap(S.nn,:),[3 1 2]),[S.r,S.c,1]);
        for i=S.nn-1:-1:1
            S.imm = (ones(size(S.mask2))-repmat(S.ampl(:,:,i),[1 1 3])).*S.imm;
            S.imm=S.imm + repmat(S.ampl(:,:,i),[1 1 3]).*repmat(permute(S.cmap(i,:),[3 1 2]),[S.r,S.c,1]);
        end
        S.imm(isnan(S.imm))=1;


        % scrollbars for Tsl:
        S.panelT=uipanel('Parent',S.fh);
        S.panelTa=uipanel('Parent',S.panelT);
        if length(S.x(:,1))>length(S.x(1,:)) % Tsl height>width
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
            'OuterPosition',[0.02 0.02 0.98 0.98],'Parent',S.panelTa,'NextPlot','add');
        % make initial plot TSL
        S.tsl=imagesc(S.ax_tsl,S.x(1,:),S.y(:,1),S.imm);
        S.line=plot(S.ax_tsl,S.xprof{1},S.yprof{1},'r','Linewidth',2); % plot radargram location
        clim([0 S.nn])
        colormap(S.ax_tsl,S.cmap)
        xlabel('x [m]')
        ylabel('y [m]')
        axis xy
        set(S.ax_tsl,'DataAspectratio',[1 1 1])
        cb=colorbar;
        pos=get(cb,'Position');
        pos(3)=pos(3)*2;
        pos(4)=pos(4)/2;
        set(cb,'Visible','off') % do not plot the standard colorbar
        S.axcb=axes('Position',pos,'Color','k'); % instead make new axes for own colorbar
        S.cb=imagesc(S.axcb,S.imcb','AlphaData',repmat(linspace(0,1,S.nn),[S.nbw 1]));
        set(S.axcb,'Color','k','YAxisLocation','right','XTick',[1 S.nn],'XTickLabel',[0 1],'YTick',[1 S.nbw/2 S.nbw],'YTickLabel',[S.tt S.tt+S.dz/2 S.tt+S.dz])
        xlabel(S.axcb,'normalized amplitude')
        if S.timedepth==1
            ylabel(S.axcb,'Time [ns]')
        else
            ylabel(S.axcb,'Depth [m]')
        end
        colormap(S.axcb,S.cmap)

        %%% RADARGRAM
        % scrollbar for radargram:
        S.panelR=uipanel('Parent',S.fh);
        S.panelRa=uipanel('Parent',S.panelR);
        if length(S.x(:,1))>length(S.x(1,:)) % Tsl height>width
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
        S.h_rdg=imagesc(S.ax_data,S.loc_coord{1},S.tr,S.data{1});
        colormap(S.ax_data,flipud(gray));
        xlabel('ProfileX [m]')
        set(S.ax_data,'Xlim',[min(S.loc_coord{1}) max(S.loc_coord{1})],'CLim',[S.cmin3R(1) S.cmax3R(1)])
        S.ax_data.NextPlot='add';
        % plot tsl-range
        S.range1=plot(S.ax_data,[min(S.loc_coord{1}) max(S.loc_coord{1})],[S.t(S.ind(1)) S.t(S.ind(1))],'r','Linewidth',1);
        S.range2=plot(S.ax_data,[min(S.loc_coord{1}) max(S.loc_coord{1})],[S.t(S.ind(end)) S.t(S.ind(end))],'r','Linewidth',1);
        if S.timedepth==1
            ylabel('t [ns]')
            axis ij
            set(S.ax_data,'DataAspectRatio',[1 1/S.ar 1]);
        elseif S.timedepth==2 % horizontal depthslices
            ylabel('z [m]')
            axis ij
            set(S.ax_data,'DataAspectRatio',[1 1 1]);
        end
        grid on
        axis tight

        
        % UI control for tsl-time text
        if S.timedepth==1
            S.tsltext=uicontrol(S.fh,'style','text','unit','pix','position',[15 580 100 15],'String',['Tsl: ',num2str(S.t(S.ind(1)),'%.1f'),' - ',num2str(S.t(S.ind(end)),'%.1f'),' ns'],'HorizontalAlignment','left','callback',@tsltextcall,'FontWeight','bold');
        else
            S.tsltext=uicontrol(S.fh,'style','text','unit','pix','position',[15 580 100 15],'String',['Tsl: ',num2str(S.t(S.ind(1)),'%.1f'),' - ',num2str(S.t(S.ind(end)),'%.1f'),' m'],'HorizontalAlignment','left','callback',@tsltextcall,'FontWeight','bold');
        end
        
        % UI control for aspect ratio
        S.aspect=uicontrol(S.fh,'Style','popupmenu','unit','pix','position',[5 195 110 25],'Value',12,'String',[{'AspectRatio'} {'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'} {'1/20'}],'callback',{@aspect_call});
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[10 225 100 15],'String','Aspect ratio','HorizontalAlignment','left','FontWeight','bold');
        if S.timedepth==2
            S.aspect.Value=7;
        end

        % UI control for tsl zoom
        S.zoom=uicontrol(S.fh,'Style','popupmenu','unit','pix','position',[5 550 110 25],'Value',5,'String',[{'10 %'} {'20 %'} {'50 %'} {'75 %'} {'100 %'} {'150 %'} {'175 %'} {'200 %'} {'300 %'}],'callback',{@zoom_call});
        

        % UI control for profile number 
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[10 80 100 85],'String',S.profileslist,'Value',1,'callback',{@profnum_call});
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[10 170 100 15],'String','Choose profile','HorizontalAlignment','left','FontWeight','bold');

        % UI control for thickness of slices 
        S.thick=uicontrol(S.fh,'Style','edit','unit','pix','position',[10 420 100 15],'String',num2str(S.dz,'%.1f'),'callback',{@thick_call});
        S.thicktext=uicontrol(S.fh,'style','text','unit','pix','position',[10 440 100 30],'String',['Thickness of slice [',S.unit,']'],'HorizontalAlignment','left','FontWeight','bold');

        % UI control for threshold 
        S.threshold=uicontrol(S.fh,'Style','edit','unit','pix','position',[10 360 100 15],'String',num2str(S.thresh,'%.2f'),'callback',{@threshold_call});
        S.threshtext=uicontrol(S.fh,'style','text','unit','pix','position',[10 380 100 30],'String','Threshold for slices (0-1)','HorizontalAlignment','left','FontWeight','bold');

        % UI control for start
        S.start=uicontrol(S.fh,'Style','edit','unit','pix','position',[10 480 100 15],'String',num2str(S.tt,'%.1f'),'callback',{@start_call});
        S.starttext=uicontrol(S.fh,'style','text','unit','pix','position',[10 500 100 15],'String',['Start of slices [',S.unit,']'],'HorizontalAlignment','left','FontWeight','bold');

        
        % UI control - radiobuttons for layout
        S.lay1=uicontrol('Style','radiobutton','String','Layout 1','Position',[10 35 100 15],'Callback',@layout_call); 
        S.lay2=uicontrol('Style','radiobutton','String','Layout 2','Position',[10 15 100 15],'Callback',@layout_call);
        if S.layout==1
            S.lay1.Value=1;
        else
            S.lay2.Value=1;
        end
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[10 365-50 100 15],'Callback',@colorscale_call);
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[10 345-50 100 15],'Callback',@colorscale_call);
        S.rb3=uicontrol('Style','radiobutton','String','3 %','Position',[10 325-50 100 15],'Callback',@colorscale_call,'Value',1);  % this button is on
        S.flag_cs=3;
               
        % UI control - button for saving of georefrenced images
        S.saveall=uicontrol('Style','pushbutton','String','Save georeferenced slices','Position',[150 10 180 30],'Value',0,'Callback',@saveall_call);

        % UI control for depth scale switch
        if S.timedepth==2
            S.dswitch=uicontrol('Style','checkbox','String','Absolute depth','Position',[400 10 180 30],'Value',0,'Callback',@dswitch_call);
        end

        % UI control for info text
        S.infotext=uicontrol('Style','text','String','For scrolling: move mouse inside plot and use arrow keys','Position',[600 10 200 30]);
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
        set(S.axcb,'OuterPosition',[0.9 0.3 0.1 0.4],'Units','normalized','Parent',S.panelTa);

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


    function mouseMove(object, eventdata)
        % get position inside radargram
        C = get (S.ax_data, 'CurrentPoint');
        S=guidata(object);
        
        % if mouse inside radargram:
        if C(1,1)>=S.ax_data.XLim(1) && C(1,1)<=S.ax_data.XLim(2) && C(1,2)>=S.ax_data.YLim(1) && C(1,2)<=S.ax_data.YLim(2)
            set(gcf,'pointer','crosshair');

            % find global coords of point:
            S.xy=[interp1(S.loc_coord{str2num(S.rdgnum)},S.xprof{str2num(S.rdgnum)},C(1,1)) interp1(S.loc_coord{str2num(S.rdgnum)},S.yprof{str2num(S.rdgnum)},C(1,1))];

            % find tsl-number:
            n=find(C(1,2)>S.t_tsl(:,1) & C(1,2)<=S.t_tsl(:,2),1,'first');

            if n~=S.tslnum % if new tsl-number
                S.tslnum=n; % set new tsl-number
                % update tsl range in radargram
                if isfield(S,'range1') && isvalid(S.range1) % is handle is there and valid -> update
                    S.range1.YData=[S.t_tsl(n,1) S.t_tsl(n,1)];
                    S.range2.YData=[S.t_tsl(n,2) S.t_tsl(n,2)];
                else % plot new
                    if length(S.ax_data.Children)>3
                        obj=S.ax_data.Children;
                        for i=1:length(obj)-1
                            delete(obj(i));
                        end
                    end
                    % plot tsl-range
                    S.range1=plot(S.ax_data,[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})],[S.t_tsl(n,1) S.t_tsl(n,1)],'r','Linewidth',1);
                    S.range2=plot(S.ax_data,[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})],[S.t_tsl(n,2) S.t_tsl(n,2)],'r','Linewidth',1);
                end
                guidata(object,S); % Update
                tsltextcall(); % change display of Tsl-times and plot new slice
            end
            % update cross in tsl-plot:
            if abs(mean(S.x(1,:))-mean(S.xprof{str2num(S.rdgnum)}))<=1000 && abs(mean(S.y(:,1))-mean(S.yprof{str2num(S.rdgnum)}))<=1000
                obj=get(S.ax_tsl,'Children');
                nobj=length(obj);
                for o=1:nobj-2
                    delete(obj(o)); % delete old lines
                end
                % plot current position in slice
                plot(S.ax_tsl,S.xy(1),S.xy(2),'r*','Linewidth',4,'Markersize',8) % red thick dot for current position in slice
            end
            
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

    function [] = start_call(varargin)
        % new start time/depth for slices
        S=guidata(gcbf);
        S.tt=str2double(S.start.String);
        s=S.tt:S.dz/2:S.t(end)-S.dz;
        e=S.tt+S.dz:S.dz/2:S.t(end);
        S=rmfield(S,'t_tsl');
        S.t_tsl(:,1)=s(1:min([length(s),length(e)])'); % start of thickslice
        S.t_tsl(:,2)=e(1:min([length(s),length(e)])'); % end of thickslice
        S.tslnum=1;
        guidata(gcbf,S); % Update
        tsltextcall();
        profnum_call();
    end

    function [] = threshold_call(varargin)
        % new threshold for slices
        S=guidata(gcbf);
        S.thresh=str2double(S.threshold.String);
        guidata(gcbf,S); % Update
        tslplot();
    end


    function [] = thick_call(varargin)
        % new thickness for slices
        S=guidata(gcbf);
        S.dz=str2double(S.thick.String);
        s=S.tt:S.dz/2:S.t(end)-S.dz;
        e=S.tt+S.dz:S.dz/2:S.t(end);
        S=rmfield(S,'t_tsl');
        S.t_tsl(:,1)=s(1:min([length(s),length(e)])'); % start of thickslice
        S.t_tsl(:,2)=e(1:min([length(s),length(e)])'); % end of thickslice
        S.tslnum=1;
        guidata(gcbf,S); % Update
        tsltextcall();
        profnum_call();
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
        for o=1:nobj
            delete(obj(o)); % delete old lines
        end
        S.ax_data.NextPlot='replace';
        S.h_rdg=imagesc(S.ax_data,S.loc_coord{str2num(S.rdgnum)},S.t,S.data{str2num(S.rdgnum)});
        colormap(S.ax_data,flipud(gray));
        S.ax_data.NextPlot='add';
        % plot tsl-range
        S.range1=plot(S.ax_data,[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})],[S.t_tsl(S.tslnum,1) S.t_tsl(S.tslnum,1)],'r','Linewidth',1);
        S.range2=plot(S.ax_data,[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})],[S.t_tsl(S.tslnum,2) S.t_tsl(S.tslnum,2)],'r','Linewidth',1);
        set(S.ax_data,'Xlim',[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})]);
        if S.timedepth==1
            ylabel(S.ax_data,'t [ns]')
        else
            ylabel(S.ax_data,'z [m]')
        end
        xlabel(S.ax_data,'ProfileX [m]')
        axes(S.ax_data);
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
            S.t_tsl=S.maxElevation-S.t_tsl;
            S.t_tsl=fliplr(S.t_tsl); % to make first row smaller than second
            axes(S.ax_data);
            axis xy
        else
            % relative depth from top of radargram
            S.t=abs(S.t-S.maxElevation);
            S.t_tsl=abs(S.t_tsl-S.maxElevation);
            S.t_tsl=fliplr(S.t_tsl); % to make first row smaller than second
            axes(S.ax_data);
            axis ij
        end
        guidata(gcbf,S); % Update
        profnum_call();
    end


    function [] = tsltextcall(varargin)
        % callback for tsl-textbox (display tsl-times) and change tsl
        % display
        S=guidata(gcbf);
        n=S.tslnum;
        if S.timedepth==1
            S.tsltext.String=['Tsl: ',num2str(S.t_tsl(n,1),'%.2f'),' - ',num2str(S.t_tsl(n,2),'%.2f'),' ns'];
        else
            S.tsltext.String=['Dsl: ',num2str(S.t_tsl(n,1),'%.2f'),' - ',num2str(S.t_tsl(n,2),'%.2f'),' m'];
        end
        guidata(gcbf,S); % Update
        % call tsl-plotting function
        tslplot(varargin);
    end

    function [] = tslplot(varargin)
        % function for plotting new thickslice
        S=guidata(gcbf);
        n=S.tslnum; % tsl-number

        S.ind=find(S.t_tsl(n,1)<=S.t & S.t_tsl(n,2)>=S.t); % indices of samples in this thickslice
        S.nn=length(S.ind)-mod(length(S.ind),2); % even number of indices
        S.ind=S.ind(1:S.nn); % even number of indices
        S.maxampl=max(abs(S.slices(:,:,S.ind(1):S.ind(end)))); % max abs(amplitude) in this thickslice
        S.ampl=abs(S.slices(:,:,S.ind(1):S.ind(end))./max(S.maxampl,[],2)); % normalize amplitudes in this thickslice
        % apply amplitude threshold:
        S.ampl(S.ampl<S.thresh)=0;

        % make new colorscale-image
        S.cmap=[[linspace(1,0,S.nn/2)'; zeros(S.nn/2,1)] [ones(S.nn/2,1); linspace(1,0,S.nn/2)'] [zeros(S.nn/2,1); linspace(0,1,S.nn/2)']];
        S.scale=linspace(0,1,S.nbw); % scaling factor for black-white
        cmap_all=[];
        for i=1:S.nbw
            cmap_all=[cmap_all; S.cmap.*S.scale(i)];
        end
        S.imcb=reshape(linspace(0,1,length(cmap_all(:,1))),[S.nn S.nbw]); % image of colorscale

        % fuse images
        S.imm=(ones(size(S.mask2))-repmat(S.ampl(:,:,S.nn),[1 1 3])).*S.mask2;
        S.imm=S.imm + repmat(S.ampl(:,:,S.nn),[1 1 3]).*repmat(permute(S.cmap(S.nn,:),[3 1 2]),[S.r,S.c,1]);
        for i=S.nn-1:-1:1
            S.imm = (ones(size(S.mask2))-repmat(S.ampl(:,:,i),[1 1 3])).*S.imm;
            S.imm=S.imm + repmat(S.ampl(:,:,i),[1 1 3]).*repmat(permute(S.cmap(i,:),[3 1 2]),[S.r,S.c,1]);
        end
        S.imm(isnan(S.imm))=1;

        % plot new fused image
        set(S.tsl,'CData',S.imm);
        clim(S.ax_tsl,[0 S.nn]);
       
        % update colorbar:
        set(S.cb,'CData',S.imcb','AlphaData',repmat(linspace(0,1,S.nn),[S.nbw 1]));
        set(S.axcb,'Color','k','YAxisLocation','right','XLim',[1 S.nn],'XTick',[1 S.nn],'XTickLabel',[0 1],'YTick',[1 S.nbw/2 S.nbw],'YTickLabel',[S.t_tsl(n,1) S.t_tsl(n,1)+S.dz/2 S.t_tsl(n,1)+S.dz]);

        drawnow;
        guidata(gcbf,S); % Update (!)
    end


    function [] = saveall_call(varargin)
        % save all slices as georeferenced png:
        S=guidata(gcbf);
        set(gcf,'pointer','watch')
        
        folder=['georef_start',S.start.String,S.unit,'_thick',S.thick.String,S.unit,'_threshold',S.threshold.String];
        if ~exist(fullfile(S.pfad,folder))
            mkdir(fullfile(S.pfad,folder))
        end
        for n=1:length(S.t_tsl(:,1)) % go through all slices
            %%% Georeferenced png:
            S.ind=find(S.t_tsl(n,1)<=S.t & S.t_tsl(n,2)>=S.t); % indices of samples in this thickslice
            S.nn=length(S.ind)-mod(length(S.ind),2); % even number of indices
            S.ind=S.ind(1:S.nn); % even number of indices
            S.maxampl=max(abs(S.slices(:,:,S.ind(1):S.ind(end)))); % max abs(amplitude) in this thickslice
            S.ampl=abs(S.slices(:,:,S.ind(1):S.ind(end))./max(S.maxampl,[],2)); % normalize amplitudes in this thickslice
            % apply amplitude threshold:
            S.ampl(S.ampl<S.thresh)=0;

            % make new colorscale-image
            S.cmap=[[linspace(1,0,S.nn/2)'; zeros(S.nn/2,1)] [ones(S.nn/2,1); linspace(1,0,S.nn/2)'] [zeros(S.nn/2,1); linspace(0,1,S.nn/2)']];
            S.scale=linspace(0,1,S.nbw); % scaling factor for black-white
%             cmap_all=[];
%             for i=1:S.nbw
%                 cmap_all=[cmap_all; S.cmap.*S.scale(i)];
%             end
%             S.imcb=reshape(linspace(0,1,length(cmap_all(:,1))),[S.nn S.nbw]); % image of colorscale
%             
            for i=1:S.nn
                S.cim(:,:,i)=zeros(S.nn,S.nn/2);
                S.cim(i,:,i)=linspace(0,1,S.nn/2);
            end
            S.cimm=zeros(S.nn,S.nn/2,3)-repmat(S.cim(:,:,S.nn),[1 1 3]);
            S.cimm=S.cimm+repmat(S.cim(:,:,S.nn),[1 1 3]).*repmat(permute(S.cmap(S.nn,:),[3 1 2]),[S.nn,S.nn/2,1]);
            for i=S.nn-1:-1:1
                S.cimm = (ones(S.nn,S.nn/2,3)-repmat(S.cim(:,:,i),[1 1 3])).*S.cimm;
                S.cimm=S.cimm + repmat(S.cim(:,:,i),[1 1 3]).*repmat(permute(S.cmap(i,:),[3 1 2]),[S.nn,S.nn/2,1]);
            end
            S.cimm(isnan(S.cimm))=1;
            imwrite(S.cimm,S.cmap,fullfile(S.pfad,folder,'Colorscale.png'));


            % fuse images
            S.imm=(ones(size(S.mask2))-repmat(S.ampl(:,:,S.nn),[1 1 3])).*S.mask2;
            S.imm=S.imm + repmat(S.ampl(:,:,S.nn),[1 1 3]).*repmat(permute(S.cmap(S.nn,:),[3 1 2]),[S.r,S.c,1]);
            for i=S.nn-1:-1:1
                S.imm = (ones(size(S.mask2))-repmat(S.ampl(:,:,i),[1 1 3])).*S.imm;
                S.imm=S.imm + repmat(S.ampl(:,:,i),[1 1 3]).*repmat(permute(S.cmap(i,:),[3 1 2]),[S.r,S.c,1]);
            end
            S.imm(isnan(S.imm))=1;

            tslname = fullfile(S.pfad,folder,['Slice_',num2str(S.t_tsl(n,1),'%.1f'),'-',num2str(S.t_tsl(n,2),'%.1f'),S.unit,'.png']);              
            imwrite(flipud(S.imm),S.cmap,tslname);
            

            % write pngw
            fname = ['Slice_',num2str(S.t_tsl(n,1),'%.1f'),'-',num2str(S.t_tsl(n,2),'%.1f'),S.unit,'.pgw'];
            write_geoPNGW(S.x,S.y,S.coordtrans,fullfile(S.pfad,folder,fname));
        end
        disp('All slices saved!')
        set(gcf,'pointer','arrow')
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
        for o=1:nobj-1
            delete(obj(o)); % delete old lines
        end
        if abs(mean(S.x(1,:))-mean(S.xprof{str2num(S.rdgnum)}))<=1000 && abs(mean(S.y(:,1))-mean(S.yprof{str2num(S.rdgnum)}))<=1000
            % plot new line
            S.line=plot(S.ax_tsl,S.xprof{str2num(S.rdgnum)},S.yprof{str2num(S.rdgnum)},'r','Linewidth',2);
            drawnow;
        end
        guidata(gcbf,S); % Update
    end
end