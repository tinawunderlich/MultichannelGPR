function [] = plot_layerpicking(radargrams,global_coords,x,tz,tzflag,folder)

% function [] = plot_layerpicking(radargrams,global_coords,x,tz,tzflag,folder)
%
% Plot radargrams and picking of layers
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% radargrams: radargrams in cells
% global_coords: global coordinates of radargrams (cell array of matrices)
% x: cell array of local profile coordinates of radargrams
% tz: vector of time or depth
% tzflag: =1 for time and =2 for depth
% folder: path to radargrams and other inputs (will be used for saving of
% picks)


% set data
S.x = x;
S.global_coords = global_coords;
if tzflag==1
    S.t = tz;
    S.z=NaN;
else
    S.t=NaN;
    S.z=tz;
end
S.dtz=tz(2)-tz(1);
S.radargrams=radargrams;
S.folder=folder;
S.col=[];
S.profilelist=num2cell([1:length(x)]');
S.allpicks=[];
S.prof_old=1;

% get screensize
S.siz = get( 0, 'Screensize' );

% Add the UI components
S=addcomponents(S);

% Make figure visible after adding components
S.fh.Visible='on';

    function S=addcomponents(S)
        
        %---------------------------------------
        % create figure handle:
        S.fh = figure('units','pixels',...
            'position',S.siz,...
            'menubar','none',...
            'name','Layer picking',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui,'Renderer','painters','KeyPressFcn',@keypress);
        
        % --
        % upper plot for comparison
        % scrollbar:
        S.panel1=uipanel('Parent',S.fh);
        S.panel1a=uipanel('Parent',S.panel1);
        set(S.panel1,'Position',[210 80+S.siz(4)/2+80 S.siz(3)-220 S.siz(4)/2-100],'Units','pixels'); % panel for cutting frame with slider
        set(S.panel1a,'Position',[210 80+S.siz(4)/2+80 2000 S.siz(4)/2-100],'Units','pixels'); % panel inside cutting frame with image

        S.slider1=uicontrol('Style','Slider','Parent',S.fh,'Units','pixels','Position',[210 80+S.siz(4)/2+80 S.siz(3)-220 0.05],'Value',0.5,'Callback',@slider_callback1);
        
        S.ax1 = axes('OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panel1a);
        % make initial plot
        if isnan(S.z)
            S.Plot1 = imagesc(S.ax1,S.x{1},S.t,S.radargrams{1});
        else
            S.Plot1 = imagesc(S.ax1,S.x{1},S.z,S.radargrams{1});
        end
        colormap(flipud(gray));
        xlabel('x [m]')
        if isnan(S.t)
            ylabel('z [m]')
            axis xy 
        else
            ylabel('t [ns]')
            axis ij
        end
        set(S.ax1,'DataAspectratio',[0.5 1 1])
        
        S.title1 = uicontrol(S.fh,'style','text',...
            'unit','pix',...
            'position',[210 900 200 25],'String','For comparison: Profile #','FontWeight','bold');
        

        %--
        % Lower plot for picking:
        % scrollbar:
        S.panel2=uipanel('Parent',S.fh);
        S.panel2a=uipanel('Parent',S.panel2);
        set(S.panel2,'Position',[210 120 S.siz(3)-220 S.siz(4)/2-100],'Units','pixels');
        set(S.panel2a,'Position',[210 120 2000 S.siz(4)/2-100],'Units','pixels');

        S.slider2=uicontrol('Style','Slider','Parent',S.fh,'Units','pixels','Position',[210 120 S.siz(3)-220 0.05],'Value',0.5,'Callback',@slider_callback2);
        
        S.ax2 = axes('OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panel2a);
        % make initial plot
        if isnan(S.z)
            S.Plot2 = imagesc(S.ax2,S.x{1},S.t,S.radargrams{1});
        else
            S.Plot2 = imagesc(S.ax2,S.x{1},S.z,S.radargrams{1});
        end
        colormap(flipud(gray));
        xlabel('x [m]')
        if isnan(S.t)
            ylabel('z [m]')
            axis xy 
        else
            ylabel('t [ns]')
            axis ij
        end
       set(S.ax2,'DataAspectratio',[0.5 1 1],'nextplot','add');

       S.title2 = uicontrol(S.fh,'style','text',...
            'unit','pix',...
            'position',[210 450 200 25],'String',['Active Picking: Profile #',int2str(1)],'FontWeight','bold');

        
        % UI control for comparison profile number
        S.comparenum=uicontrol(S.fh,'style','edit','String','1','Position',[400 900 30 25],'Callback',@compare_call);
        S.auto=uicontrol(S.fh,'style','checkbox','Value',1,'String','auto','Position',[450 900 100 25]);

        % UI control for connect sliders button
        S.connect=uicontrol(S.fh,'Style','checkbox','Value',1,'String','connect scrollbars','Position',[550 900 150 20]);
        
        % UI control aspectratio
        S.asp_String=[{'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'} {'1/20'}];
        S.aspval=1/2;
        S.asp = uicontrol(S.fh,'style','listbox','unit','pix','position',[10 120 80 60],'String',S.asp_String,'Value',7,'Callback',@asp_call);
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[10 180 80 20],'String','Aspect ratio','HorizontalAlignment','left');
        
        % Create switch button
        S.Switchbutton = uicontrol(S.fh, 'Style','listbox');
        S.Switchbutton.Position = [480 20 70 30];
        S.Switchbutton.String = {'Line','Points'};
        S.Switchbutton.Value = 1;

        
        % UI control - radiobuttons for colorscale
        S.rbtext=uicontrol(S.fh,'Style','text','string','Color scale','Position',[10 80 80 20],'HorizontalAlignment','left');
        S.rb1=uicontrol(S.fh,'Style','radiobutton','String','Auto','Position',[10 60 80 20],'Value',1,'Callback',@rb_call); % this button is on
        S.rb2=uicontrol(S.fh,'Style','radiobutton','String','1 %','Position',[10 40 80 20],'Callback',@rb_call);
        S.rb3=uicontrol(S.fh,'Style','radiobutton','String','3 %','Position',[10 20 80 20],'Callback',@rb_call);
        S.flag=1;
        
        % UI control - pushbutton for saving
        S.save=uicontrol(S.fh,'Style','pushbutton','String','Save all picks','Position',[10 220 180 25],'Callback',@save_picks);
        
        % UI control - pushbutton for loading
        S.load=uicontrol(S.fh,'Style','pushbutton','String','Load picks','Position',[10 250 180 25],'Callback',@load_picks);
        
        % UI control - pushbutton for deleting of layer
        S.dellayer=uicontrol(S.fh,'Style','pushbutton','String','Delete current layer','Position',[10 280 180 25],'Callback',@delete_layer);
        
        
        % UI control for profile number
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[100 20 90 160],'String',S.profilelist,'Value',1,'callback',@profnum_call);
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[100 180 90 20],'String','Profile','HorizontalAlignment','left');
        
        % UI control for layer list
        S.layerlist=uicontrol(S.fh,'Style','listbox','unit','pix','position',[10 320 180 160],'String','','Value',1,'callback',@layerlist_call);
        S.layertext=uicontrol(S.fh,'style','text','unit','pix','position',[10 480 90 20],'String','List of layers','HorizontalAlignment','left');
        
        
        % UI control - pushbutton for delete all picks
        S.del_all=uicontrol(S.fh,'Style','pushbutton','String','Delete all picks in profile','Position',[700 20 180 25],'Callback',@delall_call);
                
        % UI control for info text
        S.infotext = uicontrol(S.fh,'style','text',...
            'unit','pix',...
            'position',[700 50 200 25],'String','For deleting click on one line.');
        
        % UI control for layer ID & name
        S.idtext=uicontrol(S.fh,'style','text','String','Layer ID:','Position',[200 45 50 25],'HorizontalAlignment','Left');
        S.id=uicontrol(S.fh,'style','edit','String','1','Position',[250 50 40 25]);
        S.idtext2=uicontrol(S.fh,'style','text','String','Layer name:','Position',[300 45 80 25],'HorizontalAlignment','Left');
        S.name=uicontrol(S.fh,'style','edit','String','Layer','Position',[370 50 80 25],'Callback',@name_call);
                       
        % UI control - pushbutton for new layer
        S.new_layer=uicontrol(S.fh,'Style','pushbutton','String','Start picking','Position',[200 20 250 25],'Callback',@newlayer_call);
        
        % UI control for info text
        S.infotext=uicontrol('Style','text','String','For scrolling: move mouse inside plot and use arrow keys','Position',[10 550 180 50]);
    end

% save handles
guidata(S.fh,S);



%%% Callback functions:

    function keypress(hObject,event)
        pos=event.Source.CurrentPoint; % current point in figure
        pos_ax1_panel=S.panel1.Position; % position of upper radargram-panel
        pos_ax2_panel=S.panel2.Position; % position of lower radargram-panel
       
        % find out if mouse is inside a plot
        if pos(1)>=pos_ax1_panel(1) && pos(1)<=pos_ax1_panel(1)+pos_ax1_panel(3) && pos(2)>=pos_ax1_panel(2) && pos(2)<=pos_ax1_panel(2)+pos_ax1_panel(4)
            % Mouse in upper radargram-panel
            val = get(S.slider1,'Value');
            if strcmp(event.Key,'rightarrow')
                set(S.panel1a,'Position',[-(val+0.05) 0 2 1],'Units','normalized')
                S.slider1.Value=val+0.05;
                if S.connect.Value==1
                    set(S.panel2a,'Position',[-(val+0.05) 0 2 1],'Units','normalized')
                    set(S.slider2,'Value',val+0.05);
                end
            elseif strcmp(event.Key,'leftarrow')
                set(S.panel1a,'Position',[-(val-0.05) 0 2 1],'Units','normalized')
                S.slider1.Value=val-0.05;
                if S.connect.Value==1
                    set(S.panel2a,'Position',[-(val-0.05) 0 2 1],'Units','normalized')
                    set(S.slider2,'Value',val-0.05);
                end
            end
        elseif pos(1)>=pos_ax2_panel(1) && pos(1)<=pos_ax2_panel(1)+pos_ax2_panel(3) && pos(2)>=pos_ax2_panel(2) && pos(2)<=pos_ax2_panel(2)+pos_ax2_panel(4)
            % Mouse in lower radargram-panel
            val = get(S.slider2,'Value');
            if strcmp(event.Key,'rightarrow')
                set(S.panel2a,'Position',[-(val+0.05) 0 2 1],'Units','normalized')
                S.slider2.Value=val+0.05;
                if S.connect.Value==1
                    set(S.panel1a,'Position',[-(val+0.05) 0 2 1],'Units','normalized')
                    set(S.slider1,'Value',val+0.05);
                end
            elseif strcmp(event.Key,'leftarrow')
                set(S.panel2a,'Position',[-(val-0.05) 0 2 1],'Units','normalized')
                S.slider2.Value=val-0.05;
                if S.connect.Value==1
                    set(S.panel1a,'Position',[-(val-0.05) 0 2 1],'Units','normalized')
                    set(S.slider1,'Value',val-0.05);
                end
            end
        end
    end

    function slider_callback1(varargin)
        S=guidata(gcbf);
        val = get(S.slider1,'Value');
        set(S.panel1a,'Position',[-val 0 2 1],'Units','normalized')
        if S.connect.Value==1
            set(S.panel2a,'Position',[-val 0 2 1],'Units','normalized')
            set(S.slider2,'Value',val);
        end
        guidata(gcbf,S); % Update
    end

    function slider_callback2(varargin)
        S=guidata(gcbf);
        val = get(S.slider2,'Value');
        set(S.panel2a,'Position',[-val 0 2 1],'Units','normalized')
        if S.connect.Value==1
            set(S.panel1a,'Position',[-val 0 2 1],'Units','normalized')
            set(S.slider1,'Value',val);
        end
        guidata(gcbf,S); % Update
    end

    function resizeui(hObject,event)

        % get current size of figure
        Position=get(hObject,'Position');
        wid_fig=Position(3); % Figure width
        wid=wid_fig-300; % width of axes
        hei_fig=Position(4); % figure height
        hei=hei_fig-200;    % height of axes
        
        % change size of axes
        set(S.ax1,'OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panel1a) % top
        set(S.ax2,'OuterPosition',[0.02 0.02 0.98 0.98],'Units','normalized','Parent',S.panel2a) % bottom

        % change size of panels (top)
        S.panel1.Position=[210 80+hei/2+80 wid_fig-220 hei_fig/2-100];

        % change size of panels (bottom)
        S.panel2.Position=[210 120 wid_fig-220 hei_fig/2-100];

        % change size of sliders
        S.slider2.Position=[210 120 wid_fig-220 0.05]; % bottom
        S.slider1.Position=[210 80+hei/2+80 wid_fig-220 0.05]; % top

        % change position of auto and comparenum-boxes
        S.comparenum.Position=[380 170+hei 40 20];
        S.auto.Position=[440 170+hei 100 20];
        S.connect.Position=[550 170+hei 200 20];

        % change position of titles:
        S.title1.Position=[210 170+hei 200 20];
        S.title2.Position=[210 120+hei/2 200 20];

        % slider update
        val1 = get(S.slider1,'Value');
        set(S.panel1a,'Position',[-val1 0 2 1],'Units','normalized');
        val2 = get(S.slider2,'Value');
        set(S.panel2a,'Position',[-val2 0 2 1],'Units','normalized');
    end



    function [] = asp_call(varargin)
        % Callback for aspectratio
        S=guidata(gcbf);
        S.aspval=S.asp_String(S.asp.Value);
        S.aspval=str2num(S.aspval{1});
        set(S.ax1,'DataAspectratio',[S.aspval 1 1]);
        drawnow;
        set(S.ax2,'DataAspectratio',[S.aspval 1 1]);
        drawnow;

        S.slider1.Value=0.5;
        S.slider2.Value=0.5;

        guidata(gcbf,S); % Update
    end

    function [] = rb_call(varargin)
        % Callback for radiobuttons for color
        S=guidata(gcbf);
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1 && S.flag~=1
            S.rb2.Value=0;
            S.rb3.Value=0;
            set(S.ax1,'ClimMode','auto');
            drawnow;
            set(S.ax2,'ClimMode','auto');
            drawnow;
            S.flag=1;
        elseif val2==1 && S.flag~=2
            S.rb1.Value=0;
            S.rb3.Value=0;
            % ax1
            coldata=get(S.Plot1,'cdata');
            coldata=sort(unique(coldata(~isnan(coldata(:)))));
            cmin=coldata(round(length(coldata)/100*1));
            cmax=coldata(end-round(length(coldata)/100*1));
            set(S.ax1,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            % ax2
            coldata=get(S.Plot2,'cdata');
            coldata=sort(unique(coldata(~isnan(coldata(:)))));
            cmin=coldata(round(length(coldata)/100*1));
            cmax=coldata(end-round(length(coldata)/100*1));
            set(S.ax2,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            S.flag=2;
        elseif val3==1 && S.flag~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            % ax1
            coldata=get(S.Plot1,'cdata');
            coldata=sort(unique(coldata(~isnan(coldata(:)))));
            cmin=coldata(round(length(coldata)/100*2));
            cmax=coldata(end-round(length(coldata)/100*2));
            set(S.ax1,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            % ax2
            coldata=get(S.Plot2,'cdata');
            coldata=sort(unique(coldata(~isnan(coldata(:)))));
            cmin=coldata(round(length(coldata)/100*2));
            cmax=coldata(end-round(length(coldata)/100*2));
            set(S.ax2,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            S.flag=3;
        end
        guidata(gcbf,S); % Update (for flag)
    end

    function [] = profnum_call(varargin)
        % Callback for new profile
        S=guidata(gcbf);
        % delete old plotted picks
        obj=get(S.ax2,'Children');
        nobj=length(obj);
        for o=1:nobj-1
            delete(obj(o));
        end
        % plot new profile for picking
        set(S.Plot2,'CData',S.radargrams{S.prof.Value},'xdata',S.x{S.prof.Value});
        set(S.ax2,'DataAspectratio',[str2num(S.asp_String{S.asp.Value}) 1 1],'nextplot','add','Xlim',[min(S.x{S.prof.Value}) max(S.x{S.prof.Value})]);
        S.title2.String=['Active Picking: Profile #',int2str(S.prof.Value)];
        rb_call();
        if S.auto.Value==1
            S.comparenum.String=num2str(S.prof_old);
            compare_call();
        end
        S.prof_old=S.prof.Value; % set new old profile number for next change
        % check if picks are available
        if ~isempty(S.allpicks) && any(S.allpicks(:,6)==S.prof.Value)
            set(S.ax2,'DataAspectratio',[str2num(S.asp_String{S.asp.Value}) 1 1],'nextplot','add','Xlim',[min(S.x{S.prof.Value}) max(S.x{S.prof.Value})]);
            picks=S.allpicks(S.allpicks(:,6)==S.prof.Value,:); % picks of this profile
            n=unique(picks(:,5)); % layer IDs
            % get color of layers
            l_list=S.layerlist.String;
            for k=1:length(l_list)
                l_list{k}=l_list{k}(28:end-14); % extract only ID - name
            end
            for j=1:length(n)
                for i=1:length(l_list)
                    % get ID
                    ii=1;
                    while ~isempty(str2num(l_list{i}(1:ii)))
                        ii=ii+1;
                    end
                    id=l_list{i}(1:ii-2);   % id from l_list string
                    in2(i,j)=strcmp(id,num2str(n(j))); % check if layer ID already used
                end
                if any(in2(:,j)==1)
                    w(j)=find(in2(:,j)==1);
                end
            end
            for i=1:length(n)
                % get possible more lines per layer
                numlay=unique(picks(picks(:,5)==n(i),7));
                for j=1:length(numlay)
                    xtemp=picks(picks(:,5)==n(i) & picks(:,7)==numlay(j),1);
                    ytemp=picks(picks(:,5)==n(i) & picks(:,7)==numlay(j),2);
                    flag=mean(picks(picks(:,5)==n(i) & picks(:,7)==numlay(j),8));
                    if flag==1 % line
                        plot(S.ax2,xtemp,ytemp,'Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numlay(j));
                    else % points
                        plot(S.ax2,xtemp,ytemp,'*','Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numlay(j));
                    end
                    text(S.ax2,xtemp(1),ytemp(1)+10*S.dtz,num2str(n(i)),'Color',S.col(w(i),:));
                end
            end
        end
        guidata(gcbf,S); % Update
    end

    function [] = newlayer_call(varargin)
        % Callback for new line of layer
        S=guidata(gcbf);
        % check if line or point:
        if S.Switchbutton.Value==1 % line
            flag=1; % line
        else
            flag=2; % points
        end
        % make new list entry or check for old one
        IDstr=S.id.String;
        name=S.name.String;
        l_list=S.layerlist.String; % list of layers
        for k=1:length(l_list)
            l_list{k}=l_list{k}(28:end-14); % extract only "ID - name (flag)"
        end
        if isempty(l_list) % if first layer
            % Color management
            S.col = [S.col; uisetcolor([1 1 0],'Select a color')];
            laynum=1;
            numlayer=1; % number of line for this layer
            if flag==1
                S.layerlist.String={['<HTML><FONT color="',rgb2Hex(S.col(1,:).*255),'">',[IDstr,' - ',name,' (line)'],'</FONT></HTML>']}; % add first layer to list
            else
                S.layerlist.String={['<HTML><FONT color="',rgb2Hex(S.col(1,:).*255),'">',[IDstr,' - ',name,' (points)'],'</FONT></HTML>']}; % add first layer to list
            end
        else % if there are already layers
            for i=1:length(l_list)
                % get ID
                ii=1;
                while ~isempty(str2num(l_list{i}(1:ii)))
                    ii=ii+1;
                end
                id=l_list{i}(1:ii-2);   % id from l_list string
                in(i)=strcmp(id,IDstr); % check if layer ID already used
            end
            if any(in==1) % this layer ID is already in use...
                w=find(in==1);
                if flag==1 % line
                    if ~strcmp(l_list{w},[IDstr,' - ',name,' (line)']) % if not same name
                        % change name of layer
                        S.layerlist.String{w}=['<HTML><FONT color="',rgb2Hex(S.col(w,:).*255),'">',[IDstr,' - ',name,' (line)'],'</FONT></HTML>'];
                    end
                else
                    if ~strcmp(l_list{w},[IDstr,' - ',name,' (points)']) % if not same name
                        % change name of layer
                        S.layerlist.String{w}=['<HTML><FONT color="',rgb2Hex(S.col(w,:).*255),'">',[IDstr,' - ',name,' (points)'],'</FONT></HTML>'];
                    end
                end
                laynum=w;
                % check if already layer with this ID in this profile
                obj=get(S.ax2,'Children');
                for i=1:length(obj)-1
                    if ~isempty(obj(i).Tag)
                        I(i)=str2num(obj(i).Tag);
                        nn(i)=obj(i).UserData;
                    end
                end
                if exist('I','var') && any(I==str2num(IDstr))
                    numlayer=max(nn(I==str2num(IDstr)))+1;
                else
                    numlayer=1;
                end
            else % create new layer in list
                % Color management
                S.col = [S.col; uisetcolor([1 1 0],'Select a color')];
                laynum=length(S.col(:,1));
                numlayer=1;  % number of line for this layer
                if flag==1
                    S.layerlist.String=[S.layerlist.String; {['<HTML><FONT color="',rgb2Hex(S.col(laynum,:).*255),'">',[IDstr,' - ',name,' (line)'],'</FONT></HTML>']}]; % new entry in layerlist
                else
                    S.layerlist.String=[S.layerlist.String; {['<HTML><FONT color="',rgb2Hex(S.col(laynum,:).*255),'">',[IDstr,' - ',name,' (points)'],'</FONT></HTML>']}]; % new entry in layerlist
                end
            end
        end
        guidata(gcbf,S); % Update
        % picking:
        if flag==1 % line
            b=1; % initialize button
            xx=[];
            zz=[];
            while ~isempty(b)
                [xt,zt,b]=ginput(1); % picking
                xx=[xx; xt]; % append to list
                zz=[zz; zt];
                plot(S.ax2,xx,zz,'-.','Marker','+','Linewidth',2,'Color',S.col(laynum,:),'Tag','temp');
            end
            obj=get(S.ax2,'Children');
            for i=1:length(obj)-1
                if ~isempty(obj(i).Tag) && strcmp(obj(i).Tag,'temp')
                    delete(obj(i)); 
                end
            end
            plot(S.ax2,xx,zz,'Linewidth',2,'Color',S.col(laynum,:),'ButtonDownFcn',@delone_call,'Tag',IDstr,'UserData',numlayer);
            text(S.ax2,xx(1),zz(1)+10*S.dtz,IDstr,'Color',S.col(laynum,:));
            E=interp1(S.x{S.prof.Value},S.global_coords{S.prof.Value}(:,1),xx);
            N=interp1(S.x{S.prof.Value},S.global_coords{S.prof.Value}(:,2),xx);
        else % points
            b=1; % initialize button
            xx=[];
            zz=[];
            while ~isempty(b)
                [xt,zt,b]=ginput(1); % picking
                xx=[xx; xt]; % append to list
                zz=[zz; zt];
                plot(S.ax2,xx,zz,'+','Linewidth',2,'Color',S.col(laynum,:),'Tag','temp');
            end
            obj=get(S.ax2,'Children');
            for i=1:length(obj)-1
                if ~isempty(obj(i).Tag) && strcmp(obj(i).Tag,'temp')
                    delete(obj(i)); 
                end
            end
            plot(S.ax2,xx,zz,'*','Linewidth',2,'Color',S.col(laynum,:),'ButtonDownFcn',@delone_call,'Tag',IDstr,'UserData',numlayer);
            for xi=1:length(xx)
                text(S.ax2,xx(xi),zz(xi)+10*S.dtz,IDstr,'Color',S.col(laynum,:));
            end
            E=interp1(S.x{S.prof.Value},S.global_coords{S.prof.Value}(:,1),xx);
            N=interp1(S.x{S.prof.Value},S.global_coords{S.prof.Value}(:,2),xx);
        end
        % store picks: x z E N ID Profile# numlayer flag(1=line,2=point)
        S.allpicks=[S.allpicks; [xx zz E N zeros(length(E),1)+str2num(IDstr) zeros(length(E),1)+S.prof.Value zeros(length(E),1)+numlayer zeros(length(E),1)+flag]];
        
        guidata(gcbf,S); % Update
    end

    function [] = layerlist_call(varargin)
        % Callback for layer list -> choose layer
        S=guidata(gcbf);
        w=S.layerlist.String{S.layerlist.Value}(28:end-14);
        i=1;
        while ~isempty(str2num(w(1:i)))
            num=str2num(w(1:i)); % get ID number
            i=i+1;
        end
        S.id.String=num;
        % get name
        name=extractBetween(w,'- ',' (');
        S.name.String=name{1};     
        % get line/point
        flag=extractBetween(w,'(',')');
        if strcmp(flag,'line')
            S.Switchbutton.Value=1;
        else
            S.Switchbutton.Value=2;
        end
        guidata(gcbf,S); % Update
    end

    function [] = delete_layer(varargin)
        % Callback for delete current layer from list
        S=guidata(gcbf);
        w=S.layerlist.String{S.layerlist.Value}(28:end-14);
        % security question:
        opts.Default='Yes';
        opts.Interpreter='none';
        sel=questdlg('Do you really want to delete this layer and all picks of this layer?','Confirm deletion','Yes','No',opts);
        if strcmp(sel,'Yes')
            i=1;
            while ~isempty(str2num(w(1:i)))
                num=str2num(w(1:i)); % get ID number
                i=i+1;
            end
            S.layerlist.String(S.layerlist.Value)=[]; % delete from list
            S.col(S.layerlist.Value,:)=[]; % delete from colors
            S.allpicks(S.allpicks(:,5)==num,:)=[];  % delete from allpicks
            % delete line(s) in picking profile
            obj=get(S.ax2,'Children');
            nobj=length(obj);
            for o=1:nobj-1
                if strcmp(obj(o).Tag,num2str(num))
                    delete(obj(o));
                    delete(obj(o-1)); % corresponding text
                end
            end
            % delete lines in comparison profile
            obj=get(S.ax1,'Children');
            nobj=length(obj);
            for o=1:nobj-1
                if strcmp(obj(o).Tag,num2str(num))
                    delete(obj(o));
                    delete(obj(o-1)); % corresponding text
                end
            end
            guidata(gcbf,S); % Update
        end
    end



    function [] = name_call(varargin)
        % Callback for name -> change name (if layer is selected in list
        % and then only name is changed)
        S=guidata(gcbf);
        l_list=S.layerlist.String;
        for k=1:length(l_list)
            l_list{k}=l_list{k}(28:end-14); % extract only "ID - name (flag)"
            flag{k}=extractBetween(l_list{k},'(',')');
        end
        name=S.name.String;
        IDstr=S.id.String;
        for i=1:length(l_list)
            % get ID
            ii=1;
            while ~isempty(str2num(l_list{i}(1:ii)))
                ii=ii+1;
            end
            id=l_list{i}(1:ii-2);   % id from l_list string
            in(i)=strcmp(id,IDstr); % check if layer ID already used
        end
        if exist('in','var') && any(in==1)
            w=find(in==1);
            % change name of layer
            if strcmp(flag{w},'line') % line
                S.layerlist.String{w}=['<HTML><FONT color="',rgb2Hex(S.col(w,:).*255),'">',[IDstr,' - ',name,' (line)'],'</FONT></HTML>'];
            else % points
                S.layerlist.String{w}=['<HTML><FONT color="',rgb2Hex(S.col(w,:).*255),'">',[IDstr,' - ',name,' (points)'],'</FONT></HTML>'];
            end
        end
        guidata(gcbf,S); % Update
    end

    function [] = delone_call(src,varargin)
         % Callback for delete one layer (Mouse click)
         S=guidata(gcbf);
         % security question:
         opts.Default='Yes';
        opts.Interpreter='none';
        sel=questdlg('Do you really want to delete this line?','Confirm deletion','Yes','No',opts);
        if strcmp(sel,'Yes')
            % get ID and number of line
            ID=str2num(src.Tag);
            numline=src.UserData;
            % delete line handle and corresponding text
            obj=get(S.ax2,'Children');
            for i=1:length(obj)-1
                if ~isempty(obj(i).Tag) && str2num(obj(i).Tag)==ID && obj(i).UserData==numline
                    delete(obj(i)); 
                    delete(obj(i-1));
                end
            end
            % delete picks in S.allpicks
            S.allpicks(S.allpicks(:,6)==S.prof.Value & S.allpicks(:,5)==ID & S.allpicks(:,7)==numline,:)=[]; 
            % check if there exist any other picks of this layer, otherwise
            % delete layer, too:
            if all((S.allpicks(:,5)==ID)==0) % no other picks of this layer...
                for j=1:length(S.layerlist.String)
                    w=S.layerlist.String{j}(28:end-14);
                    i=1;
                    while ~isempty(str2num(w(1:i)))
                        num=str2num(w(1:i)); % get ID number
                        i=i+1;
                    end
                    if num==ID
                        S.layerlist.String(j)=[]; % delete from list
                        S.col(j,:)=[]; % delete from colors
                    end
                end
            end
        end
        guidata(S.fh,S); % Update
    end



    function [] = delall_call(varargin)
        % Callback for delete all layers
        S=guidata(gcbf);
        % security question:
        opts.Default='Yes';
        opts.Interpreter='none';
        sel=questdlg('Do you really want to delete all picks in this profile?','Confirm deletion','Yes','No',opts);
        if strcmp(sel,'Yes')
            % delete all pick lines
            obj=get(S.ax2,'Children');
            nobj=length(obj);
            for o=1:nobj-1
                if ~isempty(str2num(obj(o).Tag))
                    tags(o)=str2num(obj(o).Tag); % get tags
                else
                    tags(o)=NaN;
                end
                delete(obj(o));
            end
            % delete picks in S.allpicks
            S.allpicks(S.allpicks(:,6)==S.prof.Value,:)=[];
            
            % check if there exist any other picks of these layers, otherwise
            % delete layer, too:
            tag=unique(tags(~isnan(tags)));
            for ii=1:length(tag)
                if all((S.allpicks(:,5)==tag(ii))==0) % no other picks of this layer...
                    j=1;
                    while j<=length(S.layerlist.String)
                        w=S.layerlist.String{j}(28:end-14);
                        i=1;
                        while ~isempty(str2num(w(1:i)))
                            num=str2num(w(1:i)); % get ID number
                            i=i+1;
                        end
                        if num==tag(ii)
                            S.layerlist.String(j)=[]; % delete from list
                            S.col(j,:)=[]; % delete from colors
                        end
                        j=j+1;
                    end
                end
            end
            guidata(gcbf,S); % Update
        end
    end

    function [] = save_picks(varargin)
        % Callback for pushbutton saving of all picks
        S=guidata(gcbf);
        [file,pa]= uiputfile('*.txt','Give file name for saving',fullfile(S.folder,'picks.txt'));
        fid=fopen(fullfile(pa,file),'wt');
        if fid>0
            fprintf(fid,'%d\n',length(S.layerlist.String)); % number of layers
            l_list=S.layerlist.String;
            for k=1:length(l_list)
                l_list{k}=l_list{k}(28:end-14); % extract only ID - name
                fprintf(fid,'%s\n',l_list{k}); % layers
            end
            fprintf(fid,'%4.1f\t%4.1f\t%4.1f\n',S.col'); % layer colors
            if isnan(S.z)
                fprintf(fid,'x[m]\tt[ns]\tEasting[m]\tNorthing[m]\tID\tProfile\tlinenumber\tflag\n');
            else
                fprintf(fid,'x[m]\tz[m]\tEasting[m]\tNorthing[m]\tID\tProfile\tlinenumber\tflag\n');
            end
            fprintf(fid,'%8.2f\t%8.2f\t%8.2f\t%8.2f\t%d\t%d\t%d\t%d\n',S.allpicks');
            fclose(fid);
        end
    end

    function [] = load_picks(varargin)
        % Callback for loading of picks
        S=guidata(gcbf);
        [file,pa]=uigetfile('*.txt','Select pick file',S.folder);
        % security question:
        opts.Default='Ok';
        opts.Interpreter='none';
        sel=questdlg('All current picks in this session will be lost.','Confirm loading of new picks','Ok','Cancel',opts);
        if strcmp(sel,'Ok')
            % delete all old stuff
            S.layerlist.String=[];
            S.col=[];
            S.allpicks=[];
            % read new file
            fid=fopen(fullfile(pa,file),'r');
            anz=fscanf(fid,'%d',1);
            for i=1:anz
                temp{i}=textscan(fid,'%s',4); % ID - name (flag)
            end
            tempcol=textscan(fid,'%f%f%f',anz);
            S.col=[tempcol{1} tempcol{2} tempcol{3}];
            % make layer strings with color
            for i=1:anz
                S.layerlist.String{i}=['<HTML><FONT color="',rgb2Hex(S.col(i,:).*255),'">',[temp{i}{1}{1},' ',temp{i}{1}{2},' ',temp{i}{1}{3},' ',temp{i}{1}{4}],'</FONT></HTML>'];
            end
            temp=textscan(fid,'%f%f%f%f%d%d%d%d','Headerlines',2);
            S.allpicks=[temp{1} temp{2} temp{3} temp{4} double(temp{5}) double(temp{6}) double(temp{7}) double(temp{8})];
            fclose(fid);
            guidata(gcbf,S); % Update
            % plot picks in ax1
            compare_call();
            % plot picks in ax2
            num=S.prof.Value;
            if any(S.allpicks(:,6)==num)
                set(S.ax2,'DataAspectratio',[str2num(S.asp_String{S.asp.Value}) 1 1],'nextplot','add');
                picks=S.allpicks(S.allpicks(:,6)==num,:); % picks of this profile
                n=unique(picks(:,5)); % layer IDs
                % get color of layers
                l_list=S.layerlist.String;
                for k=1:length(l_list)
                    l_list{k}=l_list{k}(28:end-14); % extract only ID - name
                end
                for j=1:length(n)
                    for i=1:length(l_list)
                        % get ID
                        ii=1;
                        while ~isempty(str2num(l_list{i}(1:ii)))
                            ii=ii+1;
                        end
                        id=l_list{i}(1:ii-2);   % id from l_list string
                        in2(i,j)=strcmp(id,num2str(n(j))); % check if layer ID already used
                    end
                    w(j)=find(in2(:,j)==1); % which color
                end
                for i=1:length(n)
                    % get possible more lines per layer
                    numline=unique(picks(picks(:,5)==n(i),7));
                    for j=1:length(numline)
                        xtemp=picks(picks(:,5)==n(i) & picks(:,7)==numline(j),1);
                        ytemp=picks(picks(:,5)==n(i) & picks(:,7)==numline(j),2);
                        flag=mean(picks(picks(:,5)==n(i) & picks(:,7)==numline(j),8));
                        if flag==1 % line
                            plot(S.ax2,xtemp,ytemp,'Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numline(j));
                            text(S.ax2,xtemp(1),ytemp(1)+10*S.dtz,num2str(n(i)),'Color',S.col(w(i),:));
                        else % points
                            plot(S.ax2,xtemp,ytemp,'*','Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numline(j));
                            for xi=1:length(xtemp)
                                text(S.ax2,xtemp(xi),ytemp(xi)+10*S.dtz,num2str(n(i)),'Color',S.col(w(i),:));
                            end
                        end
                    end
                end
            end
        end
        guidata(gcbf,S); % Update
        
    end

    function [] = compare_call(varargin)
        % Callback for compare profile number
        S=guidata(gcbf);
        num=str2num(S.comparenum.String); % profile number
        % delete pick lines
        obj=get(S.ax1,'Children');
        nobj=length(obj);
        for o=1:nobj-1
            delete(obj(o));
        end
        % plot new profile for comparison (top plot)
        set(S.Plot1,'CData',S.radargrams{num},'xdata',S.x{num});
        set(S.ax1,'xlim',[min(S.x{num}) max(S.x{num})]);
        % plot picks if available
        if ~isempty(S.allpicks) && any(S.allpicks(:,6)==num)
            set(S.ax1,'DataAspectratio',[str2num(S.asp_String{S.asp.Value}) 1 1],'nextplot','add');
            picks=S.allpicks(S.allpicks(:,6)==num,:); % picks of this profile
            n=unique(picks(:,5)); % layer IDs
            % find which layer
            l_list=S.layerlist.String;
            for k=1:length(l_list)
                l_list{k}=l_list{k}(28:end-14); % extract only ID - name
            end
            for j=1:length(n)
                for i=1:length(l_list)
                    % get ID
                    ii=1;
                    while ~isempty(str2num(l_list{i}(1:ii)))
                        ii=ii+1;
                    end
                    id=l_list{i}(1:ii-2);   % id from l_list string
                    in2(i,j)=strcmp(id,num2str(n(j))); % check if layer ID already used
                end
                w(j)=find(in2(:,j)==1);
            end
            % plot picks
            for i=1:length(n)
                % get possible more lines per layer
                numline=unique(picks(picks(:,5)==n(i),7));
                for j=1:length(numline)
                    xtemp=picks(picks(:,5)==n(i) & picks(:,7)==numline(j),1);
                    ytemp=picks(picks(:,5)==n(i) & picks(:,7)==numline(j),2);
                    flag=mean(picks(picks(:,5)==n(i) & picks(:,7)==numline(j),8));
                    if flag==1 % line
                        plot(S.ax1,xtemp,ytemp,'Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numline(j));
                        text(S.ax1,xtemp(1),ytemp(1)+10*S.dtz,num2str(n(i)),'Color',S.col(w(i),:));
                    else % points
                        plot(S.ax1,xtemp,ytemp,'*','Linewidth',2,'Color',S.col(w(i),:),'ButtonDownFcn',@delone_call,'Tag',num2str(n(i)),'UserData',numline(j));
                        for xi=1:length(xtemp)
                            text(S.ax1,xtemp(xi),ytemp(xi)+10*S.dtz,num2str(n(i)),'Color',S.col(w(i),:));
                        end
                    end
                end
            end
        end
        guidata(gcbf,S); % Update
    end

end

function hexStr = rgb2Hex( rgbColour )
hexStr = reshape( dec2hex( round(rgbColour), 2 )',1, 6);
end