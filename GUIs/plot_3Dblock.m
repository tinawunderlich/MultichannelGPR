function [] = plot_3Dblock(data,xgrid,ygrid,tgrid,pfad)

% plot_3Dblock(data,xgrid,ygrid,tgrid,pfad)
%
% GUI for plotting of 3D block
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3D-matrix with traces along third dimension
% xgrid, ygrid: corresponding coordinate ngrids to data, dx=dy! [m]]
% tgrid: travel time [ns], same size as xgrid/ygrid with time along third dimension
% pfad: complete path 


% Path
S.pfad=pfad;

% set data
S.x = xgrid;
S.y = ygrid;
S.t = tgrid;
S.dt=abs(S.t(1,1,2)-S.t(1,1,1));
S.data=data;
S.dx=xgrid(1,2,1)-xgrid(1,1,1);

% calculate gained data
gainf=[-20 0 15 25 30];
g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
S.datagain=S.data.*(10.^(g./20)); % apply gain

% colorscale limits
coldata=sort(unique(S.datagain(~isnan(S.datagain(:)))));
S.cmingain1=coldata(round(length(coldata)/100*1));
S.cmaxgain1=coldata(end-round(length(coldata)/100*1));
S.cmingain3=coldata(round(length(coldata)/100*3));
S.cmaxgain3=coldata(end-round(length(coldata)/100*3));

coldata=sort(unique(S.data(~isnan(S.data(:)))));
S.cmin1=coldata(round(length(coldata)/100*1));
S.cmax1=coldata(end-round(length(coldata)/100*1));
S.cmin3=coldata(round(length(coldata)/100*3));
S.cmax3=coldata(end-round(length(coldata)/100*3));


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
            'name','3D-Block',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui);
        S.fh.MenuBar = 'figure';
        
        S.ax = axes('unit','pix',...
            'position',[350 150 S.siz(3)-400 S.siz(4)-200]);
        % make initial plot
        S.xloc=S.x(1,round(length(S.x(1,:,1))/2),1); % x-location of slice
        S.yloc=S.y(round(length(S.y(:,1,1))/2),1,1); % y-location of slice
        S.tloc=max(S.t(:)); % t-location of slice
        % Plot slices:
        surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
        hold on
        surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
        surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
        shading interp
        plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 S.tloc],'k') % vertical crossing line of slices
        plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'k') % horizontal line 1
        plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'k') % horizontal line 2
        colorbar('east');
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('t [ns]')
        axis xy
        set(S.ax,'DataAspectratio',[1 1 10],'Zdir','reverse')
        grid on
        axis tight
        
        
        % UI control - Slider x
        S.slx = uicontrol(S.fh,'style','slide',...
            'unit','pix',...
            'position',[30 20 120 30],...
            'min',min(S.x(1,:,1)),'max',max(S.x(1,:,1)),'val',S.xloc,...
            'sliderstep',[S.dx/(max(S.x(1,:,1))-min(S.x(1,:,1))) 1/(max(S.x(1,:,1))-min(S.x(1,:,1)))],'Callback',@slx_callback);
        S.cb_x = uicontrol(S.fh,'style','checkbox',...
            'unit','pix','Value',1,...
            'position',[30 80 150 20],'String','x-location of slice [m]','Callback',@cbx_call);
        S.text_xloc=uicontrol(S.fh,'style','edit','unit','pix','position',[30 55 100 20],'String',num2str(S.xloc,2),'Callback',@xloc_edit);
        
        % UI control - Slider y
        S.sly = uicontrol(S.fh,'style','slide',...
            'unit','pix',...
            'position',[180 20 120 30],...
            'min',min(S.y(:,1,1)),'max',max(S.y(:,1,1)),'val',S.yloc,...
            'sliderstep',[S.dx/(max(S.y(:,1,1))-min(S.y(:,1,1))) 1/(max(S.y(:,1,1))-min(S.y(:,1,1)))],'Callback',@sly_callback);
        S.cb_y = uicontrol(S.fh,'style','checkbox',...
            'unit','pix','Value',1,...
            'position',[180 80 150 20],'String','y-location of slice [m]','Callback',@cby_call);
        S.text_yloc=uicontrol(S.fh,'style','edit','unit','pix','position',[180 55 100 20],'String',num2str(S.yloc,2),'Callback',@yloc_edit);
        
        % UI control - Slider t
        S.slt = uicontrol(S.fh,'style','slide',...
            'unit','pix',...
            'position',[330 20 120 30],...
            'min',0,'max',max(S.t(1,1,:)),'val',S.tloc,...
            'sliderstep',[S.dt/max(S.t(1,1,:)) 1/max(S.t(1,1,:))],'Callback',@slt_callback);
        S.cb_t = uicontrol(S.fh,'style','checkbox',...
            'unit','pix','Value',1,...
            'position',[330 80 150 20],'String','t-location of slice [ns]','Callback',@cbt_call);
        S.text_tloc=uicontrol(S.fh,'style','edit','unit','pix','position',[330 55 100 20],'String',num2str(S.tloc,2),'Callback',@tloc_edit);
        
        
        % UI control aspectratio
        S.asp_String=[{'1/1'} {'1/5'} {'1/10'} {'1/20'} {'1/30'} ];
        S.aspval=1/20;
        S.asp = uicontrol(S.fh,'style','listbox','unit','pix','position',[30 150 60 60],'String',S.asp_String,'Value',4,'Callback',@asp_call);
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[30 215 80 20],'String','Aspect ratio t','HorizontalAlignment','left');
        

        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[30 330 150 15],'Value',1,'Callback',@colorscale_call); % this button is on
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[30 310 150 15],'Callback',@colorscale_call);
        S.rb3=uicontrol('Style','radiobutton','String','3 %','Position',[30 290 150 15],'Callback',@colorscale_call);
        S.flag_cs=1;
        
        
        % UI control: gain
        S.gain_cb=uicontrol(S.fh,'Style','checkbox','unit','pix','Position',[30 500 100 15],'String','Apply gain [dB]','Value',0,'Callback',@gain_call);
        S.gain1 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[30 480 60 15],'String','-20','Callback',@gain_call);
        S.gain2 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[30 460 60 15],'String','0','Callback',@gain_call);
        S.gain3 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[30 440 60 15],'String','15','Callback',@gain_call);
        S.gain4 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[30 420 60 15],'String','25','Callback',@gain_call);
        S.gain5 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[30 400 60 15],'String','30','Callback',@gain_call);
        
        % UI control: save picture
        S.save = uicontrol(S.fh,'Style',...
            'pushbutton','unit','pix',...
            'String','Save picture',...
            'Position',[480 50 100 20],...
            'Callback',@save_call);
        S.savefolder=uicontrol(S.fh,'Style','edit','Position',[480 20 100 20],'String','Figures');
        S.textsave=uicontrol(S.fh,'Style','text','position',[590 20 150 20],'String','<- Give folder name for saving');
        
        
        % UI control: view
        S.r1 = uicontrol(S.fh,'Style',...
            'radiobutton','unit','pix',...
            'String','3D view',...
            'Position',[125 230 100 20],...
            'HandleVisibility','on','Visible','on','Value',1,'Callback',@rb_call);
        S.r2 = uicontrol(S.fh,'Style','radiobutton','unit','pix',...
            'String','xt-plane',...
            'Position',[125 205 100 20],...
            'HandleVisibility','on','Visible','on','Callback',@rb_call);
        S.r3 = uicontrol(S.fh,'Style','radiobutton','unit','pix',...
            'String','yt-plane',...
            'Position',[125 180 100 20],...
            'HandleVisibility','on','Visible','on','Callback',@rb_call);
        S.r4 = uicontrol(S.fh,'Style','radiobutton','unit','pix',...
            'String','xy-plane (Tsl)',...
            'Position',[125 155 100 20],...
            'HandleVisibility','on','Visible','on','Callback',@rb_call);
        S.flag=1;
        
    end


% save handles
guidata(S.fh,S);


%%%---------- Callback functions ----------------


    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        wid=wid_fig-400; % width of axes
        hei_fig=S.fh.Position(4); % figure height
        hei=hei_fig-200;    % height of axes
        
        % change size of axes
        S.ax.Position=[350 150 wid hei];
    end


    function [] = save_call(varargin)
        % Callback for saving of picture
        S=guidata(gcbf);
        % create invisible figure
        f2=figure('Visible','off','position',S.siz);
        if S.cb_x.Value==1 && S.cb_y.Value==1 && S.cb_t.Value==1            % X Y T
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 max(S.t(1,1,:))],'r') % vertical crossing line of slices
            plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'r') % horizontal line 1
            plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'r') % horizontal line 2
            % name for saving
            n1=['x',num2str(S.xloc,4),'m_y',num2str(S.yloc,4),'m_t',num2str(S.tloc,4),'ns'];
        elseif S.cb_x.Value==1 && S.cb_y.Value==1 && S.cb_t.Value==0        % X Y
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            end
            shading interp
            plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 max(S.t(1,1,:))],'r') % vertical crossing line of slices
            % name for saving
            n1=['x',num2str(S.xloc,4),'m_y',num2str(S.yloc,4),'m'];
        elseif S.cb_x.Value==1 && S.cb_y.Value==0 && S.cb_t.Value==0        % X
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
            end
            shading interp
            % name for saving
            n1=['x',num2str(S.xloc,4),'m'];
        elseif S.cb_x.Value==0 && S.cb_y.Value==1 && S.cb_t.Value==0        % Y
            if S.gain_cb.Value==1
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            else
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            end
            shading interp
            % name for saving
            n1=['y',num2str(S.yloc,4),'m'];
        elseif S.cb_x.Value==0 && S.cb_y.Value==1 && S.cb_t.Value==1        % Y T
            if S.gain_cb.Value==1
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'r') % horizontal line 2
            % name for saving
            n1=['y',num2str(S.yloc,4),'m_t',num2str(S.tloc,4),'ns'];
        elseif S.cb_x.Value==0 && S.cb_y.Value==0 && S.cb_t.Value==1        % T
            if S.gain_cb.Value==1
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            % name for saving
            n1=['t',num2str(S.tloc,4),'ns'];
        elseif S.cb_x.Value==1 && S.cb_y.Value==0 && S.cb_t.Value==1        % X T
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'r') % horizontal line 1
            % name for saving
            n1=['x',num2str(S.xloc,4),'m_t',num2str(S.tloc,4),'ns'];
        end
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('t [ns]')
        axis xy
        set(gca,'DataAspectratio',[1 1 1/S.aspval],'Zdir','reverse')
        grid on
        axis tight
        % set colorscale
        if S.gain_cb.Value==1
            if S.flag_cs==1
                set(gca,'ClimMode','auto');
            elseif S.flag_cs==2
                set(gca,'ClimMode','manual','CLim',[S.cmingain1 S.cmaxgain1]);
            else
                set(gca,'ClimMode','manual','CLim',[S.cmingain3 S.cmaxgain3]);
            end
        else
            if S.flag_cs==1
                set(gca,'ClimMode','auto');
            elseif S.flag_cs==2
                set(gca,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
            else
                set(gca,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
            end
        end
        % set view
        if S.flag==1 % 3d view
            % get view
            view(S.ax.View);
            n2='3D'; % name for saving
        elseif S.flag==2 % xt plane
            view(0, 0);
            n2='xt'; % name for saving
        elseif S.flag==3 % yt plane
            view(90, 0);
            n2='yt'; % name for saving
        elseif S.flag==4 % Tsl
            view(90, 90);
            n2='Tsl'; % name for saving
        end
        % save figure
        if ~exist(fullfile(S.pfad,S.savefolder.String),'dir')
            mkdir(fullfile(S.pfad,S.savefolder.String));
        end
        f2.PaperPositionMode='auto';
        saveas(f2,fullfile(S.pfad,S.savefolder.String,['Plot_',n2,'_',n1,'.png']),'png');
    end



    function [] = colorscale_call(varargin)
        % Callback for colorscale limits
        S=guidata(gcbf);
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
            if S.gain_cb.Value==1
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain1 S.cmaxgain1]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
            end
            drawnow;
            S.flag_cs=2;
        elseif val3==1 && S.flag_cs~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            if S.gain_cb.Value==1
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain3 S.cmaxgain3]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
            end
            drawnow;
            S.flag_cs=3;
        end
        guidata(gcbf,S); % Update (for flag_cs)
    end


    function [] = gain_call(varargin)
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch');
        S=guidata(gcbf);
        % calculate gain new
        gainf=[str2num(S.gain1.String) str2num(S.gain2.String) str2num(S.gain3.String) str2num(S.gain4.String) str2num(S.gain5.String)];
        g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
        S.datagain=S.data.*(10.^(g./20)); % apply gain
        guidata(gcbf,S); % Update
        plotnew(varargin);
    end


    function [] = rb_call(varargin) % callback for view
        S=guidata(gcbf);
        if S.r1.Value==1 && S.flag~=1
            S.r2.Value=0;
            S.r3.Value=0;
            S.r4.Value=0;
            S.flag=1;
        elseif S.r2.Value==1 && S.flag~=2
            S.r1.Value=0;
            S.r3.Value=0;
            S.r4.Value=0;
            S.flag=2;
        elseif S.r3.Value==1 && S.flag~=3
            S.r1.Value=0;
            S.r2.Value=0;
            S.r4.Value=0;
            S.flag=3;
        elseif S.r4.Value==1 && S.flag~=4
            S.r1.Value=0;
            S.r3.Value=0;
            S.r2.Value=0;
            S.flag=4;
        end
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch');
        guidata(gcbf,S); % Update
        plotnew(varargin);
    end


    function [] = plotnew(varargin)
        % plot function
        S=guidata(gcbf);
        S.save.Enable='on';
        S.savefolder.Enable='on';
        % plot new
        hold off
        if S.cb_x.Value==1 && S.cb_y.Value==1 && S.cb_t.Value==1            % X Y T
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 max(S.t(1,1,:))],'k') % vertical crossing line of slices
            plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'k') % horizontal line 1
            plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'k') % horizontal line 2
        elseif S.cb_x.Value==1 && S.cb_y.Value==1 && S.cb_t.Value==0        % X Y
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            end
            shading interp
            plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 max(S.t(1,1,:))],'k') % vertical crossing line of slices
        elseif S.cb_x.Value==1 && S.cb_y.Value==0 && S.cb_t.Value==0        % X
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
            end
            shading interp
        elseif S.cb_x.Value==0 && S.cb_y.Value==0 && S.cb_t.Value==0        % nix
            S.save.Enable='off';
            S.savefolder.Enable='off';
            hold off
            plot3([S.xloc S.xloc],[S.yloc S.yloc],[0 max(S.t(1,1,:))],'k') % vertical crossing line of slices
            hold on
            plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'k') % horizontal line 1
            plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'k') % horizontal line 2
        elseif S.cb_x.Value==0 && S.cb_y.Value==1 && S.cb_t.Value==0        % Y
            if S.gain_cb.Value==1
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            else
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
            end
            shading interp
        elseif S.cb_x.Value==0 && S.cb_y.Value==1 && S.cb_t.Value==1        % Y T
            if S.gain_cb.Value==1
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.datagain(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                hold on
                surf(permute(S.x(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.y(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.t(S.y(:,1,1)==S.yloc,:,:),[2 3 1]),permute(S.data(S.y(:,1,1)==S.yloc,:,:),[2 3 1]))
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            hold on
            plot3([min(S.x(1,:,1)) max(S.x(1,:,1))],[S.yloc S.yloc],[S.tloc S.tloc],'k') % horizontal line 2
        elseif S.cb_x.Value==0 && S.cb_y.Value==0 && S.cb_t.Value==1        % T
            if S.gain_cb.Value==1
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
        elseif S.cb_x.Value==1 && S.cb_y.Value==0 && S.cb_t.Value==1        % X T
            if S.gain_cb.Value==1
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.datagain(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.datagain(:,:,S.t(1,1,:)==S.tloc))
            else
                surf(permute(S.x(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.y(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.t(:,S.x(1,:,1)==S.xloc,:),[1 3 2]),permute(S.data(:,S.x(1,:,1)==S.xloc,:),[1 3 2]))
                hold on
                surf(S.x(:,:,S.t(1,1,:)==S.tloc),S.y(:,:,S.t(1,1,:)==S.tloc),S.t(:,:,S.t(1,1,:)==S.tloc),S.data(:,:,S.t(1,1,:)==S.tloc))
            end
            shading interp
            plot3([S.xloc S.xloc],[min(S.y(:,1,1)) max(S.y(:,1,1))],[S.tloc S.tloc],'k') % horizontal line 1
        end
        colorbar('east');
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('t [ns]')
        axis xy
        set(S.ax,'DataAspectratio',[1 1 1/S.aspval],'Zdir','reverse')
        grid on
        axis tight
        % set colorscale
        if S.gain_cb.Value==1
            if S.flag_cs==1
                set(S.ax,'ClimMode','auto');
            elseif S.flag_cs==2
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain1 S.cmaxgain1]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain3 S.cmaxgain3]);
            end
        else
            if S.flag_cs==1
                set(S.ax,'ClimMode','auto');
            elseif S.flag_cs==2
                set(S.ax,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
            end
        end
        % set view
        if S.flag==1 % 3d view
            S.ax.View=[-37.5 30]; % default
        elseif S.flag==2 % xt plane
            S.ax.View=[0 0];
        elseif S.flag==3 % yt plane
            S.ax.View=[90 0];
        elseif S.flag==4 % Tsl
            S.ax.View=[90 90];
        end
        drawnow;
        % change pointer back to arrow
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'arrow')
        guidata(gcbf,S); % Update
    end

    function [] = cbx_call(varargin)
        % Callback for cb x
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        if S.cb_x.Value==0
            S.slx.Enable='Off';
        else
            S.slx.Enable='On';
        end
        guidata(gcbf,S); % Update
        plotnew(varargin);
    end

    function [] = cby_call(varargin)
        % Callback for cb y
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        if S.cb_y.Value==0
            S.sly.Enable='Off';
        else
            S.sly.Enable='On';
        end
        guidata(gcbf,S); % Update
        plotnew(varargin);
    end


    function [] = cbt_call(varargin)
        % Callback for cb t
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        if S.cb_t.Value==0
            S.slt.Enable='Off';
        else
            S.slt.Enable='On';
        end
        guidata(gcbf,S); % Update
        plotnew(varargin);
    end



    function [] = slx_callback(varargin)
        % Callback for slider x
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.x(1,min(abs(S.x(1,:,1)-S.slx.Value))==abs(S.x(1,:,1)-S.slx.Value),1);
        S.xloc=temp(1);
        S.text_xloc.String=num2str(S.xloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end

    function [] = xloc_edit(varargin)
        % Callback for x loc edit
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.x(1,min(abs(S.x(1,:,1)-str2num(S.text_xloc.String)))==abs(S.x(1,:,1)-str2num(S.text_xloc.String)),1);
        S.xloc=temp(1);
        S.slx.Value=S.xloc;
        S.text_xloc.String=num2str(S.xloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end


    function [] = sly_callback(varargin)
        % Callback for slider y
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.y(min(abs(S.y(:,1,1)-S.sly.Value))==abs(S.y(:,1,1)-S.sly.Value),1,1);
        S.yloc=temp(1);
        S.text_yloc.String=num2str(S.yloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end

    function [] = yloc_edit(varargin)
        % Callback for yloc edit
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.y(min(abs(S.y(:,1,1)-str2num(S.text_yloc.String)))==abs(S.y(:,1,1)-str2num(S.text_yloc.String)),1,1);
        S.yloc=temp(1);
        S.sly.Value=S.yloc;
        S.text_yloc.String=num2str(S.yloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end


    function [] = slt_callback(varargin)
        % Callback for slider t
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.t(1,1,min(abs(S.t(1,1,:)-S.slt.Value))==abs(S.t(1,1,:)-S.slt.Value));
        S.tloc=temp(1);
        S.text_tloc.String=num2str(S.tloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end

    function [] = tloc_edit(varargin)
        % Callback for tloc edit
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        temp=S.t(1,1,min(abs(S.t(1,1,:)-str2num(S.text_tloc.String)))==abs(S.t(1,1,:)-str2num(S.text_tloc.String)));
        S.tloc=temp(1);
        S.slt.Value=S.tloc;
        S.text_tloc.String=num2str(S.tloc,'%8.2f');
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end



    function [] = asp_call(varargin)
        % Callback for aspectratio
        set(findobj('Type','Figure','Name','3D-Block'), 'pointer', 'watch')
        S=guidata(gcbf);
        S.aspval=S.asp_String(S.asp.Value);
        S.aspval=str2num(S.aspval{1});
        guidata(gcbf,S); % Update
        % call plot function
        plotnew(varargin);
    end
end