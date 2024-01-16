function []=plot_TslRadargramm(tsl,xgrid,ygrid,t_tsl,radar,t,coords,profileslist,timedepth)

% GUI for plotting of Tsl and radargrams with connected cursor
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% (created by make_Timelices.m)
% tsl: data of timeslices 
% xgrid/ygrid: corresponding coordinate grids to data, dx=dy! [m]
% t_tsl: matrix of start/end time of each timeslice [ns]
% (created by make_Radargrams.m, Bin2Radargram.m or own creation)
% radar: cells with data of radargrams 
% t: time vector of radargrams [ns]
% coords: cells with starting and ending coordinates of profile (or
% complete profile lists for curved profiles)
% profilelist: number of profiles
% timedepth: flag: time=1, depth=2


S.timedepth=timedepth;
% set tsl
S.x = xgrid;
S.y = ygrid;
S.t_tsl = t_tsl;
S.tsl=tsl;
S.dx=xgrid(1,2)-xgrid(1,1);
S.tslnum=1; % current number of tsl
S.coordtrans=[];
S.ar=1; %aspect ratio factor

% set radargramms (=data)
S.data=radar;
S.t=t;
if timedepth==2
    S.t=abs(S.t-S.t(1)); % depth is positive down, starting at 0m -> corresponds to dsl
end
S.dt=S.t(2)-S.t(1);
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
S.rdgnum='1';


% colorscale limits Tsl
for i=1:length(S.tsl)
    coldata=sort(unique(S.tsl{i}(~isnan(S.tsl{i}(:)))));
    S.cmin1(i)=coldata(round(length(coldata)/100*1));
    S.cmax1(i)=coldata(end-round(length(coldata)/100*1));
    S.cmin3(i)=coldata(round(length(coldata)/100*3));
    S.cmax3(i)=coldata(end-round(length(coldata)/100*3));
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
        
        % create figure handle:
        S.fh = figure('units','pixels',...
            'position',S.siz,...
            'menubar','none',...
            'name','Tsl and Radargramms',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui);
        
        S.fh.MenuBar = 'figure';
        set (S.fh, 'WindowButtonMotionFcn', @mouseMove); % for connecting of Tsl and radargrams
        
        % UI control for tsl-time text
        S.tsltext=uicontrol(S.fh,'style','text','unit','pix','position',[10 400 100 15],'String',['Tsl: ',num2str(S.t_tsl(1,1)),' - ',num2str(S.t_tsl(1,2)),' ns'],'HorizontalAlignment','left','callback',@tsltextcall,'FontWeight','bold');
        
        S.aspect=uicontrol(S.fh,'Style','popupmenu','unit','pix','position',[5 200 110 25],'String',[{'AspectRatio'} {'10/1'} {'5/1'} {'4/1'} {'3/1'} {'2/1'} {'1/1'} {'1/2'} {'1/3'} {'1/4'} {'1/5'} {'1/10'} {'1/20'}],'callback',{@aspect_call});
        
        % UI control for profile number 
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[10 80 100 90],'String',S.profileslist,'Value',1,'callback',{@profnum_call});
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[10 165 100 15],'String','Choose profile','HorizontalAlignment','left','FontWeight','bold');
        
        %%% TSL
        if length(S.x(:,1))>length(S.x(1,:)) % Tsl height>width
            S.ax_tsl = axes('unit','pix',...
            'position',[150 40 (S.siz(3)-200)/3 S.siz(4)-40]);
            S.layout=1;
        else
            S.ax_tsl = axes('unit','pix',...
            'position',[150 S.siz(4)-(S.siz(4)-90)/3*2 S.siz(3)-170 (S.siz(4)-90)/3*2]);
            S.layout=2;
        end
        % make initial plot TSL
        S.h_tsl=imagesc(S.x(1,:),S.y(:,1),S.tsl{1});
        hold on
        if abs(mean(S.x(1,:))-mean(S.xprof{1}))<=1000 && abs(mean(S.y(:,1))-mean(S.yprof{1}))<=1000
            S.line=plot(S.xprof{1},S.yprof{1},'r','Linewidth',2);
        else
            S.line=plot([],[],'r','Linewidth',2);
        end
        colormap(flipud(gray));
        xlabel('x [m]')
        ylabel('y [m]')
        axis xy
        set(S.ax_tsl,'DataAspectratio',[1 1 1])
        grid on
        axis tight
        
        
        %%% RADARGRAM
        if length(S.x(:,1))>length(S.x(1,:)) % Tsl height>width
            S.ax_data=axes('unit','pix',...
                'position',[180+(S.siz(3)-200)/3 40 S.siz(3)-200-(S.siz(3)-200)/3 S.siz(4)/3]);
        else
            S.ax_data=axes('unit','pix',...
                'position',[150 40 S.siz(3)-170 (S.siz(4)-90)/3]);
        end
        % make initial plot RADARGRAM
        S.h_rdg=imagesc(S.loc_coord{1},S.t,S.data{1});
        colormap(flipud(gray));
        xlabel('ProfileX [m]')
        set(S.ax_data,'Xlim',[min(S.loc_coord{1}) max(S.loc_coord{1})])
        if S.timedepth==1
            ylabel('t [ns]')
        else
            ylabel('z [m]')
        end
        axis ij
        grid on
        axis tight
        
        % UI control - radiobuttons for layout
        S.lay1=uicontrol('Style','radiobutton','String','Layout 1','Position',[10 480 100 15],'Callback',@layout_call); 
        S.lay2=uicontrol('Style','radiobutton','String','Layout 2','Position',[10 460 100 15],'Callback',@layout_call);
        if S.layout==1
            S.lay1.Value=1;
        else
            S.lay2.Value=1;
        end
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[10 360 100 15],'Value',1,'Callback',@colorscale_call); % this button is on
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[10 340 100 15],'Callback',@colorscale_call);
        S.rb3=uicontrol('Style','radiobutton','String','3 %','Position',[10 320 100 15],'Callback',@colorscale_call);
        S.flag_cs=1;
        
        % UI control for wiggleplot
        S.cb_wiggle=uicontrol('Style','checkbox','String','Wiggle plot','Position',[10 270 100 30],'Value',0,'Callback',@wiggle_call);
        S.wigscal=uicontrol('Style','Edit','String','1','Position',[10 250 45 25],'Callback',@wiggle_call);
        S.wigscaltext=uicontrol('Style','text','String','Scale','Position',[52 250 50 25]);
        
        % UI control for coordinatetransformation
        S.cb_coordtrans=uicontrol('Style','checkbox','String','Apply coord. transf. to radargrams','Position',[10 10 180 30],'Value',0,'Callback',@coordtrans_call);
        
        % UI control - checkbox for plotting of all radargrams
        S.cb_allrdgs=uicontrol('Style','checkbox','String','Show location of all prof.','Position',[10 40 180 30],'Value',0,'Callback',@allprofiles_call);
    end


% save handles
guidata(S.fh,S);


%%%---------- Callback functions ----------------


    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        hei_fig=S.fh.Position(4); % figure height
        
        % change size of axes
        if S.layout==1 % Tsl height>width: Layout1
            S.ax_tsl.Position = [150 40 (wid_fig-200)/3 hei_fig-60]; 
            S.ax_data.Position=[180+(wid_fig-200)/3 60 wid_fig-200-(wid_fig-200)/3 hei_fig/3];
        else  % Layout2
            S.ax_tsl.Position = [150 hei_fig-(hei_fig-90)/3*2 wid_fig-170 (hei_fig-90)/3*2];
            S.ax_data.Position= [150 40 wid_fig-170 (hei_fig-90)/3];
        end
    end

    function mouseMove (object, eventdata)
        % get position inside radargram
        C = get (S.ax_data, 'CurrentPoint');
        S=guidata(object);
        % find tsl-number:
        n=find(C(1,2)>S.t_tsl(:,1) & C(1,2)<=S.t_tsl(:,2));
        % find global coords of point:
        S.xy=[interp1(S.loc_coord{str2num(S.rdgnum)},S.xprof{str2num(S.rdgnum)},C(1,1)) interp1(S.loc_coord{str2num(S.rdgnum)},S.yprof{str2num(S.rdgnum)},C(1,1))];
        % if mouse inside radargram:
        if C(1,1)>=S.ax_data.XLim(1) && C(1,1)<=S.ax_data.XLim(2) && C(1,2)>=S.ax_data.YLim(1) && C(1,2)<=S.ax_data.YLim(2)
            if S.cb_allrdgs.Value==1
                S.cb_allrdgs.Value=0; % turn "plot all lines" off
            end
            set(gcf,'pointer','crosshair');
            if n~=S.tslnum % if new tsl-number  
                S.tslnum=n;
                guidata(object,S); % Update
                tsltextcall(); % change display of Tsl-times
            end
            
            % update cross in tsl-plot:
            if abs(mean(S.x(1,:))-mean(S.xprof{str2num(S.rdgnum)}))<=1000 && abs(mean(S.y(:,1))-mean(S.yprof{str2num(S.rdgnum)}))<=1000
                obj=get(S.ax_tsl,'Children');
                nobj=length(obj);
                for o=1:nobj-2
                    delete(obj(o)); % delete old lines
                end
                dd=S.dx*15;  % half width of cross
                plot(S.ax_tsl,[S.xy(1)-dd S.xy(1)+dd],[S.xy(2) S.xy(2)],'b','Linewidth',2);
                plot(S.ax_tsl,[S.xy(1) S.xy(1)],[S.xy(2)-dd S.xy(2)+dd],'b','Linewidth',2);
            end
            
            % update title with current coordinates:
            title(S.ax_data, ['(ProfileX,T) = (', num2str(C(1,1),'%8.2f'), 'm, ',num2str(C(1,2),'%8.2f'), 'ns)','   (X,Y) = (', num2str(S.xy(1),'%8.2f'), 'm, ',num2str(S.xy(2),'%8.2f'), 'm)']); 
        else
            set(gcf,'pointer','arrow');
        end
    end

    function [] = coordtrans_call(varargin)
        % callback for plotting of all profile lines
        S=guidata(gcbf);
        if S.cb_coordtrans.Value==1    % if on, get coordtrans-info
            % get path to coordtrans:
            if ~ispc; menu('Choose folder with coordtrans.mat (e.g. 3D_Grid_R*)','OK'); end
            if exist('.cttemp.temp') % read last opened folder from temp.temp
                fid=fopen('.cttemp.temp','r');
                fn=textscan(fid,'%s');
                fclose(fid);
                if ~isempty(fn{1})
                    pfad_ct=uigetdir(fn{1}{1},'Choose folder with coordtrans.mat (e.g. 3D_Grid_R*)');
                else
                    pfad_ct=uigetdir([],'Choose folder with coordtrans.mat (e.g. 3D_Grid_R*)');
                end
            else
                pfad_ct=uigetdir([],'Choose folder with coordtrans.mat (e.g. 3D_Grid_R*)'); % path to tsl-folder
            end
            fid=fopen('.cttemp.temp','w');
            fprintf(fid,pfad_ct);
            fclose(fid);
            
            if exist(fullfile(pfad_ct,'coordtrans.mat'),'file')
                S.coordtrans=load(fullfile(pfad_ct,'coordtrans.mat')); % lokal/global
            else
                S.coordtrans.coordtrans=[1 1 1 1; 2 2 2 2];
            end
            
            % apply coordtrans to radargram (global->lokal)
            for i=1:length(S.xprof)
                temp=helmert([S.xprof{i} S.yprof{i}],S.coordtrans.coordtrans(:,3:4),S.coordtrans.coordtrans(:,1:2));
                S.xprof{i}=temp(:,1);
                S.yprof{i}=temp(:,2);
            end
        else
            % apply coordtrans to radargram (lokal->global)
            for i=1:length(S.xprof)
                temp=helmert([S.xprof{i} S.yprof{i}],S.coordtrans.coordtrans(:,1:2),S.coordtrans.coordtrans(:,3:4));
                S.xprof{i}=temp(:,1);
                S.yprof{i}=temp(:,2);
            end
        end
        guidata(gcbf,S); % Update
    end

    function [] = allprofiles_call(varargin)
        % callback for plotting of all profile lines
        S=guidata(gcbf);
        if S.cb_allrdgs.Value==1    % if on, make profile lines visible
            for i=1:length(S.xprof)
                % if in same coordinate system as Tsl, then plot
                if abs(mean(S.x(1,:))-mean(S.xprof{i}))<=1000 && abs(mean(S.y(:,1))-mean(S.yprof{i}))<=1000
                    hold on
                    S.allline=plot(S.ax_tsl,S.xprof{i},S.yprof{i},'Linewidth',2);
                end
            end
        else
            if isfield(S,'allline')
                delete(S.allline);
                drawnow;
            end
        end
        guidata(gcbf,S); % Update
    end

    function [] = tsltextcall(varargin)
        % callback for tsl-textbox (display tsl-times) and change tsl
        % display
        S=guidata(gcbf);
        n=S.tslnum;
        if S.timedepth==1
            S.tsltext.String=['Tsl: ',num2str(S.t_tsl(n,1)),' - ',num2str(S.t_tsl(n,2)),' ns'];
        else
            S.tsltext.String=['Dsl: ',num2str(S.t_tsl(n,1)),' - ',num2str(S.t_tsl(n,2)),' m'];
        end
        guidata(gcbf,S); % Update
        % call tsl-plotting function
        tslplot(varargin);
    end

    function [] = tslplot(varargin)
        % function for plotting new tsl
        S=guidata(gcbf);
        n=S.tslnum; % tsl-number
        S.h_tsl.CData=S.tsl{n};
        % set new colorscale
        if S.flag_cs==1
            set(S.ax_tsl,'ClimMode','auto');
        elseif S.flag_cs==2
            set(S.ax_tsl,'ClimMode','manual','CLim',[S.cmin1(n) S.cmax1(n)]);
        else
            set(S.ax_tsl,'ClimMode','manual','CLim',[S.cmin3(n) S.cmax3(n)]); 
        end
        drawnow;
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
        profnum_call();
        guidata(gcbf,S); % Update
    end

    function [] = colorscale_call(varargin)
        % Callback for colorscale limits
        S=guidata(gcbf);
        n=S.tslnum; % tsl-number
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1 && S.flag_cs~=1
            S.rb2.Value=0;
            S.rb3.Value=0;
            set(S.ax_tsl,'ClimMode','auto');
            set(S.ax_data,'ClimMode','auto');
            drawnow;
            S.flag_cs=1;
        elseif val2==1 && S.flag_cs~=2
            S.rb1.Value=0;
            S.rb3.Value=0;
            set(S.ax_tsl,'ClimMode','manual','CLim',[S.cmin1(n) S.cmax1(n)]);
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin1R(n) S.cmax1R(n)]);
            drawnow;
            S.flag_cs=2;
        elseif val3==1 && S.flag_cs~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            set(S.ax_tsl,'ClimMode','manual','CLim',[S.cmin3(n) S.cmax3(n)]);
            set(S.ax_data,'ClimMode','manual','CLim',[S.cmin3R(n) S.cmax3R(n)]);
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
            guidata(gcbf,S); % Update
            profnum_call();
        end
    end

    function [] = profnum_call(varargin)
        % new profile from list
        S=guidata(gcbf);
        S.rdgnum=S.profileslist(S.prof.Value,:);
        guidata(gcbf,S); % Update
        % plot new red line in tsl-plot:
        newline(varargin);
        % plot new radargram:
        axes(S.ax_data);
        delete(S.ax_data.Children)
        if S.cb_wiggle.Value==0
            %set(S.h_rdg,'XData',S.loc_coord{str2num(S.rdgnum)},'CData',S.data{str2num(S.rdgnum)});
            S.h_rdg=imagesc(S.loc_coord{str2num(S.rdgnum)},S.t,S.data{str2num(S.rdgnum)});
            colormap(flipud(gray));
        else
            % Wiggle plot
            wigglesc(S.data{str2num(S.rdgnum)},S.t,S.loc_coord{str2num(S.rdgnum)},str2num(S.wigscal.String));
        end
        set(S.ax_data,'DataAspectRatio',[1 S.ar 1],'Xlim',[min(S.loc_coord{str2num(S.rdgnum)}) max(S.loc_coord{str2num(S.rdgnum)})]);
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