function [] = Tsl_slider_plot(xgrid,ygrid,tsl,topo,t,pfad,dsl,maxElevation,coordtrans)

% function [] = Tsl_slider_plot(xgrid,ygrid,tsl,t,pfad,coordtrans)
%
% Plot timeslices according to slider location and with other options
%
% Dr. Tina Wunderlich, CAU Kiel 2019, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% xgrid, ygrid: Coordinate grids in m
% tsl: timeslices (cell array of matrices)
% topo: matrix with topography, same size as xgrid/ygrid
% t: vector with start and end times of timeslices (size=[length(tsl) 2])
% pfad: path of timeslice-folder (for saving of georeferenced pngs)
% (without '/' at the end!)
% dsl: =1 if depthslices instead of timeslices (dsl=0)
% maxElevation: if dsl=1: Elevation (a.s.l) at depth=0m (maxElevation=[] for dsl=0)
% coordtrans: if area has been rotated, here are the coordinate pairs for
% helmert transformation
%
% requires folder Subfunctions

if nargin==8
    ct=0;   % no coordinate transformation

    S.xlocal=xgrid;
    S.ylocal=ygrid;
    S.xglobal=xgrid;
    S.yglobal=ygrid;
    S.dsl=dsl;
    S.maxElevation=maxElevation;
else
    ct=1;   % apply coordinate transformation
    
    % prepare rotated coordinate grids for tsl and topo:
    neu=helmert([xgrid(:) ygrid(:)],coordtrans(:,1:2),coordtrans(:,3:4));
    
    S.xglobal=reshape(neu(:,1),size(xgrid));
    S.yglobal=reshape(neu(:,2),size(ygrid));
    S.xlocal=xgrid;
    S.ylocal=ygrid;
    S.dsl=dsl;
    S.maxElevation=maxElevation;
    
    S.coordtrans=coordtrans;
end

% set data
S.x = xgrid; % for current coordinate system
S.y = ygrid;
S.t = t;
S.tsl=tsl;
S.pfad=pfad;
S.dx=xgrid(1,2)-xgrid(1,1);
S.dy=ygrid(2,1)-ygrid(1,1);
S.ct=ct;    % flag for coordinate transformation
S.topo=topo;
if S.dsl
    S.unit = 'm';
else
    S.unit = 'ns';
end

% determine color limits
for i=1:length(S.tsl)
    coldata=S.tsl{i}(~isnan(S.tsl{i}));
    if ~isempty(coldata)
        S.coldata{i}=sort(unique(coldata));
    else
        S.coldata{i}=[0 1];
    end
end

% determine tick limits and ticks (local coordinate system):
% find min/max of coordinates with data:
minx=min(S.xlocal(:));
maxx=max(S.xlocal(:));
miny=min(S.ylocal(:));
maxy=max(S.ylocal(:));
% find tick spacing which makes approx. 5 ticks per axis:
liste=[0.5 1 2 5 10 20 50 100]; % possible tick spacings
dxtick=liste(min(abs((maxx-minx)./liste-5))==abs((maxx-minx)./liste-5));
dytick=liste(min(abs((maxy-miny)./liste-5))==abs((maxy-miny)./liste-5));
dtick=min(dxtick,dytick); % common tick spacing for both directions
% change min/max coordinate values for border around data:
minx=minx-dtick/6;
maxx=maxx+dtick/6;
miny=miny-dtick/6;
maxy=maxy+dtick/6;
% find first and last tick
test=floor(minx):-1:floor(minx)-dtick;
S.xs=test(mod(test,dtick)==0);
test=ceil(maxx):ceil(maxx)+dtick;
S.xe=test(mod(test,dtick)==0);
test=floor(miny):-1:floor(miny)-dtick;
S.ys=test(mod(test,dtick)==0);
test=ceil(maxy):ceil(maxy)+dtick;
S.ye=test(mod(test,dtick)==0);

% ticks (local):
S.xticks_l=S.xs(1):dtick:S.xe(1);
S.yticks_l=S.ys(1):dtick:S.ye(1);

% Ticklabels (local):
if round(dtick)==dtick % is integer
    for i=1:length(S.xticks_l)
        S.xticklabels_l{i}=sprintf('%d',S.xticks_l(i));
    end
    for i=1:length(S.yticks_l)
        S.yticklabels_l{i}=sprintf('%d',S.yticks_l(i));
    end
else % for floats
    for i=1:length(S.xticks_l)
        S.xticklabels_l{i}=sprintf('%8.1f',S.xticks_l(i));
    end
    for i=1:length(S.yticks_l)
        S.yticklabels_l{i}=sprintf('%8.1f',S.yticks_l(i));
    end
end


% determine tick limits and ticks (global coordinate system):
% find min/max of coordinates with data:
cf=cellfun(@(x) ~isnan(x),S.tsl,'UniformOutput',false); % getting the non-NaNs in each tsl
minx=min(cell2mat(cellfun(@(x) min(S.xglobal(x)),cf,'UniformOutput',false)')); % minx in each tsl
maxx=max(cell2mat(cellfun(@(x) max(S.xglobal(x)),cf,'UniformOutput',false)'));
miny=min(cell2mat(cellfun(@(x) min(S.yglobal(x)),cf,'UniformOutput',false)'));
maxy=max(cell2mat(cellfun(@(x) max(S.yglobal(x)),cf,'UniformOutput',false)'));
% find tick spacing which makes approx. 5 ticks per axis:
liste=[0.5 1 2 5 10 20 50 100]; % possible tick spacings
dxtick=liste(min(abs((maxx-minx)./liste-5))==abs((maxx-minx)./liste-5));
dytick=liste(min(abs((maxy-miny)./liste-5))==abs((maxy-miny)./liste-5));
dtick=min(dxtick,dytick); % common tick spacing for both directions
% change min/max coordinate values for border around data:
minx=minx-dtick/6;
maxx=maxx+dtick/6;
miny=miny-dtick/6;
maxy=maxy+dtick/6;
% find first and last tick
test=floor(minx):-1:floor(minx)-dtick;
S.xsg=test(mod(test,dtick)==0);
test=ceil(maxx):ceil(maxx)+dtick;
S.xeg=test(mod(test,dtick)==0);
test=floor(miny):-1:floor(miny)-dtick;
S.ysg=test(mod(test,dtick)==0);
test=ceil(maxy):ceil(maxy)+dtick;
S.yeg=test(mod(test,dtick)==0);

% ticks (local):
S.xticks_g=S.xsg(1):dtick:S.xeg(1);
S.yticks_g=S.ysg(1):dtick:S.yeg(1);

% Ticklabels (local):
if round(dtick)==dtick % if integer
    for i=1:length(S.xticks_g)
        S.xticklabels_g{i}=sprintf('%d',S.xticks_g(i));
    end
    for i=1:length(S.yticks_g)
        S.yticklabels_g{i}=sprintf('%d',S.yticks_g(i));
    end
else % if float
    for i=1:length(S.xticks_g)
        S.xticklabels_g{i}=sprintf('%8.1f',S.xticks_g(i));
    end
    for i=1:length(S.yticks_g)
        S.yticklabels_g{i}=sprintf('%8.1f',S.yticks_g(i));
    end
end


% get screensize
S.siz = get( 0, 'Screensize' );
S.fh.Position=S.siz;  % default: screensize


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
            'name','Timeslices',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui);
        S.fh.MenuBar = 'figure';
        
        
        S.ax = axes('unit','pix',...
            'position',[150 150 S.siz(3)-200 S.siz(4)-200]);
        % make initial plot
        S.Plot = mesh(S.xlocal,S.ylocal,S.tsl{1},'Linewidth',2);
        view(0,90)
        colormap(flipud(gray));
        colorbar
        xlabel('x [m]')
        ylabel('y [m]')
        set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
        set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
        set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
        axis xy
        set(S.ax,'DataAspectratio',[1 1 1])
        title(make_title(S,1))
        
        % UI control - Slider
        if length(tsl)>1
            S.sl = uicontrol('style','slide',...
                'unit','pix',...
                'position',[30 20 200 30],...
                'min',1,'max',length(tsl),'val',1,...
                'sliderstep',[1/(length(tsl)-1) 1/(length(tsl)-1)],'callback',@sl_callback);
        else
            S.sl = uicontrol('style','slide',...
                'unit','pix',...
                'position',[30 20 200 30],...
                'min',1,'max',1,'val',1,...
                'sliderstep',[1 1],'Enable','off');
        end
        S.text = uicontrol('style','text',...
            'unit','pix',...
            'position',[30 60 200 15],'String','Move slider to choose timeslice');
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[250 60 150 20],'Value',1,'Callback',@rb_call); % this button is on
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[250 40 150 20],'Callback',@rb_call);
        S.rb3=uicontrol('Style','radiobutton','String','2 %','Position',[250 20 150 20],'Callback',@rb_call);
        S.flag=1;
        
        % UI control - pushbutton for figure saving
        S.pb=uicontrol('Style','pushbutton','String','Save current timeslice as georeferenced PNG','Position',[400 35 280 25],'Callback',@save_tsl);
        S.pb_all=uicontrol('Style','pushbutton','String','Save all timeslice as georeferenced PNG','Position',[400 10 280 25],'Callback',@save_alltsl);
        S.text2 = uicontrol('style','text',...
            'unit','pix',...
            'position',[750 20 800 20],'String','');
        
        
        % UI control - checkbox for local/UTM coordinates
        S.cb=uicontrol('Style','checkbox','String','Display area in world coordinates (e.g. UTM)','Position',[400 60 300 30],'Callback',@cb_call);
        if S.ct==0
            S.cb.Enable='off';
        end
        
        
        % UI control - checkbox for topography
        S.cb_topo=uicontrol('Style','checkbox','String','Display topography','Position',[700 60 300 30],'Callback',@topo_call);
        
    end

% save handles
guidata(S.fh,S);



%%% Callback functions:

    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        wid=wid_fig-200; % width of axes
        hei_fig=S.fh.Position(4); % figure height
        hei=hei_fig-200;    % height of axes
        
        % change size of axes
        S.ax.Position=[150 150 wid hei];
    end

    function [] = sl_callback(varargin)
        S=guidata(gcbf);  % get guidata
        numtsl = round(S.sl.Value); % number of tsl
        set(S.Plot,'zdata',S.tsl{numtsl}); % Update timeslices plot
        drawnow;
        if S.cb.Value==0
            xlabel('x [m]')
            ylabel('y [m]')
            set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
            set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
            set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
        else
            xlabel('Easting [m]')
            ylabel('Northing [m]')
            set(S.ax,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
            set(S.ax,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
            set(S.ax,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
        end

        S.Plot.Parent.Title.String=make_title(S,numtsl);

        guidata(gcbf,S);  % Update guidata
        % Colorscale
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1
            set(S.ax,'ClimMode','auto');
            drawnow;
        elseif val2==1
            coldata=S.coldata{numtsl};
            cmin=coldata(ceil(length(coldata)/100*1));
            cmax=coldata(end-round(length(coldata)/100*1));
            set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
        elseif val3==1
            coldata=S.coldata{numtsl};
            cmin=coldata(ceil(length(coldata)/100*2));
            cmax=coldata(end-round(length(coldata)/100*2));
            set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
        end
    end

    function [] = rb_call(varargin)
        % Callback for radiobuttons
        S=guidata(gcbf);
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        numtsl=round(S.sl.Value); % number of current Tsl
        if val1==1 && S.flag~=1
            S.rb2.Value=0;
            S.rb3.Value=0;
            set(S.ax,'ClimMode','auto');
            drawnow;
            S.flag=1;
        elseif val2==1 && S.flag~=2
            S.rb1.Value=0;
            S.rb3.Value=0;
            coldata=S.coldata{numtsl};
            cmin=coldata(ceil(length(coldata)/100*1));
            cmax=coldata(end-round(length(coldata)/100*1));
            set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            S.flag=2;
        elseif val3==1 && S.flag~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            coldata=S.coldata{numtsl};
            cmin=coldata(ceil(length(coldata)/100*2));
            cmax=coldata(end-round(length(coldata)/100*2));
            set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            drawnow;
            S.flag=3;
        end
        guidata(gcbf,S); % Update (for flag)
    end

    function [] = cb_call(varargin)
        % callback for checkbox for coordinate transformation
        S=guidata(gcbf);
        state=S.cb.Value;   % on or off
        if state==0 % local coordinates of rotated area
            % check in which coordinate system currently:
            xm=mean(S.x(:));
            xm_l=mean(S.coordtrans(:,1)); % mean of local coords
            xm_w=mean(S.coordtrans(:,3)); % mean of world coords
            
            % if currently in world coordinates -> change
            if abs(xm-xm_w)<abs(xm-xm_l)
                % set new coordinates
                S.x=S.xlocal;
                S.y=S.ylocal;
            end
        elseif state==1 % area in world coordinates
            % check in which coordinate system currently:
            xm=mean(S.x(:));
            xm_l=mean(S.coordtrans(:,1)); % mean of local coords
            xm_w=mean(S.coordtrans(:,3)); % mean of world coords
            
            % if currently in local coordinates -> change
            if abs(xm-xm_w)>abs(xm-xm_l)
                % set new coordinates
                S.x=S.xglobal;
                S.y=S.yglobal;
            end
        end
        % update plot
        if S.cb_topo.Value==0 % plot tsl
            if S.flag==1
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                set(S.ax,'ClimMode','auto');
            elseif S.flag==2
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                coldata=S.coldata{round(S.sl.Value)};
                cmin=coldata(ceil(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            elseif S.flag==3
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                coldata=S.coldata{round(S.sl.Value)};
                cmin=coldata(ceil(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
            end
     
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(S.ax,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(S.ax,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            axis xy
            set(S.ax,'DataAspectratio',[1 1 1])
            colorbar
            colormap(flipud(gray))
            numtsl = round(S.sl.Value);
            title(make_title(S,numtsl));

            view(0,90)
            drawnow;
        elseif S.cb_topo.Value==1
            % plot topo alone
            S.Plot = mesh(S.x,S.y,S.topo,'Linewidth',2);
            colormap(jet)
            colorbar
            
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(S.ax,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(S.ax,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            axis xy
            set(S.ax,'DataAspectratio',[1 1 1])
            title('Topography')
            view(0,90)
            drawnow;
        end
        guidata(gcbf,S); % Update gui data
    end

    function [] = topo_call(varargin)
        % Callback for plotting topo
        S=guidata(gcbf);
        if S.cb_topo.Value==1 % plot topo
            set(S.sl,'Enable','off');
            set(S.rb1,'Enable','off');
            set(S.rb2,'Enable','off');
            set(S.rb3,'Enable','off');
            
            S.Plot = mesh(S.x,S.y,S.topo,'Linewidth',2);
            colormap(jet)
            colorbar
            
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(S.ax,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(S.ax,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            axis xy
            set(S.ax,'DataAspectratio',[1 1 1])
            title('Topography')
            view(0,90)
            drawnow;
        else % plot tsl
            set(S.sl,'Enable','on');
            set(S.rb1,'Enable','on');
            set(S.rb2,'Enable','on');
            set(S.rb3,'Enable','on');
            
            if S.flag==1
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                set(S.ax,'ClimMode','auto');
                view(0,90)
            elseif S.flag==2
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                coldata=S.coldata{round(S.sl.Value)};
                cmin=coldata(ceil(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
                view(0,90)
            elseif S.flag==3
                S.Plot = mesh(S.x,S.y,S.tsl{round(S.sl.Value)},'Linewidth',2);
                coldata=S.coldata{round(S.sl.Value)};
                cmin=coldata(ceil(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                set(S.ax,'ClimMode','manual','CLim',[cmin cmax]);
                view(0,90)
            end
            
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(S.ax,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(S.ax,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(S.ax,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(S.ax,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(S.ax,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            axis xy
            set(S.ax,'DataAspectratio',[1 1 1])
            colorbar
            colormap(flipud(gray))
            numtsl = round(S.sl.Value);
            title(make_title(S,numtsl))

            drawnow;
        end
        
        guidata(gcbf,S); % Update gui data
    end



    function [] = save_tsl(varargin)
        % Callback for pushbutton saving of current timeslice or topo
        S=guidata(gcbf);
        numtsl=round(S.sl.Value);   % number of timeslice
        set(gcf,'pointer','watch')
        % write png of tsl or topo
        if S.cb_topo.Value==0 %  tsl
            %%% Tsl with axes as png:
            fig=figure('Visible','off','PaperPosition',[0 0 50 50]);
            plo=mesh(S.x,S.y,S.tsl{numtsl},'Linewidth',2);
            view(0,90)
            colormap(flipud(gray));
            set(gca,'DataAspectratio',[1 1 1],'FontSize',20)
            axf=gca;     
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(axf,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(axf,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(axf,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(axf,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(axf,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(axf,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            axis xy
            title(make_title(S,numtsl));

            val1=S.rb1.Value;
            val2=S.rb2.Value;
            val3=S.rb3.Value;
            if val1==1
                set(axf,'ClimMode','auto');
            elseif val2==1
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                set(axf,'ClimMode','manual','CLim',[cmin cmax]);
            elseif val3==1
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                set(axf,'ClimMode','manual','CLim',[cmin cmax]);
            end
            saveas(fig,fullfile(S.pfad,make_fname(S,numtsl,'.png',1)));

            close(fig);
            
            %%% Georeferenced png:
            cdata=get(S.Plot,'zdata');
            cmin=min(cdata(:));
            cmax=max(cdata(:));
            tslname = fullfile(S.pfad,make_fname(S,numtsl,'.png',0));

            if S.flag==1
                cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
                cdata(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(cdata).*256,flipud(gray(256)),tslname,'Transparency',0);
            elseif S.flag==2
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                range=cmax-cmin;
                cdata=(cdata-cmin)/range;
                cdata(cdata<=0)=0;
                cdata(cdata>=1)=1;
                cdata(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(cdata).*256,flipud(gray(256)),tslname,'Transparency',0);
            elseif S.flag==3
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                range=cmax-cmin;
                cdata=(cdata-cmin)/range;
                cdata(cdata<=0)=0;
                cdata(cdata>=1)=1;
                cdata(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(cdata).*256,flipud(gray(256)),tslname,'Transparency',0);
            end
            % write pngw
            fname = make_fname(S,numtsl,'.pgw',0);
            if S.cb.Value==0    % local
                fid=fopen(fullfile(S.pfad,fname),'wt');
                fprintf(fid,[num2str(S.dx),'\n0\n0\n',num2str(-S.dx),'\n',num2str(min(S.xlocal(:))),'\n',num2str(max(S.ylocal(:)))]);
                fclose(fid);
            else % global
                write_geoPNGW(S.xlocal,S.ylocal,S.coordtrans,fullfile(S.pfad,fname));
            end
            set(S.text2,'String',['Timeslice saved successfully as ',make_fname(S,numtsl,'.png/pgw',0),'!'],'ForegroundColor','r','HorizontalAlignment','left');

            pause(2);
            set(S.text2,'String','');
        elseif S.cb_topo.Value==1 %  topo
            %%% Topo with axes as png:
            fig=figure('Visible','off','PaperPosition',[0 0 50 50]);
            mesh(S.x,S.y,S.topo,'Linewidth',2);
            view(0,90)
            colormap(jet);
            colorbar
            set(gca,'DataAspectratio',[1 1 1],'FontSize',20)
            
            axf=gca;
            
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(axf,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(axf,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(axf,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(axf,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(axf,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(axf,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            
            axis xy
            title('Topography')
            if S.cb.Value==0
                saveas(fig,fullfile(S.pfad,['TopoFigure_local.png']));
            else
                saveas(fig,fullfile(S.pfad,['TopoFigure_global.png']));
            end
            close(fig);
            
            %%% georeferenced PNG:
            cdata=get(S.Plot,'zdata');
            cmin=min(cdata(:));
            cmax=max(cdata(:));
            cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
            cdata(isnan(cdata))=0;  % set nan to 0
            m=ones(size(cdata));
            m(isnan(cdata))=0;
            im=cdata.*256;
            im(im<=2)=2;
            im(isnan(cdata))=0;  % set nan to 0
            % write png/pngw
            if S.cb.Value==0    % local
                imwrite(flipud(im),jet(256),fullfile(S.pfad,['Topography_local_',num2str(cmin,4),'-',num2str(cmax,4),'m.png']),'Transparency',0);
                
                fid=fopen(fullfile(S.pfad,['Topography_local_',num2str(cmin,4),'-',num2str(cmax,4),'m.pgw']),'wt');
                fprintf(fid,[num2str(S.dx),'\n0\n0\n',num2str(-S.dx),'\n',num2str(min(S.xlocal(:))),'\n',num2str(max(S.ylocal(:)))]);
                fclose(fid);
                set(S.text2,'String',['Topography saved successfully as ',['Topography_local_',num2str(cmin,4),'-',num2str(cmax,4),'m.png/pgw'],'!'],'ForegroundColor','r','HorizontalAlignment','left');
            else  % global
                imwrite(flipud(im),jet(256),fullfile(S.pfad,['Topography_global_',num2str(cmin,4),'-',num2str(cmax,4),'m.png']),'Transparency',0);
                
                write_geoPNGW(S.xlocal,S.ylocal,S.coordtrans,fullfile(S.pfad,['Topography_global_',num2str(cmin,4),'-',num2str(cmax,4),'m.pgw']));
                set(S.text2,'String',['Topography saved successfully as ',['Topography_global_',num2str(cmin,4),'-',num2str(cmax,4),'m.png/pgw'],'!'],'ForegroundColor','r','HorizontalAlignment','left');
            end
            pause(2);
            set(S.text2,'String','');
        end
        set(gcf,'pointer','arrow')
    end

    function [] = save_alltsl(varargin)
        % Callback for pushbutton saving of ALL timeslices (no topo!)
        S=guidata(gcbf);
        set(gcf,'pointer','watch')
        for numtsl=1:S.sl.Max
            %%% Tsl with axes as png:
            fig=figure('Visible','off','PaperPosition',[0 0 50 50]);
            plo=mesh(S.x,S.y,S.tsl{numtsl},'Linewidth',2);
            view(0,90)
            colormap(flipud(gray));
            set(gca,'DataAspectratio',[1 1 1],'FontSize',20)
       
            axf=gca;
            if S.cb.Value==0
                xlabel('x [m]')
                ylabel('y [m]')
                set(axf,'XTick',S.xticks_l,'YTick',S.yticks_l) % set ticks
                set(axf,'XTickLabel',S.xticklabels_l,'XTickLabelRotation',0)
                set(axf,'YTickLabel',S.yticklabels_l,'XLim',[S.xs(1) S.xe(1)],'YLim',[S.ys(1) S.ye(1)])
            else
                xlabel('Easting [m]')
                ylabel('Northing [m]')
                set(axf,'XTick',S.xticks_g,'YTick',S.yticks_g) % set ticks
                set(axf,'XTickLabel',S.xticklabels_g,'XTickLabelRotation',90)
                set(axf,'YTickLabel',S.yticklabels_g,'XLim',[S.xsg(1) S.xeg(1)],'YLim',[S.ysg(1) S.yeg(1)])
            end
            
            axis xy
            title(make_title(S,numtsl));
            val1=S.rb1.Value;
            val2=S.rb2.Value;
            val3=S.rb3.Value;
            if val1==1
                set(axf,'ClimMode','auto');
            elseif val2==1
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                set(axf,'ClimMode','manual','CLim',[cmin cmax]);
            elseif val3==1
                coldata=S.coldata{numtsl};
                cmin=coldata(ceil(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                set(axf,'ClimMode','manual','CLim',[cmin cmax]);
            end


            saveas(fig,fullfile(S.pfad,make_fname(S,numtsl,'.png',1)));
            
            close(fig);
            
            %%% Georeferenced png:
            cdata=S.tsl{numtsl};
            cmin=min(cdata(:));
            cmax=max(cdata(:));


            tslname = fullfile(S.pfad,make_fname(S,numtsl,'.png',0));
              

            if S.flag==1
                cdata=(cdata-cmin)./(cmax-cmin); % scale to 0-1
                m=ones(size(cdata));
                m(isnan(cdata))=0;
                im=cdata.*256;
                im(im<=2)=2;
                im(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(im),flipud(gray(256)),tslname,'Transparency',0);
            elseif S.flag==2
                coldata=S.coldata{numtsl};
                cmin=coldata(round(length(coldata)/100*1));
                cmax=coldata(end-round(length(coldata)/100*1));
                range=cmax-cmin;
                cdata=(cdata-cmin)/range;
                cdata(cdata<=0)=0;
                cdata(cdata>=1)=1;
                m=ones(size(cdata));
                m(isnan(cdata))=0;
                im=cdata.*256;
                im(im<=2)=2;
                im(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(im),flipud(gray(256)),tslname,'Transparency',0);
            elseif S.flag==3
                coldata=S.coldata{numtsl};
                cmin=coldata(round(length(coldata)/100*2));
                cmax=coldata(end-round(length(coldata)/100*2));
                range=cmax-cmin;
                cdata=(cdata-cmin)/range;
                cdata(cdata<=0)=0;
                cdata(cdata>=1)=1;
                m=ones(size(cdata));
                m(isnan(cdata))=0;
                im=cdata.*256;
                im(im<=2)=2;
                im(isnan(cdata))=0;  % set nan to 0
                imwrite(flipud(im),flipud(gray(256)),tslname,'Transparency',0);
            end

            % write pngw
            fname = make_fname(S,numtsl,'.pgw',0);
            if S.cb.Value==0    % local
                fid=fopen(fullfile(S.pfad,fname),'wt');
                fprintf(fid,[num2str(S.dx),'\n0\n0\n',num2str(-S.dx),'\n',num2str(min(S.xlocal(:))),'\n',num2str(max(S.ylocal(:)))]);
                fclose(fid);
            else % global
                write_geoPNGW(S.xlocal,S.ylocal,S.coordtrans,fullfile(S.pfad,fname));
            end

        end
        
        set(S.text2,'String',['All timeslices saved successfully!'],'ForegroundColor','r','HorizontalAlignment','left');
        pause(2);
        set(S.text2,'String','');
        set(gcf,'pointer','arrow')
    end
    
    function titleStr = make_title(S,numtsl)
        if S.dsl
            dslStr = [' (',num2str(S.maxElevation - S.t(numtsl,2),'%5.2f'), ' - ',num2str(S.maxElevation - S.t(numtsl,1),'%5.2f'),' m asl.)'];
            dslStr2='Depth ';
        else
            dslStr = [];
            dslStr2=[];
        end
        titleStr = [dslStr2,num2str(S.t(numtsl,1),2),' - ',num2str(S.t(numtsl,2),2),S.unit,dslStr];
    end
    
    function fnameStr = make_fname(S,numtsl,extension,fig)
        
        % is this file a Matlab-Figure or a Tsl/Dsl Image
        if fig
            fig = 'Figure';
        else
            fig = [];
        end
        
        % are the coordinates local or global
        if S.cb.Value
            cb = '_global_';
        else
            cb = '_local_';
        end
        
        % is it a depth slice or a time slice
        if S.dsl
            dslStr = ['_-_',num2str(S.maxElevation - S.t(numtsl,2),'%5.2f'), '-',...
                    num2str(S.maxElevation - S.t(numtsl,1),'%5.2f'),'m'];
        else
            dslStr = [];
        end

        % create filename
        fnameStr = ['Tsl',fig,'_',num2str(numtsl,'%2d'),cb,num2str(S.t(numtsl,1),2),...
            '-',num2str(S.t(numtsl,2),2),S.unit,dslStr,extension];
    end
end
