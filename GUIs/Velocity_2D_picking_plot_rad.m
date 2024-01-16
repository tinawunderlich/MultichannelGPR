function [] = Velocity_2D_picking_plot_rad(foldername,name,profilelist,chanlist)

% function [] = Velocity_2D_picking_plot_rad(foldername,name,profilelist,chanlist)
% FOR RADARGRAMS.MAT and correspondning files!
%
% Dr. Tina Wunderlich, CAU Kiel 2019, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% foldername: complete foldername and path of folder containing
% radargrams.mat,...
% name: Name of datafiles without '_...'
% profilelist: numbers of profiles
% chanlist: number of channels
%
% requires folders Export_Import, Processing and Subfunctions




% set data
S.foldername=foldername; % name of rSlicer-folder
S.name=name; % name of data
S.profiles=profilelist'; % list of available profiles
S.chanlist=chanlist'; % list of channels
S.currchan=1;   % current channel number
S.hyplen=0.5;   % half width of hyperbola
S.profnum=S.profiles(1);
S.liste=[];
for i=1:length(profilelist)
    for j=1:length(chanlist)
        S.liste=[S.liste; profilelist(i) chanlist(j)];
    end
end

% load data
S.radargrams=load(fullfile(foldername,'radargrams.mat'));
S.radargrams=S.radargrams.radargrams;
x=load(fullfile(foldername,'x.mat'));
S.t=load(fullfile(foldername,'t.mat'));
S.t=S.t.t;
S.global_coords=load(fullfile(foldername,'global_coords.mat'));
S.global_coords=S.global_coords.global_coords;



% set first radargram
S.traces=S.radargrams{1};
S.dt=S.t(2)-S.t(1);
S.ns=length(S.t);
S.x=S.global_coords{1}(:,1);
S.y=S.global_coords{1}(:,2);
S.coord=x.x{1};    % profile coordinate
S.dx=mean(sqrt(diff(S.x).^2+diff(S.y).^2));   % mean dx

% calculate gained data
gainf=[-20 0 15 25 30];
g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
S.tracesgain=S.traces.*repmat((10.^(g./20))',[1 length(S.traces(1,:))]); % apply gain

% colorscale limits
coldata=sort(unique(S.tracesgain(~isnan(S.tracesgain(:)))));
S.cmingain1=coldata(round(length(coldata)/100*1));
S.cmaxgain1=coldata(end-round(length(coldata)/100*1));
S.cmingain3=coldata(round(length(coldata)/100*3));
S.cmaxgain3=coldata(end-round(length(coldata)/100*3));

coldata=sort(unique(S.traces(~isnan(S.traces(:)))));
S.cmin1=coldata(round(length(coldata)/100*1));
S.cmax1=coldata(end-round(length(coldata)/100*1));
S.cmin3=coldata(round(length(coldata)/100*3));
S.cmax3=coldata(end-round(length(coldata)/100*3));

% get screensize
S.siz=get(0,'ScreenSize');
S.fh.Position=S.siz;  % default: screensize

% Add the UI components
S=addcomponents(S);

% Make figure visible after adding components
S.fh.Visible='on';

% Toolbar for zooming + pan
figureToolBar = uimenu('Label','Zoom');
uimenu(figureToolBar,'Label','Zoom In','Callback','zoom on');
uimenu(figureToolBar,'Label','Zoom Out','Callback','zoom out');
uimenu(figureToolBar,'Label','Pan','Callback','pan on');


    function S=addcomponents(S)
        
        S.fh= figure('menubar','none','Position',S.siz,...
            'name','Velocity picking',...
            'numbertitle','off',...
            'resize','on','Visible','off','SizeChangedFcn',@resizeui);
        
        %%% UI control elements:
        % UI control - listbox for channel number
        S.chan = uicontrol(S.fh,'style','listbox','unit','pix','position',[20 220 50 60],'callback',{@chan_call},'Value',1,'String',int2str(S.chanlist));
        S.text = uicontrol(S.fh,'style','text','unit','pix','position',[20 290 100 15],'String','Choose channel','HorizontalAlignment','Left');
        
        % UI control for profile number
        S.prof=uicontrol(S.fh,'Style','listbox','unit','pix','position',[20 310 50 150],'String',int2str(S.profiles),'callback',{@profnum_call});
        S.proftext=uicontrol(S.fh,'style','text','unit','pix','position',[20 470 100 15],'String','Choose profile','HorizontalAlignment','left');
        
        % UI control aspectratio
        S.asp_String=[{'1/1'} {'1/5'} {'1/10'} {'1/20'} {'1/30'} ];
        S.aspval=1/20;
        S.asp = uicontrol(S.fh,'style','listbox','unit','pix','position',[20 20 50 60],'String',S.asp_String,'Value',4,'Callback',@asp_call);
        S.asptext=uicontrol(S.fh,'style','text','unit','pix','position',[20 85 100 20],'String','Aspect ratio t','HorizontalAlignment','left');
        
        % UI control - radiobuttons for colorscale
        S.rb1=uicontrol('Style','radiobutton','String','Auto color scale','Position',[280 60 150 15],'Value',1,'Callback',@colorscale_call); % this button is on
        S.rb2=uicontrol('Style','radiobutton','String','1 %','Position',[280 40 150 15],'Callback',@colorscale_call);
        S.rb3=uicontrol('Style','radiobutton','String','3 %','Position',[280 20 150 15],'Callback',@colorscale_call);
        S.flag_cs=1;
        
        % UI control: gain
        S.gain_cb=uicontrol(S.fh,'Style','checkbox','unit','pix','Position',[150 120 100 15],'String','Apply gain [dB]','Value',1,'Callback',@gain_call);
        S.gain1 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[150 100 60 15],'String','-20','Callback',@gain_call);
        S.gain2 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[150 80 60 15],'String','0','Callback',@gain_call);
        S.gain3 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[150 60 60 15],'String','15','Callback',@gain_call);
        S.gain4 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[150 40 60 15],'String','25','Callback',@gain_call);
        S.gain5 = uicontrol(S.fh,'Style','edit','unit','pix','Position',[150 20 60 15],'String','30','Callback',@gain_call); 
        
        % UI control for new v-pick
        S.new=uicontrol(S.fh,'Style','pushbutton','unit','pix','position',[480 55 100 20],'callback',{@new_call},'String','Make new v-pick');
        
        % UI control for semblance button
        S.semb=uicontrol(S.fh,'Style','checkbox','unit','pix','position',[480 20 150 20],'callback',{@semb_call},'Value',0,'String','Auto v (via semblance)');
        
        % UI control for v-slider
        S.vSlider=uicontrol(S.fh,'Style','slider','unit','pix','position',[650 30 100 20],'callback',{@v_call},'min',3,'max',18,'Value',10,'sliderstep',[0.1/15 1/15],'Enable','off');
        S.text6=uicontrol(S.fh,'style','text','unit','pix','position',[650 55 100 20],'String','v =           cm/ns');
        S.vtext=uicontrol(S.fh,'style','text','unit','pix','position',[675 55 30 20],'String','10');
        
        % UI control for Hyperbola length
        S.hyplenedit=uicontrol(S.fh,'Style','edit','unit','pix','position',[800 20 40 20],'callback',{@hyplen_call},'Value',S.hyplen,'String',num2str(S.hyplen));
        S.text7=uicontrol(S.fh,'style','text','unit','pix','position',[850 20 160 15],'String','Half width of hyperbola [m]','HorizontalAlignment','Left');
        
        % UI control for saving of velocities
        S.save=uicontrol(S.fh,'Style','pushbutton','unit','pix','position',[800 55 120 20],'callback',{@save_call},'String','Save current velocity','Enable','off');
        
        % UI control for tmax
        S.tmax=uicontrol(S.fh,'Style','edit','unit','pix','position',[20 140 50 20],'callback',{@tmax_call},'Value',max(S.t),'String',num2str(max(S.t)));
        S.text8=uicontrol(S.fh,'style','text','unit','pix','position',[20 160 80 20],'String','tmax [ns]','HorizontalAlignment','Left');
        
        % Plot of radargram
        S.ax = axes('unit','pix',...
            'position',[200 140 S.siz(3)-300 S.siz(4)-200]);
        S.Plot = imagesc(S.coord,S.t,S.tracesgain);
        hold on
        xlabel('x [m]')
        ylabel('t [ns]')
        grid on
        title(['Profile ',int2str(S.profiles(1)),', Channel ',int2str(S.currchan)])
        colormap(flipud(gray))
        set(S.ax,'DataAspectratio',[S.aspval 1 1])
        
    end

% save handles
guidata(S.fh,S);

%--------------------------------------------------------------------------
%%% Callback functions:

    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        wid=wid_fig-300; % width of axes
        hei_fig=S.fh.Position(4); % figure height
        hei=hei_fig-200;    % height of axes
        
        % change size of axes
        S.ax.Position=[200 140 wid hei];
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


    function [] = save_call(varargin)
        % callback for hyperbola length
        S=guidata(gcbf);
        if isfield(S,'v')
            v=S.v;
            xloc=S.xpick;
            tloc=S.tpick;
            ind=find(min(abs(S.coord-xloc))==abs(S.coord-xloc));
            xy=[S.x(ind(1)) S.y(ind(1))];
            S.profnum=S.profiles(S.prof.Value);
            if exist(fullfile(S.foldername,'Velocity_picks.txt'))
                fid=fopen(fullfile(S.foldername,'Velocity_picks.txt'),'a+');
                fprintf(fid,'%d\t%d\t%10.3f\t%10.3f\t%4.2f\t%4.2f\n',[S.profnum S.currchan xy tloc v]');
                fclose(fid);
            else
                fid=fopen(fullfile(S.foldername,'Velocity_picks.txt'),'wt');
                fprintf(fid,'Profilenumber\tChannel\tx [m]\ty [m]\tTWT [ns]\tv [cm/ns]\n');
                fprintf(fid,'%d\t%d\t%10.3f\t%10.3f\t%4.2f\t%4.2f\n',[S.profnum S.currchan xy tloc v]');
                fclose(fid);
            end
            set(S.hyp,'Color','y');
        end
    end

    function [] = hyplen_call(varargin)
        % callback for hyperbola length
        S=guidata(gcbf);
        S.hyplen=str2num(S.hyplenedit.String);
        % update handles
        guidata(gcf,S);
        semb_call(varargin);
    end

    function [] = v_call(varargin)
        % callback for v-slider
        S=guidata(gcbf);
        S.v=S.vSlider.Value;
        S.vtext.String=num2str(S.v,3);
        % calculate x and t along hyperbola
        xhyp=S.coord(abs(S.coord-S.xpick)<=S.hyplen);
        thyp=2.*sqrt((S.tpick./2).^2+((xhyp-S.xpick)./(S.v/100)).^2);
        % update plot hyperbola
        set(S.hyp,'XData',xhyp,'YData',thyp);
        drawnow;
        % update handles
        guidata(gcf,S);
    end

    function [] = semb_call(varargin)
        % callback for semblance checkbox
        S=guidata(gcbf);
        if S.semb.Value==1 % automatically determine velocity
            set(S.vSlider,'Enable','off');
            if isfield(S,'xpick')
                S=calcsemblance(S); % result is in S.v
                set(S.vtext,'String',num2str(S.v,3));
            end
        else  % manually fit hyperbola
            set(S.vSlider,'Enable','on');
            % get v
            S.v=S.vSlider.Value; % in cm/ns
            set(S.vtext,'String',num2str(S.v,3));
        end
        if isfield(S,'xpick')
            % calculate x and t along hyperbola
            xhyp=S.coord(abs(S.coord-S.xpick)<=S.hyplen);
            thyp=2.*sqrt((S.tpick./2).^2+((xhyp-S.xpick)./(S.v/100)).^2);
            % plot hyperbola
            set(S.hyp,'XData',xhyp,'YData',thyp);
            drawnow;
        end
        % update handles
        guidata(gcf,S);
    end

    function [] = gain_call(varargin)
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'watch');
        S=guidata(gcbf);
        % calculate gain new
        gainf=[str2num(S.gain1.String) str2num(S.gain2.String) str2num(S.gain3.String) str2num(S.gain4.String) str2num(S.gain5.String)];
        g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
        S.tracesgain=S.traces.*repmat((10.^(g./20))',[1 length(S.traces(1,:))]); % apply gain
        guidata(gcbf,S); % Update
        newprofile_call(varargin);
    end


    function [] = new_call(varargin)
        % Callback for new v-pick-button
        S=guidata(gcbf);
        % activate pick of apex
        [S.xpick,S.tpick]=ginput(1);
        if S.semb.Value==1 % automatically determine velocity
            S=calcsemblance(S); % result is in S.v
            set(S.vtext,'String',num2str(S.v,3));
        else  % manually fit hyperbola
            set(S.vSlider,'Enable','on');
            % get v
            S.v=S.vSlider.Value; % in cm/ns
        end
        % calculate x and t along hyperbola
        xhyp=S.coord(abs(S.coord-S.xpick)<=S.hyplen);
        thyp=2.*sqrt((S.tpick./2).^2+((xhyp-S.xpick)./(S.v/100)).^2);
        % plot hyperbola
        S.hyp=plot(xhyp,thyp,'r','linewidth',2);
        drawnow;
        set(S.save,'Enable','On');
        % update handles
        guidata(gcf,S);
    end


    function S = calcsemblance(S)
        % calculate semblance
        xx=S.coord(abs(S.coord-S.xpick)<=S.hyplen); % x-values
        % change pointer
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'watch')
        drawnow;
        v=0.03:0.001:0.18; % in cm/ns
        M=15; % window width for semblance calculation
        for i=1:length(v)   % compute for every velocity
            tt=2.*sqrt((S.tpick./2).^2+((xx-S.xpick)./v(i)).^2);    % t values along hyperbola
            
            % check if sample is in aperture of antenna
            beta=atand((xx-S.xpick)./(tt/2*v(i)));
            in=(abs(beta)<=40/2);
            xx=xx(in);  % x-values in aperture
            temp=ismember(S.coord,xx);
            xx2{i}=find(temp==1);  % only use samples in aperture
            tt2{i}=round(tt(in)/S.dt);
            % xx2 and tt2 are now indices!
            
            if length(xx2{i})>1
                % find samples around hyperbola
                for j=1:length(xx2{i})
                    if xx2{i}(j)>0 && tt2{i}(j)-M>=1 && tt2{i}(j)+M<=length(S.t)
                        hyp(:,j)=S.traces(tt2{i}(j)-M:tt2{i}(j)+M,xx2{i}(j));
                    else
                        hyp(:,j)=NaN(2*M+1,1);
                    end
                end
                hyp(:,isnan(hyp(1,:)))=[];
                
                % calculate semblance
                klammer_oben=sum(hyp,2);
                oben=sum(klammer_oben.^2);
                klammer_unten=sum(hyp.^2,2);
                unten=length(xx2{i})*sum(klammer_unten);
                semblance(i)=oben/unten;    % semblance for current velocity
            else
                semblance(i)=0;
            end
        end
        % change pointer back to arrow
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'arrow')
        % index of max semblance
        maxind=find(max(semblance)==semblance);
        S.v=v(maxind(1))*100;  % best velocity in cm/ns
    end


    function [] = tmax_call(varargin)
        % Callback for edit field tmax
        S=guidata(gcbf);
        tmax=str2num(S.tmax.String);
        set(S.ax,'YLim',[0 tmax]);
        drawnow;
    end


    function [] = profnum_call(varargin)
        % Callback for profile number listbox
        S=guidata(gcbf);
        % delete old hyperbolas:
        obj=get(gca,'Children');
        nobj=length(obj);
        for o=1:nobj-1
            delete(obj(o));
        end
        set(S.save,'Enable','Off');
        % current profile number
        S.profnum=S.profiles(S.prof.Value);
        % current channel number
        S.currchan=S.chanlist(S.chan.Value);
        % change pointer
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'watch')
        drawnow;
        % read new data
        S.traces=S.radargrams{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan};
        S.x=S.global_coords{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan}(:,1);
        S.y=S.global_coords{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan}(:,2);
        S.coord=cumsum([0; sqrt(diff(S.x).^2+diff(S.y).^2)]);    % profile coordinate
        S.dx=mean(sqrt(diff(S.x).^2+diff(S.y).^2));   % mean dx
        % calculate gain new
        gainf=[str2num(S.gain1.String) str2num(S.gain2.String) str2num(S.gain3.String) str2num(S.gain4.String) str2num(S.gain5.String)];
        g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
        S.tracesgain=S.traces.*repmat((10.^(g./20))',[1 length(S.traces(1,:))]); % apply gain
        % colorscale limits
        coldata=sort(unique(S.tracesgain(~isnan(S.tracesgain(:)))));
        S.cmingain1=coldata(round(length(coldata)/100*1));
        S.cmaxgain1=coldata(end-round(length(coldata)/100*1));
        S.cmingain3=coldata(round(length(coldata)/100*3));
        S.cmaxgain3=coldata(end-round(length(coldata)/100*3));
        
        coldata=sort(unique(S.traces(~isnan(S.traces(:)))));
        S.cmin1=coldata(round(length(coldata)/100*1));
        S.cmax1=coldata(end-round(length(coldata)/100*1));
        S.cmin3=coldata(round(length(coldata)/100*3));
        S.cmax3=coldata(end-round(length(coldata)/100*3));

        % update handles
        guidata(gcf,S);
        % chnage pointer back to arrow
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'arrow')
        newprofile_call(varargin); % call plotting function
    end


    function [] = asp_call(varargin)
        % Callback for aspectratio
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'watch')
        S=guidata(gcbf);
        S.aspval=S.asp_String(S.asp.Value);
        S.aspval=str2num(S.aspval{1});
        guidata(gcbf,S); % Update
        % call plot function
        newprofile_call(varargin);
    end


    function [] = chan_call(varargin)
        % call for plotting new channel
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'watch')
        S=guidata(gcbf);
        % delete old hyperbolas:
        obj=get(gca,'Children');
        nobj=length(obj);
        for o=1:nobj-1
            delete(obj(o));
        end
        S.currchan=S.chanlist(S.chan.Value);
        % load new channel data
        S.traces=S.radargrams{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan};
        S.x=S.global_coords{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan}(:,1);
        S.y=S.global_coords{S.liste(:,1)==S.profnum & S.liste(:,2)==S.currchan}(:,2);
        S.coord=cumsum([0; sqrt(diff(S.x).^2+diff(S.y).^2)]);    % profile coordinate
        S.dx=mean(sqrt(diff(S.x).^2+diff(S.y).^2));   % mean dx
        % calculate gain new
        gainf=[str2num(S.gain1.String) str2num(S.gain2.String) str2num(S.gain3.String) str2num(S.gain4.String) str2num(S.gain5.String)];
        g=interp1(linspace(0,max(S.t(:)),length(gainf)),gainf,S.t);
        S.tracesgain=S.traces.*repmat((10.^(g./20))',[1 length(S.traces(1,:))]); % apply gain
        % colorscale limits
        coldata=sort(unique(S.tracesgain(~isnan(S.tracesgain(:)))));
        S.cmingain1=coldata(round(length(coldata)/100*1));
        S.cmaxgain1=coldata(end-round(length(coldata)/100*1));
        S.cmingain3=coldata(round(length(coldata)/100*3));
        S.cmaxgain3=coldata(end-round(length(coldata)/100*3));
        
        coldata=sort(unique(S.traces(~isnan(S.traces(:)))));
        S.cmin1=coldata(round(length(coldata)/100*1));
        S.cmax1=coldata(end-round(length(coldata)/100*1));
        S.cmin3=coldata(round(length(coldata)/100*3));
        S.cmax3=coldata(end-round(length(coldata)/100*3));
        
        set(S.save,'Enable','Off');
        % update handles
        guidata(gcf,S);
        newprofile_call(varargin); % call plotting function
    end


    function [] = newprofile_call(varargin)
        % call for plotting new radargram
        S=guidata(gcbf);

        % plot radargram
        if S.gain_cb.Value==1
            set(S.Plot,'XData',S.coord,'CData',S.tracesgain);
        else
            set(S.Plot,'XData',S.coord,'CData',S.traces);
        end
        set(S.ax,'Xlim',[0 max(S.coord)],'DataAspectratio',[S.aspval 1 1]);
        title(['Profile ',int2str(S.profnum),', Channel ',int2str(S.currchan)])
        val1=S.rb1.Value;
        val2=S.rb2.Value;
        val3=S.rb3.Value;
        if val1==1 && S.flag_cs~=1
            S.rb2.Value=0;
            S.rb3.Value=0;
            set(S.ax,'ClimMode','auto');
            S.flag_cs=1;
        elseif val2==1 && S.flag_cs~=2
            S.rb1.Value=0;
            S.rb3.Value=0;
            if S.gain_cb.Value==1
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain1 S.cmaxgain1]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmin1 S.cmax1]);
            end
            S.flag_cs=2;
        elseif val3==1 && S.flag_cs~=3
            S.rb1.Value=0;
            S.rb2.Value=0;
            if S.gain_cb.Value==1
                set(S.ax,'ClimMode','manual','CLim',[S.cmingain3 S.cmaxgain3]);
            else
                set(S.ax,'ClimMode','manual','CLim',[S.cmin3 S.cmax3]);
            end
            S.flag_cs=3;
        end
        drawnow;
        
        % chnage pointer back to arrow
        set(findobj('Type','Figure','Name','Velocity picking'), 'pointer', 'arrow')
    end
end