function [] = Check_t0_slider_plot(foldername,name,profilelist)

% function [] = Check_t0_slider_plot(foldername,name,profilelist)
%
% GUI for plotting of script Check_t0
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% foldername: complete rSlicer-foldername and path
% name: Name of datafiles without '_...'
% profilelist: numbers of profiles
%
% requires folders Export_Import, Subfunctions and Processing


% set data
S.threshhold=-2000;
S.t0=NaN;
S.t0pick=S.t0;
S.flag=2; % use threshhold criterium first
S.foldername=foldername;
S.name=name;
S.profiles=profilelist;

% load data of first profile
[S.traces,S.dt,S.ns]=readmala(S.foldername,S.name,S.profiles(1));
S.traceanz=length(S.traces{1}(1,:));
S.chan=length(S.traces);
S.t=0:S.dt:(S.ns-1)*S.dt;
% Process initial data
S.alltraces_proc=zeros(S.ns,S.traceanz*length(S.traces));
for i=1:length(S.traces)
    % Amplitude offset removal:
    S.traces_proc{i}=DCremoval(S.traces{i},S.t);
    % Bandpassfilter:
    S.traces_proc{i}=bandpass_gpr(S.traces_proc{i},S.dt,100,600);
    
    S.alltraces_proc(:,1+(i-1)*S.traceanz:i*S.traceanz) = S.traces_proc{i};
end
% S.traces_proc = processed traces (=DCremoval + bandpass)


% reference trace
S.reftrace=S.traces_proc{1}(:,1);
S.refchannel=1;
S.reftracenumber=1;
S.refprofile=S.profiles(1);



% t0 correction
[S.traces_t0,S.shiftsamples]=channelshift(S.traces_proc,40); % processed and t0-corrected traces

S.alltraces_t0=zeros(S.ns,S.traceanz*length(S.traces));
for i=1:length(S.traces)
    S.alltraces_t0(:,1+(i-1)*S.traceanz:i*S.traceanz)=S.traces_t0{i};
end

%--------------------------------------------------------------
% get screensize
S.siz = get( 0, 'Screensize' );
S.fh.Position=S.siz;  % default: screensize

% Add the UI components
S=addcomponents(S);

% Make figure visible after adding components
S.fh.Visible='on';

    function S=addcomponents(S)
        
        % create figure handle
        S.fh = figure('units','pixels',...
            'position',S.siz,...
            'menubar','none',...
            'name','Traces',...
            'numbertitle','off',...
            'SizeChangedFcn',@resizeui,'Visible','off');
        
        % Plot of reference trace and t0 (left top)
        S.ax3 = axes('unit','pix',...
            'position',[20 440 100 S.siz(4)-460]);
        S.Plot3a = plot(S.reftrace,S.t,'k','Linewidth',2);
        axis ij
        title('Ref. trace and t0')
        grid on
        set(S.ax3,'YLim',[0 10],'XTick',[0],'XTickLabel',[])

        
        % Plot of original data (middle)
        S.ax1 = axes('unit','pix',...
            'position',[320 420 S.siz(3)-320 240]);
        S.Plot1a= imagesc(1:length(S.alltraces_proc(1,:)),S.t,S.alltraces_proc);
        ylabel('t [ns]')
        colormap(flipud(gray))
        set(S.ax1,'YLim',[0 50],'XTick',[1:length(S.alltraces_proc(1,:))/length(S.traces):length(S.alltraces_proc(1,:))],'XTickLabel',1:length(S.traces))
        grid on
        title('Original data (incl. DCremoval + bandpass)')
        axis ij
        
        % Plot of t0 corrected data (bottom)
        S.ax2 = axes('unit','pix',...
            'position',[320 120 S.siz(3)-320 240]);
        S.Plot2a = imagesc(1:length(S.alltraces_t0(1,:)),S.t,S.alltraces_t0);
        set(S.ax2,'YLim',[0 50],'XTick',[1:length(S.alltraces_proc(1,:))/length(S.traces):length(S.alltraces_proc(1,:))],'XTickLabel',1:length(S.traces))
        ylabel('t [ns]')
        colormap(flipud(gray))
        grid on
        title('t0 corrected data (incl. DCremoval + bandpass before)')
        axis ij
        
        % wiggle plots:
        S.axw1=axes('unit','pix','position',[200 420 100 240]);
        S.Wiggle1=plot(S.traces_proc{1}(:,1),S.t);
        hold on
        for ii=2:length(S.traces)
            S.Wiggle1=plot(S.traces_proc{ii}(:,1),S.t);
        end
        axis ij
        set(S.axw1,'YLim',[0 50])
        
        S.axw2=axes('unit','pix','position',[200 120 100 240]);
        S.Wiggle2=plot(S.traces_t0{1}(:,1),S.t);
        hold on
        for ii=2:length(S.traces)
            S.Wiggle2=plot(S.traces_t0{ii}(:,1),S.t);
        end
        axis ij
        set(S.axw2,'YLim',[0 50])
        
        %%% UI control elements:
        % UI control for profile number
        S.prof=uicontrol('Style','listbox','unit','pix','position',[20 120 50 240],'String',int2str(S.profiles),'callback',{@profnum_call});
        S.proftext=uicontrol('style','text','unit','pix','position',[20 380 100 20],'String','Choose profile','HorizontalAlignment','left');
        
        % UI control edit field for threshhold
        S.thr=uicontrol('Style','edit','unit','pix','position',[300 30 80 20],'Enable','off','callback',{@thr_call},'Value',S.threshhold,'String',num2str(S.threshhold));
        S.text3=uicontrol('style','text','unit','pix','position',[300 60 100 20],'String','Choose threshold','HorizontalAlignment','Left');
        
        % UI control edit field for tmax
        S.tmax=uicontrol('Style','edit','unit','pix','position',[450 30 100 20],'callback',{@tmax_call},'Value',50,'String',num2str(50));
        S.text4=uicontrol('style','text','unit','pix','position',[450 60 200 20],'String','Choose tmax for plotting','HorizontalAlignment','Left');
        
        % UI control for popupmenu for different possibilities for t0 determination
        S.text5=uicontrol('style','text','unit','pix','position',[600 60 300 20],'String','Choose method for t0-correction','HorizontalAlignment','Left');
        S.pop=uicontrol('Style','popupmenu','unit','pix','position',[600 30 300 20],'callback',{@pop_call},'String',{'Use channelshift for each profile','Use channelshift for all profiles','Use t0-pick on reference trace and cross-correlation','Use amplitude threshold criterium'},'Value',1);
        
        % UI control for normalization
        S.normcheck=uicontrol('Style','checkbox','unit','pix','position',[950 30 200 20],'String','Normalize data before t0 correction?','Value',0,'callback',{@norm_call});
        S.norm=0;
        
        % UI control for reference trace
        S.refch=uicontrol('Style','edit','unit','pix','position',[20 (S.siz(4)-200)/2+180+20 20 20],'Enable','off','callback',{@ref_call},'Value',S.refchannel,'String',int2str(S.refchannel));
        S.refchtext=uicontrol('style','text','unit','pix','position',[40 (S.siz(4)-200)/2+180+20 100 20],'String','Channel','HorizontalAlignment','Left');
        
        S.refprof=uicontrol('Style','edit','unit','pix','position',[20 (S.siz(4)-200)/2+180 20 20],'Enable','off','callback',{@ref_call},'Value',S.refprofile,'String',int2str(S.refprofile));
        S.refproftext=uicontrol('style','text','unit','pix','position',[40 (S.siz(4)-200)/2+180 100 20],'String','Profile','HorizontalAlignment','Left');
        
        S.reftrnum=uicontrol('Style','edit','unit','pix','position',[20 (S.siz(4)-200)/2+180+40 20 20],'Enable','off','callback',{@ref_call},'Value',S.reftracenumber,'String',int2str(S.reftracenumber));
        S.reftrnumtext=uicontrol('style','text','unit','pix','position',[40 (S.siz(4)-200)/2+180+40 100 20],'String','Tracenumber','HorizontalAlignment','Left');
        
        S.t0field=uicontrol('Style','edit','unit','pix','position',[20 (S.siz(4)-200)/2+180+60 20 20],'Enable','off','callback',{@t0call},'Value',0.5,'String','0.5');
        S.t0text=uicontrol('style','text','unit','pix','position',[40 (S.siz(4)-200)/2+180+60 100 20],'String','t0 [ns]','HorizontalAlignment','Left');
        
        % UI control for saving of settings
        S.save=uicontrol('Style','pushbutton','String','Save settings','Enable','off','unit','pix','position',[20 30 200 20],'callback',{@save_call});
       
        % UI control for info  text
        S.info=uicontrol('Style','text','String',['Shiftsamples = [',int2str(S.shiftsamples'),'] for current profile'],'unit','pix','position',[600 5 600 20],'Foregroundcolor','r');
    end
% save handles
guidata(S.fh,S);

%--------------------------------------------------------------------------
%%% Callback functions:
    function resizeui(hObject,event)
        % get current size of figure
        wid_fig=S.fh.Position(3); % Figure width
        wid=wid_fig-220; % width of axes
        hei_fig=S.fh.Position(4); % figure height
        hei=hei_fig;    % height of axes
        
        % change size of axes
        set(S.ax3,'position',[20 (hei-200)/2+180+85 100 (hei-200)/2-100]);    % ref trace
        set(S.refch,'position',[20 (hei-200)/2+180+20 30 20]);
        set(S.refchtext,'position',[55 (hei-200)/2+180+20 100 20]);
        set(S.refprof,'position',[20 (hei-200)/2+180 30 20]);
        set(S.refproftext,'position',[55 (hei-200)/2+180 100 20]);
        set(S.reftrnum,'position',[20 (hei-200)/2+180+40 30 20]);
        set(S.reftrnumtext,'position',[55 (hei-200)/2+180+40 100 20]);
        set(S.t0field,'position',[20 (hei-200)/2+180+60 30 20]);
        set(S.t0text,'position',[55 (hei-200)/2+180+60 100 20]);
        
        set(S.ax1,'position',[320 (hei-200)/2+180 wid-120 (hei-200)/2]);    % original data
        set(S.ax2,'position',[320 130 wid-120 (hei-200)/2]);    % t0-corr data
        set(S.axw1,'position',[200 (hei-200)/2+180 100 (hei-200)/2]);    % wiggle original
        set(S.axw2,'position',[200 130 100 (hei-200)/2]);    % wiggle t0 corr
        set(S.prof,'position',[20 130 50 (hei-200)/2]);    % profile list
        set(S.proftext,'position',[20 (hei-200)/2+135 100 15]);    % profile list text
    end


    function [] = norm_call(varargin)
        % callback for normalization checkbox
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        S.norm=S.normcheck.Value;
        % update handles
        guidata(gcf,S);
        % call new plot
        plotnew_call(varargin);
    end

    function [] = t0call(varargin)
        % callback for t0 
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        S.t0pick=str2num(S.t0field.String);
        if S.t0pick<0
            S.t0pick=0;
        elseif S.t0pick>max(S.t)
            S.t0pick=max(S.t);
        end
        % update t0-line in reference plot
        if isfield(S,'Plot3b')
            set(S.Plot3b,'Ydata',[S.t0pick S.t0pick]);
        else
            hold(S.ax3,'on');
            S.Plot3b=plot(S.ax3,[min(S.reftrace) max(S.reftrace)],[S.t0pick S.t0pick],'r','Linewidth',2);
        end
        drawnow;
        % update handles
        guidata(gcf,S);
        % call new plot
        plotnew_call(varargin);
    end

    function [] = ref_call(varargin)
        % Callback for choosing reference trace
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        % load new profile
        prnum=str2num(S.refprof.String);
        test=S.profiles(min(abs(prnum-S.profiles))==abs(prnum-S.profiles));
        prnum=test(1);
        S.refprof.String=int2str(prnum);
        temp=readmala(S.foldername,S.name,prnum);
        chnum=str2num(S.refch.String);
        if chnum>S.chan
            chnum=S.chan;
            S.refch.String=int2str(S.chan);
        elseif chnum<1
            chnum=1;
            S.refchn.String='1';
        end
        trnum=str2num(S.reftrnum.String);
        if trnum>length(temp{1}(1,:))
            trnum=length(temp{1}(1,:));
            S.reftrnum.String=num2str(trnum);
        elseif trnum<1
            trnum=1;
            S.reftrnum.String='1';
        end
        % Process data
        % Amplitude offset removal:
        tempdata=DCremoval(temp{chnum}(:,trnum),S.t);
        % Bandpassfilter:
        tempdata=bandpass_gpr(tempdata,S.dt,100,600);
        % update reference trace plot:
        set(S.Plot3a,'Xdata',tempdata);
        % update handles
        guidata(gcf,S);
        % call new plot
        plotnew_call(varargin);
    end

    function [] = save_call(varargin)
        % Callback for saving of settings
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        % Settings speichern
        S.flag=S.pop.Value;
        fid=fopen(fullfile(S.foldername,[S.name,'_t0check_Settings.txt']),'wt');
        fprintf(fid,'Settings determined with Check_t0.m\n');
        if S.flag==1 % channelshift for each profile
            fprintf(fid,'Channelshift for each profile individually\n');
        elseif S.flag==2 % chhanelshift for all profiles
            fprintf(fid,'The same channelshift for all profiles\n');
            fprintf(fid,['Shiftsamples = [',int2str(S.shiftsamples'),'] for reference profile ',S.refprof.String,'\n']);
        elseif S.flag==3 % t0pick and cross-correlation
            fprintf(fid,'t0 pick and cross-correlation with reference trace\n');
            fprintf(fid,['Reference trace: Profile #',S.refprof.String,', Channel #',S.refch.String,', Trace #',S.reftrnum.String,'\n']);
            fprintf(fid,['t0 = ',S.t0field.String,' ns\n']);
        elseif S.flag==4 % threshold
            fprintf(fid,'Use amplitude threshold\n');
            fprintf(fid,['Threshold = ',S.thr.String,'\n']);
        end
        if S.norm==1
            fprintf(fid,'Use normalized data for t0 correction');
        else
            fprintf(fid,'Use non-normalized data for t0 correction');
        end
        fclose(fid);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'arrow');
        % update handles
        guidata(gcf,S);
    end

    function [] = plotnew_call(varargin)
        % Callback for applying t0 correction and new plotting
        S=guidata(gcbf);
        % PLease wait message
        S.info.String=['PLEASE WAIT...'];
        % get method:
        S.method=S.pop.Value;
        
        % t0 correction
        S.alltraces_t0=zeros(S.ns,S.traceanz*S.chan);
        if S.norm==1
            for ii=1:S.chan
                normtraces{ii}=normalize2d(S.traces_proc{ii},0.98);
            end
        end
        if S.method==1 % channelshift for each profile individually
        	if S.norm==1
                [S.traces_t0,S.shiftsamples]=channelshift(normtraces,40);
                for ii=1:S.chan
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            else
                [S.traces_t0,S.shiftsamples]=channelshift(S.traces_proc,40);
                for ii=1:S.chan
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            end   
            S.info.String=['Shiftsamples = [',int2str(S.shiftsamples'),'] for current profile'];
        elseif S.method==2  % channelshift for all profiles the same
            [refprof]=readmala(S.foldername,S.name,str2num(S.refprof.String));
            for ii=1:S.chan
                % Amplitude offset removal:
                refprof{ii}=DCremoval(refprof{ii},S.t);
                % Bandpassfilter:
                refprof{ii}=bandpass_gpr(refprof{ii},S.dt,100,600);
            end
            if S.norm==1
                [temp,S.shiftsamples]=channelshift(normtraces,40); % determine shiftsamples
                [S.traces_t0]=channelshift(normtraces,40,S.shiftsamples);
                for ii=1:S.chan
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            else
                [temp,S.shiftsamples]=channelshift(refprof,40); % determine shiftsamples
                [S.traces_t0]=channelshift(S.traces_proc,40,S.shiftsamples);
                for ii=1:S.chan
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end      
            end
            S.info.String=['Reference profile #',S.refprof.String,': Shiftsamples=[',int2str(S.shiftsamples'),']'];
        elseif S.method==3  % pick t0 on reference trace and crosscorrelation
            S.t0pick=str2num(S.t0field.String);
            if S.norm==1
                for ii=1:S.chan
                    S.traces_t0{ii}=t0correction(normtraces{ii},normalize2d(S.reftrace,0.98),S.t0pick,S.dt);
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            else
                for ii=1:S.chan
                    S.traces_t0{ii}=t0correction(S.traces_proc{ii},S.reftrace,S.t0pick,S.dt);
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            end
            S.info.String=['Picked t0 = ',num2str(S.t0pick,2),' ns'];
        elseif S.method==4  % using amplitude threshold
            S.threshhold=str2num(S.thr.String);
            if S.norm==1
                for ii=1:S.chan
                    [S.traces_t0{ii}]=t0corr_thresh(normtraces{ii},S.threshhold);
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            else
                for ii=1:S.chan
                    [S.traces_t0{ii}]=t0corr_thresh(S.traces_proc{ii},S.threshhold);
                    S.alltraces_t0(:,1+(ii-1)*S.traceanz:ii*S.traceanz) = S.traces_t0{ii};
                end
            end
            S.info.String=['Amplitude threshold = ',num2str(S.threshhold,1)];
        end
        
        % update plots
        set(S.Plot2a,'CData',S.alltraces_t0,'XData',1:S.traceanz*S.chan);
        set(S.Plot1a,'CData',S.alltraces_proc,'XData',1:S.traceanz*S.chan);
        set(S.ax1,'YLim',[0 str2num(S.tmax.String)],'XLim',[1 length(S.alltraces_proc(1,:))],'XTick',[1:length(S.alltraces_proc(1,:))/length(S.traces):length(S.alltraces_proc(1,:))],'XTickLabel',1:length(S.traces))
        set(S.ax2,'YLim',[0 str2num(S.tmax.String)],'XLim',[1 length(S.alltraces_proc(1,:))],'XTick',[1:length(S.alltraces_t0(1,:))/length(S.traces):length(S.alltraces_t0(1,:))],'XTickLabel',1:length(S.traces))
        for ii=1:S.chan
            set(S.axw1.Children(ii),'XData',S.traces_proc{ii}(:,1))
            set(S.axw2.Children(ii),'XData',S.traces_t0{ii}(:,1))
        end
        drawnow;
        
        % chnage pointer back to arrow
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'arrow')
        
        % update handles
        guidata(gcf,S);
    end


    function [] = thr_call(varargin)
        % Callback for edit field for threshhold
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        S.threshhold=str2num(S.thr.String);
        % determine t0 on reference trace
        k=1;
        if S.threshhold<0
            while k<=length(S.reftrace) && S.reftrace(k)>S.threshhold
                k=k+1;
            end
        else
            while k<=length(S.reftrace) && S.reftrace(k)<S.threshhold
                k=k+1;
            end
        end
        S.t0=k*S.dt;
        % update reference trace plot
        if isfield(S,'Plot3b')
            set(S.Plot3b,'Ydata',[S.t0 S.t0]);
        else
            S.Plot3b=plot(S.ax3,[min(S.reftrace) max(S.reftrace)],[S.t0 S.t0],'r','Linewidth',2);
        end
        drawnow;
        % update handles
        guidata(gcf,S);
        plotnew_call(varargin); % call plotting function
    end


    function [] = tmax_call(varargin)
        % Callback for edit field tmax
        S=guidata(gcbf);
        tmax=str2num(S.tmax.String);
        set(S.ax2,'YLim',[0 tmax]);
        set(S.ax1,'YLim',[0 tmax]);
        set(S.axw1,'YLim',[0 tmax]);
        set(S.axw2,'YLim',[0 tmax]);
        drawnow;
    end


    function [] = pop_call(varargin)
        % Callback for t0-method popupmenu
        S=guidata(gcbf);
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch');
        S.flag=S.pop.Value;
        S.save.Enable='on';
        % Update t0 on reference trace plot (top)
        if S.flag==1 % channelshift individually
            set(S.thr,'Enable','off');
            set(S.reftrnum,'Enable','off');
            set(S.refch,'Enable','off');
            set(S.refprof,'Enable','off');
            set(S.t0field,'Enable','off');
            if length(S.ax3.Children)>1
                delete(S.ax3.Children(1))
            end
        elseif S.flag==2 %channelshift for all the same
            set(S.thr,'Enable','off');
            set(S.reftrnum,'Enable','off');
            set(S.refch,'Enable','off');
            set(S.refprof,'Enable','on');
            set(S.t0field,'Enable','off');
            if length(S.ax3.Children)>1
                delete(S.ax3.Children(1))
            end
        elseif S.flag==3    % t0 and reference trace
            set(S.thr,'Enable','off');
            set(S.reftrnum,'Enable','on');
            set(S.refch,'Enable','on');
            set(S.refprof,'Enable','on');
            set(S.t0field,'Enable','on');
            S.t0pick=str2num(S.t0field.String);
            if length(S.ax3.Children)>1
                set(S.Plot3b,'Ydata',[S.t0pick S.t0pick]);
            else
                hold(S.ax3,'on');
                S.Plot3b=plot(S.ax3,[min(S.reftrace) max(S.reftrace)],[S.t0pick S.t0pick],'r','Linewidth',2);
            end
            drawnow;
        elseif S.flag==4 % amplitude threshold
            set(S.thr,'Enable','on');
            set(S.reftrnum,'Enable','off');
            set(S.refch,'Enable','off');
            set(S.refprof,'Enable','off');
            set(S.t0field,'Enable','off');
            S.threshhold=str2num(S.thr.String);
            % determine t0 on reference trace
            k=1;
            if S.threshhold<0
                while k<=length(S.reftrace) && S.reftrace(k)>S.threshhold
                    k=k+1;
                end
            else
                while k<=length(S.reftrace) && S.reftrace(k)<S.threshhold
                    k=k+1;
                end
            end
            S.t0=k*S.dt;
            % update reference trace plot
            if length(S.ax3.Children)>1
                set(S.Plot3b,'Ydata',[S.t0 S.t0]);
            else
                hold(S.ax3,'on');
                S.Plot3b=plot(S.ax3,[min(S.reftrace) max(S.reftrace)],[S.t0 S.t0],'r','Linewidth',2);
            end
            drawnow;
        end
        % update handles
        guidata(gcf,S);
        plotnew_call(varargin); % call plotting function
    end


    function [] = profnum_call(varargin)
        % Callback for profile number listbox
        S=guidata(gcbf);
        S.profnum=S.profiles(S.prof.Value); % current profile number
        % change pointer
        set(findobj('Type','Figure','Name','Traces'), 'pointer', 'watch')
        drawnow;
        
        % read new data
        [S.traces]=readmala(S.foldername,S.name,S.profnum);
        S.traceanz=length(S.traces{1}(1,:));
        % Process data
        S.alltraces_proc=zeros(S.ns,S.traceanz*length(S.traces));
        for i=1:length(S.traces)
            % Amplitude offset removal:
            S.traces_proc{i}=DCremoval(S.traces{i},S.t);
            % Bandpassfilter:
            S.traces_proc{i}=bandpass_gpr(S.traces_proc{i},S.dt,100,600);
            
            S.alltraces_proc(:,1+(i-1)*S.traceanz:i*S.traceanz) = S.traces_proc{i};
        end
        
        % update handles
        guidata(gcf,S);
        plotnew_call(varargin); % call plotting function
    end
end