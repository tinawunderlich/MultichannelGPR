classdef CMPproc_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        GridLayout               matlab.ui.container.GridLayout
        LeftPanel                matlab.ui.container.Panel
        LoadCMPdataButton        matlab.ui.control.Button
        HyperbolafittingLabel    matlab.ui.control.Label
        SavevelocitiesButton     matlab.ui.control.Button
        t0nsSlider               matlab.ui.control.Slider
        t0nsSliderLabel          matlab.ui.control.Label
        vcmnsSlider              matlab.ui.control.Slider
        vcmnsSliderLabel         matlab.ui.control.Label
        PickApex                 matlab.ui.control.Button
        CorrectT0_Button         matlab.ui.control.Button
        dxoffset_EditField       matlab.ui.control.EditField
        dxoffsetmonesidedLabel   matlab.ui.control.Label
        Startoffset_EditField    matlab.ui.control.EditField
        StartoffsetTxRxmLabel    matlab.ui.control.Label
        ProfileListBox           matlab.ui.control.ListBox
        ProfileListBoxLabel      matlab.ui.control.Label
        PlotButtonBox            matlab.ui.container.ButtonGroup
        AspectRatioListBox       matlab.ui.control.ListBox
        AspectRatioListBoxLabel  matlab.ui.control.Label
        ScaleEditField           matlab.ui.control.EditField
        ScaleEditFieldLabel      matlab.ui.control.Label
        GrayscaleButton          matlab.ui.control.RadioButton
        WiggleButton             matlab.ui.control.RadioButton
        RightPanel               matlab.ui.container.Panel
        HowTo_Label              matlab.ui.control.Label
        SemblanceanalysisPanel   matlab.ui.container.Panel
        calculateSemblanceButton matlab.ui.control.Button
        SavevelocitiesButton_2   matlab.ui.control.Button
        ManualPick_Button        matlab.ui.control.Button
        ThresholdEditField       matlab.ui.control.EditField
        ThresholdEditFieldLabel  matlab.ui.control.Label
        epsEditField       matlab.ui.control.EditField
        epsEditFieldLabel  matlab.ui.control.Label
        minptsEditField       matlab.ui.control.EditField
        minptsEditFieldLabel  matlab.ui.control.Label
        AutoPick_Button          matlab.ui.control.Button
        vMax_EditField           matlab.ui.control.EditField
        v_maxcmnsEditFieldLabel  matlab.ui.control.Label
        vMin_EditField           matlab.ui.control.EditField
        v_mincmnsEditFieldLabel  matlab.ui.control.Label
        UISemblance              matlab.ui.control.UIAxes
        UIData                   matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        data
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadCMPdataButton
        function LoadData_cb(app, event)

            % get folder name with data
            if ispc
                if exist('temp.temp','file') % read last opened folder from temp.temp
                    fid=fopen('temp.temp','r');
                    fn=textscan(fid,'%s');
                    fclose(fid);
                    if ~isempty(fn{1})
                        pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
                    else
                        pfad=uigetdir([],'Choose folder with radargrams');
                    end
                    fileattrib('temp.temp','-h');
                    fid=fopen('temp.temp','wt');
                    fprintf(fid,'%s',pfad);
                    fclose(fid);
                    fileattrib('temp.temp','+h');
                else
                    pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

                    fid=fopen('temp.temp','wt');
                    fprintf(fid,'%s',pfad);
                    fclose(fid);
                    fileattrib('temp.temp','+h');
                end
            else
                if exist('.temp.temp','file') % read last opened folder from temp.temp
                    fid=fopen('.temp.temp','r');
                    fn=textscan(fid,'%s');
                    fclose(fid);
                    if ~isempty(fn{1})
                        pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
                    else
                        pfad=uigetdir([],'Choose folder with radargrams');
                    end
                else
                    pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
                end

                fid=fopen('.temp.temp','wt');
                fprintf(fid,'%s',pfad);
                fclose(fid);
            end
            focus(app.UIFigure);
            app.data.path=pfad; % path to data

            % Load data
            load(fullfile(pfad,'radargrams.mat'));
            load(fullfile(pfad,'x.mat'));
            load(fullfile(pfad,'global_coords.mat'));
            load(fullfile(pfad,'t.mat'));
            app.data.radargrams=radargrams;
            app.data.x=x;
            app.data.global_coords=global_coords;
            app.data.t=t;

            % create profile list
            app.ProfileListBox.Items=cellfun(@(x) num2str(x),num2cell(1:length(radargrams)),'UniformOutput',false);
            app.ProfileListBox.Value='1';

            % make plot of data
            plot_new_data(app);

            % enable all other buttons
            app.AspectRatioListBox.Enable='on';
            app.AspectRatioListBoxLabel.Enable='on';
            app.CorrectT0_Button.Enable='on';
            app.calculateSemblanceButton.Enable='on';
            app.GrayscaleButton.Enable='on';
            app.HyperbolafittingLabel.Enable='on';
            app.PickApex.Enable='on';
            app.PlotButtonBox.Enable='on';
            app.ProfileListBox.Enable='on';
            app.ProfileListBoxLabel.Enable='on';
            app.SemblanceanalysisPanel.Enable='on';
            app.ScaleEditField.Enable='on';
            app.ScaleEditFieldLabel.Enable='on';
            app.StartoffsetTxRxmLabel.Enable='on';
            app.Startoffset_EditField.Enable='on';
            app.t0nsSlider.Enable='on';
            app.t0nsSliderLabel.Enable='on';
            app.WiggleButton.Enable='on';
            app.dxoffset_EditField.Enable='on';
            app.dxoffsetmonesidedLabel.Enable='on';
            app.vcmnsSlider.Enable='on';
            app.vcmnsSliderLabel.Enable='on';
            app.vMax_EditField.Enable="on";
            app.vMin_EditField.Enable="on";
            app.v_maxcmnsEditFieldLabel.Enable="on";
            app.v_mincmnsEditFieldLabel.Enable="on";

            % make new lint
            app.HowTo_Label.Text = 'Check offsets for your data';

        end

        % plot new data
        function plot_new_data(app, event)
            app.HowTo_Label.Text='';
            % if new profile or different offsets -> delete semblance plot
            if exist('event','var') && (strcmp(event.Source.Tag,'Profilelist') || strcmp(event.Source.Tag,'Startoffset') || strcmp(event.Source.Tag,'dxoffset'))
                c=app.UISemblance.Children;
                for i=1:length(c)
                    delete(c(i));
                end
            end
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get t0:
            if isfield(app.data,'t0') && length(app.data.t0)>=num
                t0=app.data.t0(num);
                % make new t axis
                app.data.tnew{num}=app.data.t(app.data.t>=t0);
                app.data.tnew{num}=app.data.tnew{num}-min(app.data.tnew{num}); % start again with 0
                % correct radargram
                app.data.radargrams_corr{num}=app.data.radargrams{num}(app.data.t>=t0,:);
                flag=1;
            else
                t0=0; % no correction yet
                app.data.tnew{num}=app.data.t;
                flag=0;
            end
            % get x axis values:
            start=str2double(app.Startoffset_EditField.Value)/2;
            dx=str2double(app.dxoffset_EditField.Value);
            ende=str2double(app.Startoffset_EditField.Value)/2+str2double(app.dxoffset_EditField.Value)*(size(app.data.radargrams{num},2)-1);
            app.data.xprof{num}=start:dx:ende;
            % scale:
            scale=str2double(app.ScaleEditField.Value);
            % delete old plots:
            h=app.UIData.Children; % handles to children
            for i=1:length(h)
                delete(h(i));
            end
            % plot data:
            if app.WiggleButton.Value==true
                if flag==1 % t0-corrected data
                    wigglesc(app.data.radargrams_corr{num},app.data.tnew{num},app.data.xprof{num},scale,app.UIData);
                else
                    wigglesc(app.data.radargrams{num},app.data.t,app.data.xprof{num},scale,app.UIData);
                end
            else
                if flag==1 % t0-corrected data
                    imagesc(app.UIData,app.data.xprof{num},app.data.tnew{num},app.data.radargrams_corr{num});
                    coldata=sort(unique(app.data.radargrams_corr{num}));
                else
                    imagesc(app.UIData,app.data.xprof{num},app.data.t,app.data.radargrams{num});
                    coldata=sort(unique(app.data.radargrams{num}));
                end
                % color settings:
                colormap(app.UIData,flipud(gray))
                coldata(isnan(coldata))=[]; % delete nans
                cmin=coldata(round(length(coldata)/100*scale));
                cmax=coldata(end-round(length(coldata)/100*scale));
                app.UIData.CLim=[cmin cmax];
            end
            % axes settings:
            app.UIData.YDir='reverse';
            % get aspect ratio:
            temp=app.AspectRatioListBox.Value;
            temp2=strsplit(temp,'/');
            app.UIData.DataAspectRatio=[str2double(temp2{1}) str2double(temp2{2}) 1];
            if flag==1 % t0 corrected data:
                app.UIData.YLim=[0 max(app.data.tnew{num})];
            else
                app.UIData.YLim=[0 max(app.data.t)];
            end
            app.UIData.XLim=[0 max(app.data.xprof{num})];
        end

        % Callback for calculating semblance
        function semblance_cb(app, event)
            app.HowTo_Label.Text='Calculating semblance... this may take a while!';
            pause(0.01);
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % semblance
            vmin=str2double(app.vMin_EditField.Value);
            vmax=str2double(app.vMax_EditField.Value);
            vs=vmin:0.2:vmax;
            if isfield(app.data,'t0') && length(app.data.t0)>=num
                t0s=app.data.tnew{num};
                dat=app.data.radargrams_corr{num}; % corrected t0 data
            else
                t0s=app.data.t;
                dat=app.data.radargrams{num};
            end
            dt=t0s(2)-t0s(1);
            for i=1:length(vs) % loop over velocities
                % nmo correction:
                nmogather=zeros(size(dat));
                nmo=zeros(size(dat,2),size(dat,2));
                % for all t0 calculate hyperbola
                thyp=2*sqrt((t0s'/2).^2+4*app.data.xprof{num}.^2./(vs(i)/100)^2); % t for each xprof
                datin=NaN(4,size(thyp,2));
                t0sind=NaN(size(datin));
                for j=1:length(t0s)
                    for k=1:size(thyp,2) % for all xprof
                        ind=find(t0s>=thyp(j,k)-2*dt & t0s<=thyp(j,k)+2*dt);  % indices of neighboring amplitude values (in most cases 4 values)
                        if ~isempty(ind)
                            ind=ind(1:min(size(ind,2),4)); % make sure that no longer than 4 elements
                            temp=t0s(ind)';
                            t0sind(1:size(temp,1),k)=temp; % times of these samples
                            datin(1:size(temp,1),k)=dat(ind,k); % amplitude values for these samples
                        end
                    end
                    % reduce samples to t=0ns for all xprof -> only one input for interp1
                    thyp(j,:)=thyp(j,:)-min(t0sind,[],1);
                    t0sind=t0sind-min(t0sind,[],1); % also reduce query points
                    % make nmo correction
                    [a,b]=unique(t0sind(:,1)); % only get unique sample points for interpolation
                    if any(~isfinite(a))
                        b=b(isfinite(a));
                        a=a(isfinite(a));
                    end
                    nmo=interp1(a,datin(b,:),thyp(j,:)); % interpolate for each xprof between the 4 samples
                    nmogather(j,:)=diag(nmo);
                end
                % semblance calculation:
                M=size(dat,2); % number of traces
                N=7; % number of samples for smoothing
                innersumtop=sum(nmogather,2,'omitnan').^2;
                innersumbottom=sum(nmogather.^2,2,'omitnan');
                for j=1:length(t0s)
                    ind=j-N:j+N;
                    ind=ind(ind>0 & ind<=size(dat,1)); % only positive indices and smaller than tmax
                    S(j,i)=sum(innersumtop(ind))/(M*sum(innersumbottom(ind)));
                end
            end
            % interpolate semblance, make sampling of v smaller
            [xg,yg]=meshgrid(vs,t0s);
            app.data.S=interp2(xg,yg,S,repmat(vmin:0.05:vmax,[length(t0s) 1]),repmat(t0s',[1 length(vmin:0.05:vmax)]));
            % plot semblance
            imagesc(app.UISemblance,vs,t0s,app.data.S);
            set(app.UISemblance,'XLim',[vmin vmax],'YLim',[0 max(t0s)])
            colorbar(app.UISemblance);
            % set same height as other plot
            app.UISemblance.DataAspectRatio=[1 3.5 1];

            % enable other buttons:
            app.AutoPick_Button.Enable='on';
            app.ManualPick_Button.Enable='on';
            app.ThresholdEditField.Enable='on';
            app.ThresholdEditFieldLabel.Enable='on';
            app.vMax_EditField.Enable='on';
            app.vMin_EditField.Enable='on';
            app.v_maxcmnsEditFieldLabel.Enable="on";
            app.v_mincmnsEditFieldLabel.Enable="on";
            app.epsEditField.Enable="on";
            app.epsEditFieldLabel.Enable="on";
            app.minptsEditField.Enable='on';
            app.minptsEditFieldLabel.Enable="on";
            app.HowTo_Label.Text='';
            pause(0.01);
        end

        % Callback function for correcting t0
        function correctT0_cb(app, event)
            app.HowTo_Label.Text = 'Pick two points along first break';
            % Change mouse pointer
            set(app.UIFigure,'Pointer','cross');
            drawnow;
            % make picks:
            temp=drawpoint(app.UIData);
            picks=[temp.Position(1:2)];
            temp=drawpoint(app.UIData);
            picks=[picks; temp.Position(1:2)];
            % change mouse pointer back:
            set(app.UIFigure,'Pointer','arrow');
            drawnow;
            % fit linear function:
            p=polyfit(picks(:,1),picks(:,2),1); 
            % check for correctness of velocity (offsets wrong?):
            vdirect=1/p(1); % m/ns
            if vdirect<0.03 || vdirect>0.2 % strange velocity
                uialert(app.UIFigure,['Velocity for ground wave is ',num2str(app.data.t0(num),2),'m/ns. Please check if the offsets are correct and start again!'],'Velocity error!','Icon','warning');
            else
                % Continue with t0 shift
                % get radargram number:
                num=str2double(app.ProfileListBox.Value);
                % set t0 for this profile
                app.data.t0(num)=p(2); % t0 for this profile
                % message for user:
                uialert(app.UIFigure,['Data of this profile is shifted for t0 = ',num2str(app.data.t0(num),2),'ns.'],'t0 corrected!','Icon','success');
                plot_new_data(app);
            end
        end

        % Callback for picking apex
        function pickApex_cb(app, event)
            app.HowTo_Label.Text='Pick apex point in data plot';
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % pick apex point
            temp=drawpoint(app.UIData);
            t0=temp.Position(2);
            delete(app.UIData.Children(1)); % delete point
            % get current v
            v=app.vcmnsSlider.Value; % cm/ns
            % set t0 slider:
            app.t0nsSlider.Value=t0;
            % calculate hyperbola
            thyp=2*sqrt((t0/2)^2+4*[0 app.data.xprof{num}].^2./(v/100)^2);
            % plot hyperbola
            app.UIData.NextPlot='add';
            app.data.hyp=plot(app.UIData,[0 app.data.xprof{num}],thyp,'r','Linewidth',2);
            % enable save button:
            app.SavevelocitiesButton.Enable='on';
            % make new text:
            app.HowTo_Label.Text=['Current hyperbola: t0 = ',num2str(t0,'%.2f'),'ns, v = ',num2str(v,'%.1f'),'cm/ns'];
        end

        % Callback for v-slider
        function vslider_cb(app, event)
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get current v
            v=app.vcmnsSlider.Value; % cm/ns
            % get t0
            t0=app.t0nsSlider.Value;
            % calculate hyperbola
            thyp=2*sqrt((t0/2)^2+4*[0 app.data.xprof{num}].^2./(v/100)^2);
            app.data.hyp.YData=thyp;
            % update text:
            app.HowTo_Label.Text=['Current hyperbola: t0 = ',num2str(t0,'%.2f'),'ns, v = ',num2str(v,'%.1f'),'cm/ns'];
        end

        % Callback for t0-slider
        function t0slider_cb(app, event)
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get current v
            v=app.vcmnsSlider.Value; % cm/ns
            % get t0
            t0=app.t0nsSlider.Value;
            % calculate hyperbola
            thyp=2*sqrt((t0/2)^2+4*[0 app.data.xprof{num}].^2./(v/100)^2);
            app.data.hyp.YData=thyp;
            % update text:
            app.HowTo_Label.Text=['Current hyperbola: t0 = ',num2str(t0,'%.2f'),'ns, v = ',num2str(v,'%.1f'),'cm/ns'];
        end
        
        % Callback for moving sliders with arrow keys
        function arrowkeys(app, event)
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get current v
            v=app.vcmnsSlider.Value; % cm/ns
            % get t0
            t0=app.t0nsSlider.Value;
            if strcmp(event.Key,'leftarrow')
                event.Source.CurrentObject.Value=event.Source.CurrentObject.Value-0.1;
            elseif strcmp(event.Key,'rightarrow')
                event.Source.CurrentObject.Value=event.Source.CurrentObject.Value+0.1;
            end
            if strcmp(event.Source.CurrentObject.Tag,'t0-slider')
                t0slider_cb(app);
            elseif strcmp(event.Source.CurrentObject.Tag,'v-slider')
                vslider_cb(app);
            end
        end

        % Callback for saving velocity from hyperbola
        function save_v_cb(app, event)
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get current v
            v=app.vcmnsSlider.Value; % cm/ns
            % get t0
            t0=app.t0nsSlider.Value;
            % file identifier
            if ~exist(fullfile(app.data.path,'CMP_data.txt'),'file')
                fid=fopen(fullfile(app.data.path,'CMP_data.txt'),'wt');
                fprintf(fid,'CMP data analysis\n');
                fprintf(fid,'Methods: 1=Manual hyperbola fitting, 2=Semblance manual picking, 3=Semblance auto picking\n')
                fprintf(fid,'Profile#\tt0[ns]\tv[cm/ns]\tMethod\n');
            else
                fid=fopen(fullfile(app.data.path,'CMP_data.txt'),'a');
            end
            fprintf(fid,'%d\t%.2f\t%.2f\t%d\n',[num t0 v 1]);
            fclose(fid);
        end

        % Callback for saving velocity from semblance
        function save_v_S_cb(app, event)
            % get radargram number:
            num=str2double(app.ProfileListBox.Value);
            % get v and t picks
            if isfield(app.data,'S_auto_v') && ~isempty(app.data.S_auto_v)
                % auto picks
                method=3;
                v=app.data.S_auto_v;
                t0=app.data.S_auto_t;
            elseif isfield(app.data,'S_manu_v') && ~isempty(app.data.S_manu_v)
                % manual picks
                method=2;
                v=app.data.S_manu_v;
                t0=app.data.S_manu_t;
            end
            % file identifier
            if ~exist(fullfile(app.data.path,'CMP_data.txt'),'file')
                fid=fopen(fullfile(app.data.path,'CMP_data.txt'),'wt');
                fprintf(fid,'CMP data analysis\n');
                fprintf(fid,'Methods: 1=Manual hyperbola fitting, 2=Semblance manual picking, 3=Semblance auto picking\n')
                fprintf(fid,'Profile#\tt0[ns]\tv[cm/ns]\tMethod\n');
            else
                fid=fopen(fullfile(app.data.path,'CMP_data.txt'),'a');
            end
            for i=1:length(v)
                fprintf(fid,'%d\t%.2f\t%.2f\t%d\n',[num t0(i) v(i) method]);
            end
            fclose(fid);
        end

        % Callback for Auto pick semblance
        function autopick_cb(app, event)
            app.HowTo_Label.Text='Semblance maxima are automatically picked. Adjust threshold, eps and minpts according to your data.';
            % delete possible picks:
            c=app.UISemblance.Children;
            for i=1:length(c)-1
                delete(c(i));
            end
            app.data.S_auto_v=[];
            app.data.S_auto_t=[];
            app.data.S_manu_v=[];
            app.data.S_manu_t=[];
            % get semblance:
            S=app.data.S;
            % get velocity range:
            vmin=str2double(app.vMin_EditField.Value);
            vmax=str2double(app.vMax_EditField.Value);
            vs=linspace(vmin,vmax,size(S,2));
            % get threshold:
            thresh=str2double(app.ThresholdEditField.Value);
            % get semblance above this value:
            if max(S,[],'all')<thresh
                uialert(app.UIFigure,'No data above threshold. Decrease threshold!','Error');
            else
                S(S<thresh)=0;
                S(isnan(S))=0;
                % get maxima: (for each time step)
                mv=NaN(length(app.data.t),1);
                tv=NaN(length(app.data.t),1);
                for i=1:size(S,1) % go through rows of S (=times)
                    if ~all(S(i,:)==0)
                        ind=islocalmax(S(i,:),'MaxNumExtrema',1);
                        mv(i)=vs(ind); % v at maximum
                        tv(i)=app.data.t(i); % corresponding time
                    end
                end
                % combine maxima:
                eps=str2double(app.epsEditField.Value);
                minpts=str2double(app.minptsEditField.Value);
                idx=dbscan([mv,tv],eps,minpts); % apply dbscan to find cluster
                for i=1:max(idx) % for each cluster
                    % find points in this cluster
                    vcluster=mv(idx==i);
                    tcluster=tv(idx==i);
                    % get max(S) in this area
                    Sbox=S(app.data.t>=min(tcluster) & app.data.t<=max(tcluster),vs>=min(vcluster) & vs<=max(vcluster)); % semblance in this range
                    vtemp=vs(vs>=min(vcluster) & vs<=max(vcluster)); % v for this box
                    ttemp=app.data.t(app.data.t>=min(tcluster) & app.data.t<=max(tcluster)); % t for this box
                    [r,c]=find(max(Sbox,[],'all')==Sbox); % inidces of maximum in this box
                    vmaxC(i)=vtemp(c);
                    tmaxC(i)=ttemp(r);
                end
                app.data.S_auto_v=vmaxC;
                app.data.S_auto_t=tmaxC;
                % plot maximum points:
                app.UISemblance.NextPlot='add';
                plot(app.UISemblance,vmaxC,tmaxC,'k*','Linewidth',2)
            end
            % enable saving
            app.SavevelocitiesButton_2.Enable='on';
        end

        % Callback for manually pick semblance
        function manualpick_cb(app, event)
            app.HowTo_Label.Text='Pick maxima in the semblance plot';
            % delete possible picks:
            c=app.UISemblance.Children;
            for i=1:length(c)-1
                delete(c(i));
            end
            app.data.S_auto_v=[];
            app.data.S_auto_t=[];
            app.data.S_manu_v=[];
            app.data.S_manu_t=[];
            % Change mouse pointer
            set(app.UIFigure,'Pointer','cross');
            drawnow;
            % make picks:
            yes=1;
            while yes==1
                temp=drawpoint(app.UISemblance);
                temp.Color='k';
                app.data.S_manu_v=[app.data.S_manu_v; temp.Position(1)];
                app.data.S_manu_t=[app.data.S_manu_t; temp.Position(2)];
                selection = uiconfirm(app.UIFigure,'Do you want to pick another maximum?','Continue?', ...
                    "Options",["Yes","No"], ...
                    "DefaultOption",1);
                yes=strcmp(selection,'Yes');
            end
            % change mouse pointer back:
            set(app.UIFigure,'Pointer','arrow');
            drawnow;
            % plot maximum points:
            app.UISemblance.NextPlot='add';
            plot(app.UISemblance,app.data.S_manu_v,app.data.S_manu_t,'k*','Linewidth',2)
            % enable saving
            app.SavevelocitiesButton_2.Enable='on';
            app.HowTo_Label.Text='';
        end


        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {597, 597};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {220, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 661 597];
            app.UIFigure.Name = 'CMP Processing';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
            app.UIFigure.KeyPressFcn = createCallbackFcn(app, @arrowkeys, true);
            app.UIFigure.WindowKeyPressFcn = createCallbackFcn(app, @arrowkeys, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {220, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create PlotButtonBox (ButtonGroup for wiggle/grayscale!)
            app.PlotButtonBox = uibuttongroup(app.LeftPanel);
            app.PlotButtonBox.Enable = 'off';
            app.PlotButtonBox.Title = 'Plot options';
            app.PlotButtonBox.Position = [14 311 192 114];
            app.PlotButtonBox.SelectionChangedFcn = createCallbackFcn(app, @plot_new_data, true);

            % Create WiggleButton
            app.WiggleButton = uiradiobutton(app.PlotButtonBox);
            app.WiggleButton.Enable = 'off';
            app.WiggleButton.Text = 'Wiggle';
            app.WiggleButton.Position = [11 64 59 22];
            app.WiggleButton.Value = true;

            % Create GrayscaleButton
            app.GrayscaleButton = uiradiobutton(app.PlotButtonBox);
            app.GrayscaleButton.Enable = 'off';
            app.GrayscaleButton.Text = 'Grayscale';
            app.GrayscaleButton.Position = [11 36 76 22];

            % Create ScaleEditFieldLabel
            app.ScaleEditFieldLabel = uilabel(app.PlotButtonBox);
            app.ScaleEditFieldLabel.HorizontalAlignment = 'right';
            app.ScaleEditFieldLabel.Enable = 'off';
            app.ScaleEditFieldLabel.Position = [103 51 35 22];
            app.ScaleEditFieldLabel.Text = 'Scale';

            % Create ScaleEditField
            app.ScaleEditField = uieditfield(app.PlotButtonBox, 'text');
            app.ScaleEditField.Enable = 'off';
            app.ScaleEditField.Position = [148 52 39 21];
            app.ScaleEditField.Value='1';
            app.ScaleEditField.ValueChangedFcn = createCallbackFcn(app, @plot_new_data, true);

            % Create AspectRatioListBoxLabel
            app.AspectRatioListBoxLabel = uilabel(app.PlotButtonBox);
            app.AspectRatioListBoxLabel.HorizontalAlignment = 'right';
            app.AspectRatioListBoxLabel.Enable = 'off';
            app.AspectRatioListBoxLabel.Position = [5 6 74 22];
            app.AspectRatioListBoxLabel.Text = 'Aspect Ratio';

            % Create AspectRatioListBox
            app.AspectRatioListBox = uilistbox(app.PlotButtonBox);
            app.AspectRatioListBox.Items = {'1/1', '1/2', '1/5', '1/10', '1/20', '1/30', '1/40', '1/50', '1/60'};
            app.AspectRatioListBox.Enable = 'off';
            app.AspectRatioListBox.Position = [94 6 93 24];
            app.AspectRatioListBox.Value = '1/30';
            app.AspectRatioListBox.ValueChangedFcn = createCallbackFcn(app, @plot_new_data, true);

            % Create ProfileListBoxLabel
            app.ProfileListBoxLabel = uilabel(app.LeftPanel);
            app.ProfileListBoxLabel.HorizontalAlignment = 'right';
            app.ProfileListBoxLabel.Enable = 'off';
            app.ProfileListBoxLabel.Position = [16 512 39 22];
            app.ProfileListBoxLabel.Text = 'Profile';

            % Create ProfileListBox
            app.ProfileListBox = uilistbox(app.LeftPanel);
            app.ProfileListBox.Items = {};
            app.ProfileListBox.Enable = 'off';
            app.ProfileListBox.Position = [70 497 136 39];
            app.ProfileListBox.Value = {};
            app.ProfileListBox.ValueChangedFcn = createCallbackFcn(app, @plot_new_data, true);
            app.ProfileListBox.Tag='Profilelist';

            % Create StartoffsetTxRxmLabel
            app.StartoffsetTxRxmLabel = uilabel(app.LeftPanel);
            app.StartoffsetTxRxmLabel.Enable = 'off';
            app.StartoffsetTxRxmLabel.Position = [18 466 118 22];
            app.StartoffsetTxRxmLabel.Text = 'Start offset Tx-Rx [m]';

            % Create Startoffset_EditField
            app.Startoffset_EditField = uieditfield(app.LeftPanel, 'text');
            app.Startoffset_EditField.Enable = 'off';
            app.Startoffset_EditField.Position = [166 465 39 23];
            app.Startoffset_EditField.Value = '0.3';
            app.Startoffset_EditField.ValueChangedFcn = createCallbackFcn(app, @plot_new_data, true);
            app.Startoffset_EditField.Tag='Startoffset';

            % Create dxoffsetmonesidedLabel
            app.dxoffsetmonesidedLabel = uilabel(app.LeftPanel);
            app.dxoffsetmonesidedLabel.Enable = 'off';
            app.dxoffsetmonesidedLabel.Position = [17 436 135 22];
            app.dxoffsetmonesidedLabel.Text = 'dx offset [m] (one-sided)';

            % Create dxoffset_EditField
            app.dxoffset_EditField = uieditfield(app.LeftPanel, 'text');
            app.dxoffset_EditField.Enable = 'off';
            app.dxoffset_EditField.Position = [166 437 39 21];
            app.dxoffset_EditField.Value = '0.05';
            app.dxoffset_EditField.ValueChangedFcn = createCallbackFcn(app, @plot_new_data, true);
            app.dxoffset_EditField.Tag='dxoffset';
            
            % Create CorrectT0_Button
            app.CorrectT0_Button = uibutton(app.LeftPanel, 'push');
            app.CorrectT0_Button.Enable = 'off';
            app.CorrectT0_Button.Position = [14 274 192 27];
            app.CorrectT0_Button.Text = 'Pick direct wave and correct t0';
            app.CorrectT0_Button.ButtonPushedFcn = createCallbackFcn(app, @correctT0_cb, true);

            % Create PickApex
            app.PickApex = uibutton(app.LeftPanel, 'push');
            app.PickApex.Enable = 'off';
            app.PickApex.Position = [12 164 194 22];
            app.PickApex.Text = 'Pick apex of hyperbola';
            app.PickApex.ButtonPushedFcn = createCallbackFcn(app, @pickApex_cb, true);

            % Create vcmnsSliderLabel
            app.vcmnsSliderLabel = uilabel(app.LeftPanel);
            app.vcmnsSliderLabel.Enable = 'off';
            app.vcmnsSliderLabel.Position = [14 130 54 22];
            app.vcmnsSliderLabel.Text = 'v [cm/ns]';

            % Create vcmnsSlider
            app.vcmnsSlider = uislider(app.LeftPanel);
            app.vcmnsSlider.Limits = [3 18];
            app.vcmnsSlider.MajorTicks = [4 6 8 10 12 14 16 18];
            app.vcmnsSlider.MinorTicks = [3 5 7 9 11 13 15 17];
            app.vcmnsSlider.Enable = 'off';
            app.vcmnsSlider.Position = [89 139 109 7];
            app.vcmnsSlider.Value = 10;
            app.vcmnsSlider.Tag='v-slider';
            app.vcmnsSlider.ValueChangedFcn = createCallbackFcn(app, @vslider_cb, true);

            % Create t0nsSliderLabel
            app.t0nsSliderLabel = uilabel(app.LeftPanel);
            app.t0nsSliderLabel.Enable = 'off';
            app.t0nsSliderLabel.Position = [13 67 38 22];
            app.t0nsSliderLabel.Text = 't0 [ns]';

            % Create t0nsSlider
            app.t0nsSlider = uislider(app.LeftPanel);
            app.t0nsSlider.MinorTicks = [];
            app.t0nsSlider.Enable = 'off';
            app.t0nsSlider.Position = [72 76 126 7];
            app.t0nsSlider.Value = 20;
            app.t0nsSlider.Tag = 't0-slider';
            app.t0nsSlider.ValueChangedFcn = createCallbackFcn(app, @t0slider_cb, true);

            % Create SavevelocitiesButton
            app.SavevelocitiesButton = uibutton(app.LeftPanel, 'push');
            app.SavevelocitiesButton.Enable = 'off';
            app.SavevelocitiesButton.Position = [10 10 196 23];
            app.SavevelocitiesButton.Text = 'Save velocity';
            app.SavevelocitiesButton.ButtonPushedFcn = createCallbackFcn(app, @save_v_cb, true);

            % Create HyperbolafittingLabel
            app.HyperbolafittingLabel = uilabel(app.LeftPanel);
            app.HyperbolafittingLabel.FontWeight = 'bold';
            app.HyperbolafittingLabel.Enable = 'off';
            app.HyperbolafittingLabel.Position = [15 194 100 22];
            app.HyperbolafittingLabel.Text = 'Hyperbola fitting';

            % Create LoadCMPdataButton
            app.LoadCMPdataButton = uibutton(app.LeftPanel, 'push');
            app.LoadCMPdataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadData_cb, true);
            app.LoadCMPdataButton.Position = [14 547 192 40];
            app.LoadCMPdataButton.Text = 'Load CMP data';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIData
            app.UIData = uiaxes(app.RightPanel);
            title(app.UIData, 'Data')
            xlabel(app.UIData, 'x [m]')
            ylabel(app.UIData, 't [ns]')
            app.UIData.Position = [17 165 200 380];

            % Create UISemblance
            app.UISemblance = uiaxes(app.RightPanel);
            title(app.UISemblance, 'Semblance')
            xlabel(app.UISemblance, 'v [cm/ns]')
            ylabel(app.UISemblance, 't [ns]')
            app.UISemblance.Position = [231 165 200 380];

            % Create SemblanceanalysisPanel
            app.SemblanceanalysisPanel = uipanel(app.RightPanel);
            app.SemblanceanalysisPanel.Enable = 'off';
            app.SemblanceanalysisPanel.Title = 'Semblance analysis';
            app.SemblanceanalysisPanel.FontWeight = 'bold';
            app.SemblanceanalysisPanel.Position = [16 10 413 142];

            % Create CalculateSemblanceButton
            app.calculateSemblanceButton = uibutton(app.SemblanceanalysisPanel);
            app.calculateSemblanceButton.Position = [10 89 130 22];
            app.calculateSemblanceButton.Text = 'Calculate semblance';
            app.calculateSemblanceButton.Enable = 'off';
            app.calculateSemblanceButton.ButtonPushedFcn = createCallbackFcn(app, @semblance_cb, true);

            % Create v_mincmnsEditFieldLabel
            app.v_mincmnsEditFieldLabel = uilabel(app.SemblanceanalysisPanel);
            app.v_mincmnsEditFieldLabel.HorizontalAlignment = 'right';
            app.v_mincmnsEditFieldLabel.Enable = 'off';
            app.v_mincmnsEditFieldLabel.Position = [155 89 79 22];
            app.v_mincmnsEditFieldLabel.Text = 'v_min [cm/ns]';

            % Create vMin_EditField
            app.vMin_EditField = uieditfield(app.SemblanceanalysisPanel, 'text');
            app.vMin_EditField.Enable = 'off';
            app.vMin_EditField.Position = [239 89 30 23];
            app.vMin_EditField.Value = '5';
            app.vMin_EditField.ValueChangedFcn =createCallbackFcn(app, @semblance_cb, true);

            % Create v_maxcmnsEditFieldLabel
            app.v_maxcmnsEditFieldLabel = uilabel(app.SemblanceanalysisPanel);
            app.v_maxcmnsEditFieldLabel.HorizontalAlignment = 'right';
            app.v_maxcmnsEditFieldLabel.Enable = 'off';
            app.v_maxcmnsEditFieldLabel.Position = [286 89 83 22];
            app.v_maxcmnsEditFieldLabel.Text = 'v_max [cm/ns]';

            % Create vMax_EditField
            app.vMax_EditField = uieditfield(app.SemblanceanalysisPanel, 'text');
            app.vMax_EditField.Enable = 'off';
            app.vMax_EditField.Position = [381 89 30 22];
            app.vMax_EditField.Value = '15';
            app.vMax_EditField.ValueChangedFcn =createCallbackFcn(app, @semblance_cb, true);

            % Create AutoPick_Button
            app.AutoPick_Button = uibutton(app.SemblanceanalysisPanel, 'push');
            app.AutoPick_Button.Enable = 'off';
            app.AutoPick_Button.Position = [10 48 150 27];
            app.AutoPick_Button.Text = 'Auto pick 1D-v function';
            app.AutoPick_Button.ButtonPushedFcn = createCallbackFcn(app, @autopick_cb, true);

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.SemblanceanalysisPanel);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.Enable = 'off';
            app.ThresholdEditFieldLabel.Position = [160 51 60 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.SemblanceanalysisPanel, 'text');
            app.ThresholdEditField.Enable = 'off';
            app.ThresholdEditField.Position = [225 49 30 25];
            app.ThresholdEditField.Value='0.4';

            % Create epsEditFieldLabel
            app.epsEditFieldLabel = uilabel(app.SemblanceanalysisPanel);
            app.epsEditFieldLabel.HorizontalAlignment = 'right';
            app.epsEditFieldLabel.Enable = 'off';
            app.epsEditFieldLabel.Position = [245 51 60 22];
            app.epsEditFieldLabel.Text = 'eps';

            % Create epsEditField
            app.epsEditField = uieditfield(app.SemblanceanalysisPanel, 'text');
            app.epsEditField.Enable = 'off';
            app.epsEditField.Position = [310 49 20 25];
            app.epsEditField.Value='3';

            % Create minptsEditFieldLabel
            app.minptsEditFieldLabel = uilabel(app.SemblanceanalysisPanel);
            app.minptsEditFieldLabel.HorizontalAlignment = 'right';
            app.minptsEditFieldLabel.Enable = 'off';
            app.minptsEditFieldLabel.Position = [325 51 60 22];
            app.minptsEditFieldLabel.Text = 'minpts';

            % Create minptsEditField
            app.minptsEditField = uieditfield(app.SemblanceanalysisPanel, 'text');
            app.minptsEditField.Enable = 'off';
            app.minptsEditField.Position = [390 49 20 25];
            app.minptsEditField.Value='5';

            % Create ManualPick_Button
            app.ManualPick_Button = uibutton(app.SemblanceanalysisPanel, 'push');
            app.ManualPick_Button.Enable = 'off';
            app.ManualPick_Button.Position = [13 9 160 27];
            app.ManualPick_Button.Text = 'Manually pick 1D-v-function';
            app.ManualPick_Button.ButtonPushedFcn = createCallbackFcn(app, @manualpick_cb, true);

            % Create SavevelocitiesButton_2
            app.SavevelocitiesButton_2 = uibutton(app.SemblanceanalysisPanel, 'push');
            app.SavevelocitiesButton_2.Enable = 'off';
            app.SavevelocitiesButton_2.Position = [245 9 163 27];
            app.SavevelocitiesButton_2.Text = 'Save velocities';
            app.SavevelocitiesButton_2.ButtonPushedFcn = createCallbackFcn(app, @save_v_S_cb, true);

            % Create HowTo_Label
            app.HowTo_Label = uilabel(app.RightPanel);
            app.HowTo_Label.FontWeight = 'bold';
            app.HowTo_Label.Position = [16 565 800 22];
            app.HowTo_Label.Text = 'Load CMP data and check offsets';
            app.HowTo_Label.WordWrap="on";

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CMPproc_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end