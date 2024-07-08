function [dataz,z]=topomig2d_varV(data,x,t,topo,v,aperture,flag,verbose,zmin,zmax)

% [dataz,z]=topomig2d_varV(data,x,t,topo,v,aperture,flag,verbose,zmin,zmax)
%
% Attention: all space units in m!
%
% Topographic Migration for GPR-Data based on Semicircle superposition,
% using a 2D velocity function
% based on Allroggen et al. 2014 and Wilken et al. 2016
%
% code by Tina Wunderlich, CAU Kiel, 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data,x,t: Radartraces with coordinates along profile (x in m) and time vector
% in ns (t) (x with constant trace distance in m!!)
% topo: height values along profile for each x (in m)
% v: RMS velocity in m/ns either as grid with same size as data (vertically and horizontally variable) or constant
% value
% aperture: Aperture angle in ° (e.g. 30°)
% flag: =1: Topocorrection, =2: Topomigration
% verbose: display progress on=1 or off=0
% zmin: minimum of vector z. If zmin shall be determined automatically, give [], otherwise zmin
% in m
% zmax: maximal height (should be above maximum height of profile, rounded)
% in m. If zmax shall be determined automatically, give [].
%
% Output:
% dataz, z: migrated data and corresponding depth-vector (z in m) below zmax


% if necessary, convert v [cm/ns] to [m/ns]
if all(v>=3) % if given in cm/ns:
    v=v./100;    % convert to m/ns (which is required by the following code)
end

if length(t(1,:))>length(t(:,1))
    t=t'; % column vector
end
if length(x(1,:))<length(x(:,1))
    x=x'; % row vector
end

% Parameters
dx=x(2)-x(1); %[m] Trace distance
if round(dx*1000)~=round(mean(diff(x))*1000)
    disp('!!Make constant trace distance before migration! Aborting!')
    return;
end
dt=t(2)-t(1); %[ns] sampling interval
ns=length(t);   % number of samples
nt=length(x); % number of traces

dz=dt*min(v(:))/2; % z-grid interval
if isempty(zmin)
    gzmin=floor(min(topo)-(t(end)./2.*max(v(:))));
else
    gzmin=zmin;
end
if isempty(zmax)
    gzmax=ceil(max(topo));
else
    gzmax=zmax;
end
z=[gzmax:-dz:gzmin+dz]; % depth vector in m
nsz=length(z); % number of z samples

if verbose==1
    disp(['Number of traces: ',int2str(nt)])
end


%Standard topographic correction:
if flag==1
    if verbose==1
        disp('Topographic correction')
    end

    % initialize data matrix
    dataz=NaN(nsz,nt);

    % for constant v:
    if all(size(v)==1)
        for j=1:nt % go through all traces
            dataz(round((gzmax-topo(j))/dz):round((gzmax-topo(j))/dz)+ns-1,j)=data(:,j);
        end
        % for variable v:
    else
        for j=1:nt % go through all traces
            ztemp=t./2.*v(:,j); % use v at this trace for depth conversion
            % interpolate amplitudes from trace to general zgrid
            neu=0:dz:max(ztemp);
            dataz(round((gzmax-topo(j))/dz):round((gzmax-topo(j))/dz)+length(neu)-1,j)=interp1(ztemp,data(:,j),neu);
        end
    end

else % Topomigration
    if verbose==1
        disp('Topo migration')
    end

    if length(topo(:,1))<length(topo(1,:))
        topo=topo'; % make topo a column vector
    end

    if all(size(v)==1) % for constant v -> make vgrid
        v=zeros(size(data))+v;
    end

    % Estimate width of aperture to reduce matrix size later
    maxi_z=t(end)./2.*max(v(:)); % max. penetration depth [m]
    l=2*pi*maxi_z/360*aperture; % length of semicircle for aperture width and maximum depth [m]
    width=ceil(max(maxi_z,l)).*2; % width of area for matrix reduction [m]

    % if there are missing (NaN) traces inbetween -> split radargram here
    ind=find(isnan(data(1,:))); % find NaN-Traces
    diffind=diff(ind); % find tracenumber-differences between NaN-traces
    if ~isempty(ind) % there are nan-traces
        di=find(diffind~=1);    % jumps
        part=[1 ind(1)-1]; % start-end of part with data
        for i=1:length(di)
            part=[part; ind(di(i))+1 ind(di(i)+1)-1];
        end
        part=[part; ind(end)+1 nt];
    else
        part=[1 nt]; % all traces
    end

    dataz=NaN(nsz,nt); % initialize big matrix

    % make topomigration for all parts of radargram
    for k=1:length(part(:,1))
        % pad radargram (and corresponding matrices) with zeros in the beginning and end to reduce edge
        % effects
        data_p=[zeros(length(t),round((width/2/dx))) data(:,part(k,1):part(k,2)) zeros(length(t),round((width/2/dx)))];
        v_p=[zeros(length(t),round((width/2/dx))) v(:,part(k,1):part(k,2)) zeros(length(t),round((width/2/dx)))];
        x_p=[zeros(1,round((width/2/dx))) x(part(k,1):part(k,2)) zeros(1,round((width/2/dx)))];
        topo_p=[zeros(1,round((width/2/dx))) topo(part(k,1):part(k,2))' zeros(1,round((width/2/dx)))];
        los=find(data_p(10,:)~=0,1,'first'); % index for starting of original radargram...
        bis=find(data_p(10,:)~=0,1,'last'); % ...index for ending of original radargram

        nt2=length(x_p); % number of traces for padded radargram

        % initialize big matrix
        mat=zeros(nsz,nt2); % matrix for amplitude summation
        count=zeros(nsz,nt2); % counter matrix

        for j=los:bis-1 % go through traces (without padded area)
            if verbose==1
                if mod(j-los+part(k,1),100)==0
                    fprintf('Trace: %d\n',j-los+part(k,1))
                end
            end

            % current gradient of topo
            dtopo=topo_p(j+1)-topo_p(j);
            % Angle of radiation (0° = vertically down)
            if dtopo<0
                alpha=atand(dtopo/dx);
            elseif dtopo==0
                alpha=0;
            else
                alpha=atand(dtopo/dx);
            end

            % calculate radius for each Sample on this trace (depending on v
            % in this trace)
            rad=t./2.*v_p(:,j); % in m

            ind_z=find(abs(z-topo_p(j))==min(abs(z-topo_p(j))),1,'first'); % get z-index of surface point in this trace
            dxgrid=[0:dx:dx*(nt2-1)]-(j-1)*dx; % x-distances to this point
            dzgrid=[0:dz:dz*(nsz-1)]'-(ind_z-1)*dz; % z-distances to this point

            % reduce width of dxgrid (to reduce runtime)
            s=find(dxgrid<=-width/2,1,'last'); % start index
            if isempty(s); s=1; end
            e=find(dxgrid>=width/2,1,'first'); % end index
            if isempty(e); e=length(dxgrid); end
            % only use this part for calculation
            dist=sqrt(dxgrid(s:e).^2+dzgrid.^2); % distances to surface point on this trace

            % set all points out of radiation pattern to 0
            angles=atand(dxgrid(s:e)./dzgrid); % angles in degree
            distw=dist.*(angles<=alpha+aperture/2 & angles>=alpha-aperture/2); % only angles in radiation pattern
            % set all angles above topo to 0
            distw(1:ind_z-1,:)=0;
            % set all points below last sample to 0
            distw(distw>max(rad))=0;
            % delete columns with only zeros in distw and reset s and e
            colzero=all(distw==0,1);
            e=s+find(colzero==0,1,'last'); % new end index
            s=s+find(colzero==0,1,'first')-1; % new start index
            distw(:,colzero)=[];

            % set first sample of this trace at surface position
            mat(ind_z,j)=data_p(1,j);
            count(ind_z,j)=1;

            % get all possible distances in radiation pattern
            [un,~,b]=unique(distw);
            b(:,2)=[1:length(b)]';  % linear index in matrix distw
            % first column in b is number of un -> first value in un is 0 ->
            % delete all rows with b(:,1)==1
            b(b(:,1)==1,:)=[];

            % get bins along trace with edges between samples:
            grenzen=[rad(1:end-1)+diff(rad)/2; rad(end)];
            % force grenzen to be monotonically increasing:
            for i=2:length(grenzen)
                if grenzen(i)<=grenzen(i-1)
                    grenzen(i)=grenzen(i-1)+mean(diff(grenzen));
                end
            end

            % get number of distances in these bins -> determine how many points
            % along semicircle get part of this amplitude
            [~,~,bins]=histcounts(un,grenzen); % anz= number of points per bin

            % get amplitudes of current trace
            ampl=data_p(2:end,j);

            % set parts of amplitudes at right place:
            wb=bins(b(:,1)); % in which bin lies current value of b?
            wb(wb==0)=1; % if somewhere 0 -> set 1
            gA=ampl(wb); % parts of amplitude for all points
            mat_temp=zeros(nsz,e-s+1); % initialise small area of matrix so that linear indices of b fit
            mat_temp(b(:,2))=gA; % put amplitudes to right place in small matrix
            mat(:,s:e)=mat(:,s:e)+mat_temp; % put small matrix into big matrix
            count(:,s:e)=count(:,s:e)+(mat_temp~=0); % counter+1
        end

        % divide amplitudes by number of counts
        dataz_p=mat./count;

        % last trace has not been migrated -> replace with the last but one
        % trace
        dataz_p(:,end)=dataz_p(:,end-1);

        % set part of radargram in big matrix and delete padded areas
        dataz(:,part(k,1):part(k,2))=dataz_p(:,los:bis);
    end
end



