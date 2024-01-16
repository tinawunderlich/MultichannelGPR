function [dataz,z]=isochrone_mig_2d_varV(data,x,t,v,aperture,verbose)

% [dataz,z]=isochrone_mig_2d_varV(data,x,t,v,aperture,verbose)
%
% Attention: all space units in m!
%
% Isochrone Migration for GPR-Data based on Semicircle superposition,
% using a 2D velocity function
% based on Wilken et al. 2016, but without topography, i.e. antenna
% directivity is vertically down
%
% code by Tina Wunderlich, CAU Kiel, 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data,x,t: Radartraces with coordinates along profile (x in m) and time vector
% in ns (t) (x with constant trace distance in m)
% v: RMS velocity in m/ns either as grid with same size as data (vertically and horizontally variable) or constant
% value
% aperture: Aperture angle in ° (e.g. 30°)
% verbose: display progress on=1 or off=0
%
% Output:
% dataz, z: migrated data and corresponding depth-vector (z in m)


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
dt=t(2)-t(1); %[ns] sampling interval
ns=length(t);   % number of samples
nt=length(x); % number of traces

dz=dt*min(v(:))/2; % z-grid interval
gzmin=floor(0-(t(end)./2.*max(v(:)))); % min z
gzmax=0; % max z
z=[gzmax:-dz:gzmin+dz]; % depth vector in m
nsz=length(z); % number of z samples

if verbose==1
    disp(['Number of traces: ',int2str(nt)])
    disp('Start migration...')
end


if all(size(v))==1 % for constant v -> make vgrid
    v=zeros(size(data))+v;
    vflag=1; % constant v
elseif sum(sum(diff(v,1,2)))==0 % 1d v model
    vflag=2; % 1d v model
else
    vflag=3; % 2d v model
end

% Estimate width of aperture to reduce matrix size later
maxi_z=t(end)./2.*max(v(:)); % max. penetration depth
l=2*pi*maxi_z/360*aperture; % length of semicircle for aperture width and maximum depth
width=ceil(max(maxi_z,l)); % width of area for matrix reduction in m

% pad radargram (and corresponding matrices) with zeros in the beginning and end to reduce edge
% effects
data=[zeros(length(t),round((width/2/dx))) data zeros(length(t),round((width/2/dx)))];
v=[zeros(length(t),round((width/2/dx))) v zeros(length(t),round((width/2/dx)))];
x=[zeros(1,round((width/2/dx))) x zeros(1,round((width/2/dx)))];
los=find(data(10,:)~=0,1,'first'); % index for starting of original radargram...
bis=find(data(10,:)~=0,1,'last'); % ...index for ending of original radargram

nt2=length(x); % number of traces for padded radargram

% initialize big matrix
mat=zeros(nsz,nt2); % matrix for amplitude summation
count=zeros(nsz,nt2); % counter matrix

if vflag==3 % 2d v model
    for j=los:bis-1 % go through traces (without padded area)
        if verbose==1
            if mod(j-los+1,100)==0
                fprintf('Trace: %d\n',j-los+1)
            end
        end
        
        % Angle of radiation (0° = vertically down)
        alpha=0;
        
        % calculate radius for each Sample on this trace (depending on v
        % in this trace)
        rad=t./2.*v(:,j); % in m
        
        dxgrid=[0:dx:dx*(nt2-1)]-(j-1)*dx; % x-distances to surface point
        dzgrid=[0:dz:dz*(nsz-1)]'; % z-distances to surface point
        
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
        % set all points below last sample to 0
        distw(distw>max(rad))=0;
        % delete columns with only zeros in distw and reset s and e
        colzero=all(distw==0,1);
        e=s+find(colzero==0,1,'last'); % new end index
        s=s+find(colzero==0,1,'first')-1; % new start index
        distw(:,colzero)=[];
        
        % set first sample of this trace at surface position
        mat(1,j)=data(1,j);
        count(1,j)=1;
        
        % get all possible distances in radiation pattern
        [un,~,b]=unique(distw);
        b(:,2)=[1:length(b)]';  % linear index in matrix distw
        % first column in b is number of un -> first value in un is 0 ->
        % delete all rows with b(:,1)==1
        b(b(:,1)==1,:)=[];
        
        % get bins along trace with edges between samples:
        grenzen=[rad(1:end-1)+diff(rad)/2; rad(end)];
        
        % get number of distances in these bins -> determine how many points
        % along semicircle get part of this amplitude
        [~,~,bins]=histcounts(un,grenzen); % anz= number of points per bin
        
        % get amplitudes of current trace
        ampl=data(2:end,j);
        
        % set parts of amplitudes at right place:
        wb=bins(b(:,1)); % in which bin lies current value of b?
        wb(wb==0)=1; % if somewhere 0 -> set 1
        gA=ampl(wb); % parts of amplitude for all points
        mat_temp=zeros(nsz,e-s+1); % initialise small area of matrix so that linear indices of b fit
        mat_temp(b(:,2))=gA; % put amplitudes to right place in small matrix
        mat(:,s:e)=mat(:,s:e)+mat_temp; % put small matrix into big matrix
        count(:,s:e)=count(:,s:e)+(mat_temp~=0); % counter+1
    end
    
else % for constant or 1d v model
    % calculate radius for each sample on this trace (the same for
    % all traces
    rad=t./2.*v(:,los); % in m
    
    for j=los:bis-1 % go through traces (without padded area)
        if verbose==1
            if mod(j-los+1,100)==0
                fprintf('Trace: %d\n',j-los+1)
            end
        end
        
        % Angle of radiation (0° = vertically down)
        alpha=0;
        
        dxgrid=[0:dx:dx*(nt2-1)]-(j-1)*dx; % x-distances to surface point
        dzgrid=[0:dz:dz*(nsz-1)]'; % z-distances to surface point
        
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
        % set all points below last sample to 0
        distw(distw>max(rad))=0;
        % delete columns with only zeros in distw and reset s and e
        colzero=all(distw==0,1);
        e=s+find(colzero==0,1,'last'); % new end index
        s=s+find(colzero==0,1,'first')-1; % new start index
        distw(:,colzero)=[];
        
        % set first sample of this trace at surface position
        mat(1,j)=data(1,j);
        count(1,j)=1;
        
        if j==los % only for first trace, then stays constant
            % get all possible distances in radiation pattern
            [un,~,b]=unique(distw);
            b(:,2)=[1:length(b)]';  % linear index in matrix distw
            % first column in b is number of un -> first value in un is 0 ->
            % delete all rows with b(:,1)==1
            b(b(:,1)==1,:)=[];
            
            % get bins along trace with edges between samples:
            grenzen=[rad(1:end-1)+diff(rad)/2; rad(end)];
            
            % get number of distances in these bins -> determine how many points
            % along semicircle get part of this amplitude
            [~,~,bins]=histcounts(un,grenzen); % anz= number of points per bin
        end
        
        % get amplitudes of current trace
        ampl=data(2:end,j);
        
        % set parts of amplitudes at right place:
        wb=bins(b(:,1)); % in which bin lies current value of b?
        wb(wb==0)=1; % if somewhere 0 -> set 1
        gA=ampl(wb); % parts of amplitude for all points
        mat_temp=zeros(nsz,e-s+1); % initialise small area of matrix so that linear indices of b fit
        mat_temp(b(:,2))=gA; % put amplitudes to right place in small matrix
        mat(:,s:e)=mat(:,s:e)+mat_temp; % put small matrix into big matrix
        count(:,s:e)=count(:,s:e)+(mat_temp~=0); % counter+1
    end
end

% divide amplitudes by number of counts
dataz=mat./count;

% last trace has not been migrated -> replace with the last but one
% trace
dataz(:,end)=dataz(:,end-1);

% delete padded areas
dataz=dataz(:,los:bis);
end


