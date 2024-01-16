function [dataz,zmig]=isochrone_mig_3d_varV(data,x,y,t,v,vt,aperture,interp,verbose)

% [dataz,zmig]=isochrone_mig_3d_varV(data,x,y,t,v,vt,aperture,interp,verbose)
%
% Attention: all space units in m!
%
% Isochrone Migration for GPR-Data based on Semicircle superposition,
% using a 1D velocity function
% based on Wilken et al. 2016, but without topography, i.e. antenna
% directivity vertically down
%
% code by Tina Wunderlich, CAU Kiel, 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data,x,y,t: 3D data cube with traces along third dimension, coordinate grids (in m) and time vector
% in ns (t)
% v, vt: v=RMS velocity in m/ns and vt=corresponding time vector
% aperture: Aperture angle in ° (e.g. 30°)
% interp: if =1: interpolation of missing traces (in the lower part this is
% always done due to the aperture of the antenna, but in the upper part
% this interpolation is done by extrapolating the pattern for 1 sample more
% in each direction). if=0, no smoothing/interpolation in the upper part is
% done.
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

% Parameters
dx=x(1,2)-x(1,1); %[m] Trace distance
dy=y(2,1)-y(1,1); %[m] Trace distance
dt=t(2)-t(1); %[ns] sampling interval
ns=length(t);   % number of samples
nx=length(x(1,:)); % number of traces in x-dir
ny=length(y(:,1)); % number of traces in y-dir

dz=dt*min(v(:))/2; % z-grid interval
gzmin=floor(0-(t(end)./2.*max(v(:)))); % min z
gzmax=0; % max z
zmig=[gzmax:-dz:gzmin]; % depth vector in m
nsz=length(zmig); % number of z samples

if verbose==1
    disp(['Number of traces: ',int2str(nx*ny)])
    disp('Start migration...')
end


if all(size(v))==1 % for constant v
    vflag=1; % constant v
    v=v(1);
elseif sum(sum(diff(v,1,2)))==0 % 1d v model
    vflag=2; % 1d v model
    v=interp1(vt,v,t);% 1d-v for all traces
end

% Estimate width of aperture to reduce matrix size later
maxi_z=t(end)./2.*max(v(:)); % max. penetration depth
l=2*pi*maxi_z/360*aperture; % length of semicircle for aperture width and maximum depth
width=ceil(max(maxi_z,l)); % width of area for matrix reduction in m

% calculate radius for each sample on this trace (the same for
% all traces
rad=t./2.*v; % in m

%% make pattern
plusx=round((width/2/dx)); % padding in x dir in samples
plusy=round((width/2/dy));  % padding in y dir in samples
dxgrid=-plusx:plusx; % x-distances to surface point at (0,0) (in indices)
dygrid=-plusy:plusy; % y-distances to surface point (in indices)
dzgrid=[0:nsz-1]'; % z-distances to surface point (in indices)

indexblock=zeros(length(dygrid),length(dxgrid),length(dzgrid));
for i=1:nsz % for each depth step in resulting matrix
    % calculate distance to surface point for each point
    dist=sqrt((dxgrid.*dx).^2+(dygrid'.*dy).^2+((i-1)*dz)^2); % in m
    % set dist larger than max. penetration to nan
    dist(dist>maxi_z)=NaN;
    % set all points out of radiation pattern to 0
    angles=atand(sqrt((dxgrid.*dx).^2+(dygrid'.*dy).^2)./((i-1).*dz)); % angles in degree
    dist(angles>=aperture/2)=NaN; % only angles in radiation pattern
    % find valid points
    valid=find(~isnan(dist));
    for j=1:length(valid)
        [r,c]=ind2sub(size(dist),valid(j));
        indexblock(r,c,i)=round(interp1(rad,[1:ns],dist(valid(j)))); % set sample-index from original trace
    end
    if length(valid)==1 && interp==1 % if only one sample for this time step and interpolation is on... smooth to the neighboring samples
        indexblock(r+1,c,i)=indexblock(r,c,i);
        indexblock(r-1,c,i)=indexblock(r,c,i);
        indexblock(r,c+1,i)=indexblock(r,c,i);
        indexblock(r,c-1,i)=indexblock(r,c,i);
    end
end
% reduce size (remove unused borders)
xind=find(all(permute(indexblock(plusx+1,:,:),[3,2,1])==0)==0,1,'first');
yind=find(all(permute(indexblock(:,plusy+1,:),[3,1,2])==0)==0,1,'first');
indexblock=indexblock(yind:2*plusy-yind+2,xind:2*plusx-xind+2,:);
plusx=plusx-xind+1; % border in x dir (in samples)
plusy=plusy-yind+1; % border in y dir (in samples)

losx=plusx+1; % index for starting of original radargram...
bisx=plusx+nx; % ...index for ending of original radargram
losy=plusy+1; % index for starting of original radargram...
bisy=plusy+ny; % ...index for ending of original radargram

% pad 3D data block (and corresponding matrices) with zeros in the beginning and end to reduce edge
% effects
for i=1:ns
    temp(:,:,i)=[zeros(ny,plusx) data(:,:,i) zeros(ny,plusx)]; % in x dir
    temp2(:,:,i)=[zeros(plusy,nx+2*plusx); temp(:,:,i); zeros(plusy,nx+2*plusx)]; % in y dir
end
data=temp2;

nx2=length(data(1,:,1)); % number of x for padded data
ny2=length(data(:,1,1)); % number of x for padded data

% initialize big matrix
mat=zeros(ny+2*plusy,nx+2*plusx,nsz); % matrix for amplitude summation
count=zeros(ny+2*plusy,nx+2*plusx,nsz); % counter matrix

%% migration
anz=1;
for j=losx:bisx % go through traces (without padded area)
    for i=losy:bisy
        
        % get amplitudes of current trace
        ampl=data(i,j,:);
        
        if ~isnan(ampl(10))
            
            if verbose==1
                if mod(anz,500)==0
                    fprintf('Trace: %d\n',anz)
                end
            end
            
            % set first sample of this trace at surface position
            mat(i,j,1)=data(i,j,1);
            count(i,j,1)=1;
            
          
            % get part of data block with same size as dxgrid/dygrid/dzgrid with
            % current trace as center:
            mattemp=mat(i-plusy:min(i+plusy,ny2),j-plusx:min(j+plusx,nx2),:);
            ctemp=count(i-plusy:min(i+plusy,ny2),j-plusx:min(j+plusx,nx2),:);
            for k=1:numel(mattemp)
                if ~isnan(indexblock(k)) && indexblock(k)~=0
                    mattemp(k)=mattemp(k)+ampl(indexblock(k)); % put amplitudes to right place in small matrix
                    ctemp(k)=ctemp(k)+1; % counter+1
                end
            end
            % replace part in big matrix
            mat(i-plusy:min(i+plusy,ny2),j-plusx:min(j+plusx,nx2),:)=mattemp;
            count(i-plusy:min(i+plusy,ny2),j-plusx:min(j+plusx,nx2),:)=ctemp;
            
        end
        anz=anz+1;
    end
end

% divide amplitudes by number of counts and remove borders
dataz=mat(losy:bisy,losx:bisx,:)./count(losy:bisy,losx:bisx,:);

