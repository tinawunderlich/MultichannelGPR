function [neu]=idw3dblock_topo(data,x,y,r,p)


%%% IDW in 2D for every sample-slice in 3D-block with topography
%%% corrected/migrated data (=every slice has different sample points)
%
% [neu]=idw3dblock_topo(data,x,y,r,p)
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3D-block
% x,y: unknown points (grids with same size)
% r: Radius for interpolation between neighboring points in m
% p: power of inverse distance weighting
%
% Output:
% neu: new interpolated 3d-block on same grid as data


neu=data; % initialize interpolated matrix (set given values)

dx=abs(x(1,2)-x(1,1));
n=length(data(1,1,:));     % number of sampleslices

h=waitbar(0,['Inverse Distance weighting interpolation is running ...']); % initialize waitbar

for ii=1:n
    if ~mod(ii,10); disp([int2str(ii),'/',int2str(n)]); end

    if ii==1; c1=clock; end

    % make euclidian distance map
    eucmap=chamfer_DT(double(isnan(data(:,:,ii)))).*dx;

    in=find(eucmap<=r & isnan(data(:,:,ii))); % linear indices of points in radius and with no given value -> interpolation only at these points

    % points with given values
    xy=[x(~isnan(data(:,:,ii))) y(~isnan(data(:,:,ii)))];   % x y for good grid cells
    v=~isnan(data(:,:,ii)); % which values are given?
    val=data(:,:,ii); % data in this slice
    valid=val(v); % valid data points only


    for i=1:length(in)

        dist=sqrt((xy(:,1)-x(in(i))).^2+(xy(:,2)-y(in(i))).^2);   % distance to all points

        temp=zeros(size(data(1,1,:)));
        w=(dist<=r);
        if any(w)
            d=dist(dist<=r).^p;
            s=sum(1./dist(dist<=r).^p);
            temp=sum(valid(w)./d)/s;
            [a,b]=ind2sub(size(x),in(i));
            neu(a,b,ii)=temp;
        end

    end

    if ii==1
        c2=clock;
        diff=minutes(datetime(c2)-datetime(c1));    % time for one run in minutes
    end

    waitbar(ii/n,h,['Approx. ',int2str(diff*n-diff*ii),' minutes remaining']);

end
close(h);
