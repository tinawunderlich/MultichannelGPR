function [neu]=linearInterpolation_3Dcube(data,x,y,radius)


%%% Linear Interpolation in 2D for every sample-slice in 3D-block
%
% [neu]=linearInterpolation_3Dcube(data,x,y,radius)
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3D-block
% x, y: matrices of same size as data(:,:,1) with coordinates
% radius: radius [m] to be interpolated next to existing data
%
% Output:
% neu: new interpolated 3d-block on same grid as data


neu=data; % initialize interpolated matrix (set given values)

dx=abs(x(1,2)-x(1,1));
n=length(data(1,1,:));     % number of sampleslices

h=waitbar(0,['Linear interpolation is running ...']); % initialize waitbar

for ii=1:n
    if ~mod(ii,10); disp([int2str(ii),'/',int2str(n)]); end

    if ii==1; c1=clock; end

    % make euclidian distance map
    eucmap=chamfer_DT(double(isnan(data(:,:,ii)))).*dx;

    in=find(eucmap<=radius & isnan(data(:,:,ii))); % linear indices of points in radius and with no given value -> interpolation only at these points

    % points with given values
    xy=[x(~isnan(data(:,:,ii))) y(~isnan(data(:,:,ii)))];   % x y for good grid cells
    v=~isnan(data(:,:,ii)); % which values are given?
    val=data(:,:,ii); % data in this slice
    valid=val(v); % valid data points only

    F=scatteredInterpolant(xy(:,1),xy(:,2),valid,'linear');

    [a,b]=ind2sub(size(x),in);

    temp=F(x(in),y(in));
    for j=1:length(temp)
        neu(a(j),b(j),ii)=temp(j);
    end

    if ii==1
        c2=clock;
        diff=minutes(datetime(c2)-datetime(c1));    % time for one run in minutes
    end

    waitbar(ii/n,h,['Approx. ',int2str(diff*n-diff*ii),' minutes remaining']);

end
close(h);
