function [neu]=idw3dblock(data,x,y,r,p)


%%% IDW in 2D for every sample-slice in 3D-block
%
% [neu]=idw3dblock(data,x,y,r,p)
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3D-block
% x,y: unknown points (grids with same size)
% r: Radius for interpolation between neighboring points in m
% p: power of inverse distance weighting
%
% Output:
% neu: new interpolated 3d-block on same grid as data


% make euclidian distance map
dx=abs(x(1,2)-x(1,1));
eucmap=chamfer_DT(double(isnan(data(:,:,1)))).*dx;

in=find(eucmap<=r & isnan(data(:,:,1))); % linear indices of points in radius and with no given value -> interpolation only at these points

n=length(data(1,1,:));     % number of sampleslices

% points with given values
xy=[x(~isnan(data(:,:,1))) y(~isnan(data(:,:,1)))];   % x y for good grid cells
v=~isnan(data(:,:,1)); % which values are given?
for j=1:n
    val=data(:,:,j); % data in this slice
    valid{j}=val(v); % valid data points only
end

neu=data; % initialize interpolated matrix (set given values)

h=waitbar(0,'Inverse Distance weighting interpolation is running...'); % initialize waitbar

for i=1:length(in)
    if i==1; c1=clock; end
    
    dist=sqrt((xy(:,1)-x(in(i))).^2+(xy(:,2)-y(in(i))).^2);   % distance to all points
    
    temp=zeros(size(data(1,1,:)));
    w=(dist<=r);
    if any(w)
        d=dist(dist<=r).^p;
        s=sum(1./dist(dist<=r).^p);
        for j=1:n
            temp(1,1,j)=sum(valid{j}(w)./d)/s;
        end
        [a,b]=ind2sub(size(x),in(i));
        neu(a,b,:)=temp;
    end
    
    if i==1
        c2=clock;
        diff=minutes(datetime(c2)-datetime(c1));    % time for one run in minutes
    end
    
    waitbar(i/length(in),h,['Approx. ',int2str(diff*length(in)-diff*i),' minutes remaining']);
end
close(h);

