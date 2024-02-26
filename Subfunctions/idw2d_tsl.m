function [neu]=idw2d_tsl(data,x,y,eucmap,r,p)


%%% IDW in 2D for all timeslices
%
% [neu]=idw2d_tsl(data,x,y,eucmap,r,p)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: Timeslices: Cell array mit Matrizen mit bekannten Werten (gleiche Größe wie x, y), Fehlstellen
% sind NaN
% x,y: unknown points (grids with same size)
% eucmap: approximated euclidian distance map (Distance to next neighbor in bins)
% r: Radius innerhalb dessen die Punkte zum Interpolieren genutzt werden
% sollen (in bins)
% p: Potenz
%
% Output:
% neu: Interpolierte Werte an den Stellen wie in den x-y-Grids (auch als Grids in Cell array)

% for topo (=one layer):
if ~iscell(data)
    temp=data;
    data=cell(1);
    data{1}=temp; % put data in cell
end
    

in=find(eucmap<=r & isnan(data{1})); % indices of points in radius and with no given value -> interpolation only at these points

n=length(data);     % number of timeslices

% points with given values
xy=[x(~isnan(data{1}(:))) y(~isnan(data{1}(:))) zeros(length(x(~isnan(data{1}(:)))),n)];   % x y value_tsl1 value_tsl2 ...
for i=1:n
    xy(:,i+2)=data{i}(~isnan(data{i}(:))); % use given values
    
    neu{i}=data{i};   % initialize interpolated matrix (set given values)
end


h=waitbar(0,'Inverse Distance weighting interpolation is running...'); % initialize waitbar

for i=1:length(in)
    if i==1; c1=clock; end
    
    dist=sqrt((xy(:,1)-x(in(i))).^2+(xy(:,2)-y(in(i))).^2);   % distance to all points
    w=dist<=r; % these points will be used for interpolation

    for j=1:n
        neu{j}(in(i))=sum(xy(w,j+2)./dist(w).^p)/sum(1./dist(w).^p);
    end
    
    if i==1
        c2=clock;
        diff=minutes(datetime(c2)-datetime(c1));    % time for one run in minutes
    end
    
    waitbar(i/length(in),h,['Approx. ',int2str(diff*length(in)-diff*i),' minutes remaining']);
end
close(h);

