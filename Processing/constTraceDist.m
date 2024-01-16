function [traces,x,global_coords]=constTraceDist(traces,dist,x,global_coords)

% [traces,x,global_coords]=constTraceDist(traces,dist,x,global_coords)
%
% Make constant trace distance
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% dist: required trace distance in m
% x: profile coordinates in m (length(x)==length(traces(1,:)))
% global_coords: global_coordinates (xyz) for traces
%
% Output:
% traces: traces with constant trace distance
% x/global_coords: new coordinates for each trace with constant trace
% distance
%

ns=length(traces(:,1));
xnew=0:dist:max(x);
% bin traces into new x values
n=discretize(x,xnew); % n is bin number for each trace

% bin traces along new x-vector:
temp=zeros(ns,length(xnew));
gc=zeros(length(xnew),length(global_coords(1,:)));
for i=1:length(xnew)
    if sum(n==i)>0 % if there are traces in this bin
        temp(:,i)=mean(traces(:,n==i),2);
        gc(i,:)=mean(global_coords(n==i,:),1);
    else
        temp(:,i)=NaN(ns,1);
    end
end

% interpolate coordinates for missing traces:
a=isnan(temp(1,:));
gc(a,1)=interp1(xnew(~a),gc(~a,1),xnew(a),'linear','extrap');
gc(a,2)=interp1(xnew(~a),gc(~a,2),xnew(a),'linear','extrap');
if length(global_coords(1,:))==3
    gc(a,3)=interp1(xnew(~a),gc(~a,3),xnew(a),'linear','extrap');
end


% replace x and global_coords and traces:
x=xnew;
global_coords=gc;
traces=temp;