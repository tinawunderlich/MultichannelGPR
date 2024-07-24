function data_out=scanline(data,xgrid,ygrid,polygon,ind)

%%% Scanline-Algorithm for filling polygons
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 2D grid with data points with corresponding coordinate grids
% xgrid/ygrid
% polygon: points of polygon in same coordinate system as xgrid/ygrid (size [n,2])
% ind: number to fill in the polygon


data_out=zeros(size(data)); % initialize output data

%%% append first point to the end
polygon=[polygon; polygon(1,:)];

%%% go through every row of the grid
for i=1:length(ygrid(:,1))
    
    xs{i}=[];   % Array for intersection points
    
    % go through neighboring points in polygon and check if they are over and above the current row
    for j=1:length(polygon(:,1))-1
        if (polygon(j,2)<=ygrid(i,1) && polygon(j+1,2)>=ygrid(i,1)) || (polygon(j,2)>=ygrid(i,1) && polygon(j+1,2)<=ygrid(i,1))
            % if yes, make a linear polyfit through points and determine intersection point with current row
            p=polyfit([polygon(j,1) polygon(j+1,1)],[polygon(j,2) polygon(j+1,2)],1);
            % intersection point
            xs{i}=[xs{i}; (ygrid(i,1)-p(2))/p(1)];
        end
    end
    
    % sort intersection points
    xs_sort{i}=sort(xs{i});
    
    if ~isempty(xs_sort{i}) % if there are intersection points in this row...      
        for k=1:2:length(xs_sort{i})-1
            % go through every x-value in this row
            for j=1:length(xgrid(1,:))
                if xgrid(1,j)>=xs_sort{i}(k) && xgrid(1,j)<=xs_sort{i}(k+1)
                    data_out(i,j)=ind; % set number in this grid position
                end
            end
        end
    end
end
