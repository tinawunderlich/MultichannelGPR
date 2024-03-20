function [traces,x,global_coords]=turnProfiles(traces,x,global_coords)

% [traces]=turnProfiles(traces)
%
% Turn profiles so that they run all from west to east
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% global_coords: original coordinates for each trace
% x: profile coordinates starting from 0 to length(profile)
% 
% Output:
% traces: traces that are turned around
% global_coords: coordinates that are turned around
% x: new profile coordinates starting from 0 to length(profile)
%

if global_coords(1,1)>global_coords(end,1)
    global_coords=flipud(global_coords);
    traces=fliplr(traces);
    dist=fliplr(diff(x));
    x=[0 reshape(cumsum(dist),[1 length(dist)])];
end