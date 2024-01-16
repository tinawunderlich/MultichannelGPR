function [global_coords]=exchange_x_y(global_coords)

% [global_coords]=exchange_x_y(global_coords)
%
% Exchange local coordinates x and y (if you have a left-handed coordinate
% system and want to change it into a right-handed, which is needed for
% helmert transformation)
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% global_coords: original coordinates for each trace (for one profile only)
% 
% Output:
% global_coords: coordinates with exchanged x and y (for one profile only)
%

temp=global_coords(:,1);
global_coords(:,1)=global_coords(:,2);
global_coords(:,2)=temp;
