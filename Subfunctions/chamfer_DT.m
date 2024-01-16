function [dist] = chamfer_DT(im)
%% Chamfer distance transform using a 3x3 chamfer mask (=approximation of euclidian distance transform)
%
% Dr.Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% [dist] = chamfer_DT(im)
%
% Input:
% im: Input image: (non-feature,feature) = (1,0) 
%
% Output:
% dist: image with same size as input image with distance to nearest
% neighbor of each pixel

% weights for neighboring cells for approximation of euclidian distance
a = 3;   
b = 4;

v=im([1,1:end,end],[1,1:end,end]); % pad matrix
v(v~=0) = inf; % Replace 1's by a suitably large number

rows = size(v,1); % Number of lines (rows)
cols = size(v,2); % Number of columns


%%% Forward run:
for k1=2:rows-1
    for k2=2:cols-1
           v(k1,k2) = min([v(k1-1,k2-1)+b, v(k1-1,k2)+a,... 
           v(k1-1,k2+1)+b, v(k1,k2-1)+a, v(k1,k2)]);
    end
end
forward = v(2:end-1,2:end-1);

%%% backward run:
for k1=rows-1:-1:2
    for k2=cols-1:-1:2
        v(k1,k2) = min([v(k1,k2), v(k1,k2+1)+a, v(k1+1,k2-1)+b,...
                       v(k1+1,k2)+a, v(k1+1,k2+1)+b]);
    end
end
backward = v(2:end-1,2:end-1);


% find minimum and divide result by 4
dist=min(forward, backward)./4;
