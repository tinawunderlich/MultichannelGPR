function [traces]=normalize3d(traces,qclip)

% function [traces]=normalize3d(traces,qclip)
%
% Normalize traces of 3D-cube
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: 3D-cube with traces along third dimension
% qclip: cumulative probability in interval [0,1], default is 0.98
%
% Output: 
% traces: normalized traces in 3d cube

if nargin==1
    qclip=0.98;
end

ns=length(traces(1,1,:));

% matrix with quantiles (constant for each trace):
referscale=repmat(quanti(abs(traces),qclip),[1 1 ns]);
% scale data:
traces=traces./referscale;
% replace nans by zeros:
traces(isnan(traces))=0;