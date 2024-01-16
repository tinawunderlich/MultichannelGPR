function [traces]=normalize2d(traces,qclip)

% Normalize traces of 2D-matrix [traces]=normalize2d(traces,qclip)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: 2d matrix with traces along columns
% qclip: cumulative probability in interval [0,1], default is 0.98
%
% Output:
% traces: normalized traces in matrix

if nargin==1
    qclip=0.98;
end

ns=length(traces(:,1));

% matrix with quantiles (constant for each trace):
%referscale=repmat(quantile(abs(traces),qclip,1),[ns 1]);
referscale=repmat(quanti(abs(traces),qclip),[ns 1]);
% scale data:
traces=traces./referscale;
% replace nans by zeros:
traces(isnan(traces))=0;