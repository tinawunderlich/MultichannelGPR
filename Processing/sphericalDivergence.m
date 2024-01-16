function [traces]=sphericalDivergence(traces,t)

% [traces]=sphericalDivergence(traces,t)
%
% correct for spherical divergence by multiplication of each trace with t
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces in columns
% t: corresponding time vector in ns
%
% Output:
% traces: Matrix with corrected traces in same order as traces



% check if t is column vector
if length(t(1,:))>1
    t=transpose(t);
end
% Check if t in in ns
if max(t)<1
    t=t.*1e9; % now in ns
end

traces=traces.*repmat(t,[1 length(traces(1,:))]);