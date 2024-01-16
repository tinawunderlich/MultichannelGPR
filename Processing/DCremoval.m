function [traces]=DCremoval(traces,t,start,ende)

% [traces]=DCremoval(traces,t,start,ende)
%
% removal of amplitude offsets of single traces: calculates mean of trace
% between start and ende and subtracts mean from trace (individually for
% each trace)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: time vector with same length as columns (in ns)
% start, ende: amplitude offset is determined in interval [start ende] (in
% ns)
% if no start&ende is given, the whole trace will be used
%
% Output:
% traces: traces in same order with subtracted amplitude offset



if nargin==2
    start=0;
    ende=t(end);
elseif nargin==3
    ende=t(end);
end
if ende>t(end)
    ende=t(end);
end

traces=traces-repmat(mean(traces(t>=start & t<=ende,:)),[length(t) 1]);