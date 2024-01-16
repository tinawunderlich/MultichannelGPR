function [traces,t]=reduceNumberOfSamples(traces,t,n)

% [traces,t]=reduceNumberOfSamples(traces,t,n)
%
% Reduce the number of samples by taking only each n-th sample (e.g. if
% n=2: traces_out=traces(1:2:end,:) reduces the sample number by 2)
%
% Dr. Tina Wunderlich, CAU Kiel 2022, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: time vector with same length as columns (in ns)
% n: number of samples to take
%
% Output:
% traces: traces in same order with reduces number of samples
% t: time vector with reduced number of samples

t=t(1:n:end);
traces=traces(1:n:end,:);