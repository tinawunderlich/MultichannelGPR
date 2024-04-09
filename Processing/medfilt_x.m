function [traces]=medfilt_x(traces,t,m,tstart)

% [traces]=medfilt_x(traces,t,m,tstart)
%
% Apply median filter along x over m samples (only for part below tstart)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: corresponding time vector for traces
% m: number of samples for moving window
% tstart: will be applied only for t>=tstart
%
% Output:
% traces: filtered traces 
%

% x-wise median filter
% only deeper part:
temp=traces(t>=tstart,:);
traces(t>=tstart,:)=medfilt(temp',m)';