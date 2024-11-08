function [traces,t]=reduceNumberOfSamples(traces,t,n)

% [traces,t]=reduceNumberOfSamples(traces,t,n)
%
% Reduce the number of samples by taking only each n-th sample (e.g. if
% n=2: traces_out=traces(1:2:end,:) reduces the sample number by 2)
% Can be also used to increase number of samples by piecewise cubic spline
% interpolation. In this case give e.g. n=0.5 -> double number of samples.
%
% Dr. Tina Wunderlich, CAU Kiel 2022-2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: time vector with same length as columns (in ns)
% n: number of samples to take
%
% Output:
% traces: traces in same order with reduces number of samples
% t: time vector with reduced number of samples

if n>1 % reduce number
    t=t(1:n:end);
    traces=traces(1:n:end,:);
else % increase number of samples
    tnew=interp1(1:length(t),t,1:n:length(t)); % linear for time
    for i=1:size(traces,2)
        temp(:,i)=interp1(t,traces(:,i),tnew,'spline');
    end
    traces=temp;
    t=tnew;
end