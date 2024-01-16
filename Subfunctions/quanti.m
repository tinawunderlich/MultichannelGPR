function [q]=quanti(x,p)

% calculates the p-quantil from values in every column of x-matrix (2D) or for every trace in 3Dmatrix
% (similar to function quantile from MATLAB-Statistics toolbox)

if length(size(x))==2 % 2D
    x=sort(x,1);
    n=length(x(:,1));
    if p~=0 && p~=1
        q=x(floor(n*p+1),:);
    else
        q=x(end,:);
    end
else  % 3D
    x=sort(x,3);
    n=length(x(1,1,:));
    if p~=0 && p~=1
        q=x(:,:,floor(n*p+1));
    else
        q=x(:,:,end);
    end
end