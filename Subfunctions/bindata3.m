function [zout,xgrid,ygrid] = bindata3(z,x,y,xrg,yrg)

% [zout] = bindata3(z,x,y,xrg,yrg)
%
% based on 2d-code by Patrick Mineault
%    %Refs: https://xcorr.net/?p=3326
%    %      http://www-pord.ucsd.edu/~matlab/bin.htm
% modified/commented by Tina Wunderlich, CAU Kiel 2020-2024, tina.wunderlich@ifg.uni-kiel.de
%
% x, y: vectors of length N with x/y coordinates
% z: traces of size (ns,N) corresponding to x/y
% xrg, yrg: vectors with edges of bins in x and y direction
% Output: zm: Binned grid of size(length(yrg)-1,length(xrg)-1) with time along the third dimension
%         xgrid, ygrid: coordinate grids corresponding to zm (midpoints of bins)
%
% This code takes advantage of finding the bins only once and them applying
% it to all time-sample-slices

    numsamples=length(z(:,1)); % number of samples

    dx=xrg(2)-xrg(1);
    dy=yrg(2)-yrg(1);

    a=xrg+dx/2;
    b=yrg+dy/2;
    [xgrid,ygrid]=meshgrid(a(1:end-1),b(1:end-1));
    [r,c]=size(xgrid);

    % delete all values outside xrg/yrg
    weg=(x>max(xrg) | x<min(xrg) | y>max(yrg) | y<min(yrg));
    x(weg)=[];
    y(weg)=[];
    z(:,weg)=[];

    % finding the x and y bins for each coordinate
    [~,whichbinx] = histc(x,xrg);   % returns the bin number in x direction (left edge is included in bin)
    [~,whichbiny] = histc(y,yrg);   % returns the bin number in y direction (left edge is included in bin)
 
    % corrected bin numbers (if there is a coordinate == last right bin edge it is set into last bin)
    binsx = min(max(whichbinx,1),c);    
    binsy = min(max(whichbiny,1),r);
 
    bins = (binsy-1).*(c)+binsx; % linear index of bin in matrix
 
    xpos = ones(size(bins,1),1);
    ns = sparse(bins,xpos,1,c*r,1); % sparse matrix with number of counts in each bin
    zout=zeros(r,c,numsamples,'single');
    for i=1:numsamples
        ysum = sparse(bins,xpos,double(z(i,:)),r*c,1);   % sparse matrix with sum of z values in each bin
        zm = full(ysum)./(full(ns));    % vector with mean values (sum/number) for each bin in vector bins (if no hit, then NaN)
        zout(:,:,i) = reshape(zm,c,r)';  % output matrix
    end

end
