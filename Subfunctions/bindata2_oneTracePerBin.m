function [zm,xgrid,ygrid,ind,vals] = bindata2_oneTracePerBin(z,x,y,xrg,yrg)

% [zm] = bindata2_oneTracePerBin(z,x,y,xrg,yrg)
% Bins data onto a regular grid, but just puts the first occuring trace per
% bin in it, no mean value!
%
% based on code by Patrick Mineault
%    %Refs: https://xcorr.net/?p=3326
%    %      http://www-pord.ucsd.edu/~matlab/bin.htm
% modified/commented by Tina Wunderlich, CAU Kiel 2020-2022, tina.wunderlich@ifg.uni-kiel.de
%
% x, y, z: vectors of same length with x/y coordinates and amplitude values z
% xrg, yrg: vectors with edges of bins in x and y direction
% Output: zm: Binned grid of size(length(yrg)-1,length(xrg)-1)
%         xgrid, ygrid: coordinate grids corresponding to zm (midpoints of
%         bins)


    dx=xrg(2)-xrg(1);
    dy=yrg(2)-yrg(1);

    % delete all values outside xrg/yrg
    weg=(x>max(xrg) | x<min(xrg) | y>max(yrg) | y<min(yrg));
    x(weg)=[];
    y(weg)=[];
    z(weg)=[];

    % finding the x and y bins for each coordinate
    [~,whichbinx] = histc(x,xrg);   % returns the bin number in x direction (left edge is included in bin)
    [~,whichbiny] = histc(y,yrg);   % returns the bin number in y direction (left edge is included in bin)
 
    % corrected bin numbers (if there is a coordinate == last right bin edge it is set into last bin)
    binsx = min(max(whichbinx,1),length(xrg)-1);    
    binsy = min(max(whichbiny,1),length(yrg)-1);
 
    bins = (binsy-1).*(length(xrg)-1)+binsx; % linear index of bin in matrix
    [test,a,b]=unique(bins);
    bins=bins(a); % only the first occuring trace per bin
    z=z(a); % just takes the first occuring trace per bin
 
    xpos = ones(size(bins,1),1);
    ns = sparse(bins,xpos,1,(length(xrg)-1)*(length(yrg)-1),1); % sparse matrix with number of counts in each bin
    ysum = sparse(bins,xpos,z,(length(xrg)-1)*(length(yrg)-1),1);   % sparse matrix with sum of z values in each bin
    zm = full(ysum)./(full(ns));    % vector with mean values (sum/number) for each bin in vector bins (if no hit, then NaN)
    
    zm = reshape(zm,length(xrg)-1,length(yrg)-1)';  % output matrix
    [xgrid,ygrid]=meshgrid([xrg(1)+dx/2:dx:xrg(end)-dx/2],[yrg(1)+dy/2:dy:yrg(end)-dy/2]);

    ind = ~isnan(zm); % indices of bin values
    vals = zm(ind); % bin values vector
end
