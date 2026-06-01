function [mig,zaxis] = topoMig3D(data,topo,x,y,dt,dz,v,apertur,nCores)
    % helper function for 3D semi circle migration
    % by Morten Harms 2026
    
    % get path to semiCircle.cc and compile
    tmp = pwd;
    tmp = split(tmp,'/');
    tmp(end) = {'Migration'};
    mexPath = char(join(tmp,'/'));
    compileMEX(mexPath);

    % preparations
    [Nx, Ny, ns] = size(data);
    mask = double(any(data,3));

    % number of bins containing data
    nTraces = nnz(mask);
    

    % prepare topoMig
    % reference height = highest point in survey area
    href = max(topo(:));

    % topoMig is depth below reference height
    topoMig = href - topo;
    topoMig(isnan(topoMig)) = 0;

    % prepare Nz and zaxis
    toporange = max(topo(:)) - min(topo(:));
    NzTopo    = round(toporange/dz);
    tmax     = ns * dt;
    zmaxData = tmax/2 * v;
    NzData   = round(zmaxData/dz);
    Nz       = NzData + NzTopo;
    zaxis = double(0:Nz-1) * dz;

    % prepare 2D axis
    xaxis = x(1,:);
    yaxis = y(:,1);

    % migration
    mig = semiCircle(data, mask, topoMig, Nx, Ny, Nz, ns, nTraces, dz, apertur, xaxis, yaxis, zaxis,v, dt, mask, nCores);
    

    % change zaxis to elevation based on topo
    zaxis = href - zaxis;

end
