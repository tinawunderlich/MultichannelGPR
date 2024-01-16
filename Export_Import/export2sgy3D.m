function []=export2sgy3D(pathfilename,dt,ns,x,y,data,topo,coordtrans)

% []=export2sgy3D(pathfilename,dt,ns,x,y,data,topo,coordtrans)
% Export sgy-data (4-byte IEEE floating point, big endian, Revision 1) of 3D binned
% area and coordinate info file for Kingdom
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
%
% Input:
% pathfilename: Path and filename of sgy-file
% dt: Sampleinterval in ns
% ns: Number of samples
% x, y: x and y grid corresponding to data(:,:,1) in m
% data: 3D binned data matrix with time in third dimension
% topo: gridded topography corresponding to x and y grids in m
% coordtrans: (Optional!) matrix for coordinate transformation, if area has
% been rotated
%
% requires helmert.m (in folder Subfunctions)



if nargin<8
    ct=0;   % no coordinate transformation
elseif nargin==8
    ct=1;   % coordinate transformation has been applied -> use for sgy export
end


% all traces matrix
n=length(data(:,1,1))*length(data(1,:,1));
datamat=zeros(ns,n);
I=zeros(1,n);
J=zeros(1,n);
x3=zeros(1,n);
y3=zeros(1,n);
z3=zeros(1,n);
for i=1:n
    [I(i),J(i)]=ind2sub(size(x),i);  % bin indices
    if ~isnan(data(I(i),J(i),1))
        datamat(:,i)=data(I(i),J(i),:);
    else
        datamat(:,i)=zeros(size(data(I(i),J(i),:)));
    end
    x3(i)=x(I(i),J(i));
    y3(i)=y(I(i),J(i));
    z3(i)=topo(I(i),J(i));
end

% Write sgy
writesegy3D(pathfilename,datamat,'dt',dt,'cdpX',x3,'cdpY',y3, 'Inline',I,'Crossline',J,'topo',z3);

% Write coordinate info file for Kingdom
[pathstr,name,ext] = fileparts(pathfilename);
fid=fopen(fullfile(pathstr,[name,'_CoordinateInfo.dat']),'wt');
if ct==0
    fprintf(fid,'%d\t%d\t%8.2f\t%8.2f\n',[1 1 x(1,1) y(1,1); length(x(:,1)) 1 x(end,1) y(end,1); 1 length(x(1,:)) x(1,end) y(1,end)]'); % inline Crossline x y
elseif ct==1
    coordtemp=[x(1,1) y(1,1); x(end,1) y(end,1); x(1,end) y(1,end)];
    newtemp=helmert(coordtemp,coordtrans(:,1:2),coordtrans(:,3:4));
    fprintf(fid,'%d\t%d\t%8.2f\t%8.2f\n',[1 1 newtemp(1,:); length(x(:,1)) 1 newtemp(2,:); 1 length(x(1,:)) newtemp(3,:)]'); % inline Crossline x y
end
fclose(fid);

% Write general info file for Kingdom
fid=fopen(fullfile(pathstr,[name,'_Info.txt']),'wt');
fprintf(fid,'Inline3D Minimum Maximum: %d\t%d\n',[1 max(I)]);
fprintf(fid,'Crossline3D Minimum Maximum: %d\t%d\n',[1 max(J)]);
fprintf(fid,'Bin-interval x y: %4.2f\t%4.2f\n',[x(1,2)-x(1,1) y(2,1)-y(1,1)]);
fprintf(fid,'dt: %6.4f ns\n',dt);
fprintf(fid,'SGY: Revision 1, 4-byte IEEE floating point, big endian');
fclose(fid);