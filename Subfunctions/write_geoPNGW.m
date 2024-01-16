function []=write_geoPNGW(x,y,coordtrans,filename)

%%% []=write_geoPNGW(x,y,coordtrans,filename)
% function for writing of pgw-file for georeferenced ongs
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% x,y: meshs with local coordinates (created by meshgrid)
% coordtrans: matrix with local coordinates in columns 1:2 and
% corresponding global coordinates in columns 3:4
% filename: path and filename for PGW-file
%
% PNG has to be written separately!


dx=x(1,2)-x(1,1);
dy=y(2,1)-y(1,1);

% determine global coords of upper left and upper right pixel
pix_ol=helmert([min(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
pix_or=helmert([max(x(:)) max(y(:))],coordtrans(:,1:2),coordtrans(:,3:4));
alpha=atand(abs(pix_ol(2)-pix_or(2))/abs(pix_ol(1)-pix_or(1))); % angle against west
if pix_ol(2)>pix_or(2)
    alpha=-alpha;
end
% determine pixel-lengths in all directions
A=dx*cosd(alpha);
D=dx*sind(alpha);
E=dy*cosd(alpha);
B=dy*sind(alpha);


    
% write pngw
fid=fopen(filename,'wt');
fprintf(fid,[num2str(A),'\n',num2str(D),'\n',num2str(B),'\n',num2str(-E),'\n',num2str(pix_ol(1)),'\n',num2str(pix_ol(2))]);
fclose(fid);
