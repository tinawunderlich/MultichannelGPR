function []=export2sgy2D(traces,dt,x,y,pathfilename,z,constoff)

% []=export2sgy2D(traces,dt,x,y,pathfilename,z)
% Save sgy-file for radargram in matrix traces and info txt-file with settings
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: matrix with traces of one radargram
% dt: sampling interval in ns
% x,y: vectors with coordinates in m (same length as number of traces)
% pathfilename: complete path and filename
% z: height (optional, otherwise set to zero)
% constoff: if==1: constant offset of coordinates will be subtracted, and will be saved in Inline3D and
% Crossline3D (optional, if =1: mm accuracy of coordinates will be used)
%
% Mean trace spacing is used as trace spacing!

if nargin==5
    z=zeros(size(x));
    constoff=0;
elseif nargin==6
    constoff=0;
end

% check for nans in data and replace with zero
traces(isnan(traces))=0;

% force x and y and z to be row vectors
if length(x(:,1))>length(x(1,:))
    x=x';
end
if length(y(:,1))>length(y(1,:))
    y=y';
end
if length(z(:,1))>length(z(1,:))
    z=z';
end

% determine mean trace spacing
dist=diff(sqrt((x(1)-x).^2+(y(1)-y).^2));
meandist=mean(dist);
len=sum(dist);
xprof=0:meandist:len;
if length(xprof)<length(traces(1,:))
    xprof=zeros(1,length(traces(1,:))); % set all coordinates to zero
end

% Write sgy
if constoff==0 % no constant coordinate offset, coordinates in cm accuracy
    writesegy(pathfilename,traces,'dt',dt,'cdpX',xprof,'cdpY',zeros(size(xprof)),'SourceX',x,'SourceY',y,'GroupX',x,'GroupY',y,'gelev',z,'selev',z);
else
    % determine constant coordinate offset:
    offX=round(x(1)/1000)*1000;
    offY=round(y(1)/1000)*1000;
    % subtract offset from x and y and multiply with 10:
    x=(x-offX).*10;
    y=(y-offY).*10;
    writesegy(pathfilename,traces,'dt',dt,'cdpX',xprof,'cdpY',zeros(size(xprof)),'SourceX',x,'SourceY',y,'GroupX',x,'GroupY',y,'gelev',z,'selev',z,'Inline3D',offX,'Crossline3D',offY);
end

% Write info-file
[pathstr,fname,ext]=fileparts(pathfilename);
fid=fopen(fullfile(pathstr,[fname,'.txt']),'wt');
fprintf(fid,['dt ',num2str(dt,'%8.4f'),' ns\n']);
fprintf(fid,['ns ',int2str(length(traces(:,1))),'\n']);
fprintf(fid,['ntraces ',int2str(length(x)),'\n']);
if constoff==0
    fprintf(fid,['x_start ',num2str(x(1),'%8.2f'),' m\n']);
    fprintf(fid,['x_end ',num2str(x(end),'%8.2f'),' m\n']);
    fprintf(fid,['y_start ',num2str(y(1),'%8.2f'),' m\n']);
    fprintf(fid,['y_end ',num2str(y(end),'%8.2f'),' m\n']);
else
    fprintf(fid,['x_start ',num2str(x(1)./10+offX,'%8.2f'),' m\n']);
    fprintf(fid,['x_end ',num2str(x(end)./10+offX,'%8.2f'),' m\n']);
    fprintf(fid,['y_start ',num2str(y(1)./10+offY,'%8.2f'),' m\n']);
    fprintf(fid,['y_end ',num2str(y(end)./10+offY,'%8.2f'),' m\n']);
end
fprintf(fid,['mean_trace_spacing ',num2str(meandist,'%8.4f'),' m\n']);
fprintf(fid,['length_of_profile ',num2str(len,'%8.2f'),' m\n']);
fprintf(fid,['profile coordinates in cdpX\n']);
if constoff==0
    fprintf(fid,['x/y coordinates in SourceX/Y & GroupXY in cm\n']);
else
    fprintf(fid,['reduced x/y coordinates in SourceX/Y & GroupXY in mm\n']);
    fprintf(fid,['Constant offsets of x/y coordinates in Inline3D (x) and Crossline3D (y) in m\n']);
end
fprintf(fid,['z coordinate in gelev & selev\n']);
fclose(fid);
