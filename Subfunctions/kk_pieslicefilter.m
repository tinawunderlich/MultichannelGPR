function [data_new]=kk_pieslicefilter(data,xgrid,ygrid,dx,angle,opening,plotflag,colorclip,sname)

%% function for kk_pieslicefilter to reduce stripe noise
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: gridded data with corresponding xgrid/ygrid
% dx: sample spacing in same unit as x/y
% angle: middle angle of pieslice in degree (0° is east-west, -45° is northeast/southwest, 45° is northwest/southeast, 90° is north/south)
% opening: opening angle of pieslice in degree
% plotflag: =1 -> plot, =0 -> no not plot
% colorclip: percent value for colour clipping
% sname: path+name of figure for spectrum saving
% 
% Output:
% data_new: gridded filtered data, same size as data


data(isnan(data))=0;


%%% 2D-FFT
Fx=1/dx;    % Samplingfrequency in x-direction
Fy=1/dx;    % Samplingfrequency in y-direction
Lx=length(data(1,:));   % number of samples in x-direction
Ly=length(data(:,1));   % number os samples in y-direction
nfftx=2^nextpow2(Lx);   
nffty=2^nextpow2(Ly);   

data_f=fft2(data,nffty,nfftx);  % data in kk-domain
data_f=fftshift(data_f);    % Shift quadrants, so that k=0 is in the middle

ky=Fy/2*linspace(0,1,nffty/2+1); % wavenumber vector y (onesided)
kx=Fx/2*linspace(0,1,nfftx/2+1); % wavenumber vector x (onesided)
kx2=[-fliplr(kx) kx(2:end-1)];    % wavenumber vector x (twosided)
ky2=[-fliplr(ky) ky(2:end-1)];    % wavenumber vector y (twosided)
[kx2g,ky2g]=meshgrid(kx2,ky2);      % wavenumber grids


%%% plot spectrum
if plotflag==1
    fh_spectrum=figure('Position',[0 0 1400 700]);
    subplot(1,2,1)
    imagesc(ky2,-kx2,abs(data_f'))
    % determine color limits
    coldata=sort(unique(abs(data_f')));
    coldata(isnan(coldata))=[]; % delete nans
    cmaxs=coldata(end-round(length(coldata)/100*1));
    set(gca,'FontSize',20,'Clim',[0 cmaxs])
    title('Spectrum before filtering')
    xlabel('ky [1/m]')
    ylabel('kx [1/m]')
    axis xy
    axis equal
    axis tight
    set(gca,'Dataaspectratio',[1 1 1])
    cb=colorbar;
    set(cb,'fontSize',20)
end



%%% create filter mask in kk-domain with cosine tapering
filtermask=ones(size(data_f));
mina=angle-opening/2;
maxa=angle+opening/2;
if mina<90 && maxa>90
    filtermask(mina<=atan2d(kx2g,ky2g) & maxa>=atan2d(kx2g,ky2g))=0;
    filtermask=[fliplr(flipud(filtermask(:,nfftx/2+1:end))) filtermask(:,nfftx/2+1:end)];
    cosmask=0.5+0.5*cosd(360/(maxa-mina)*(atan2d(kx2g,ky2g)-mina));
    cosmask=[fliplr(flipud( cosmask(:,nfftx/2+1:end)))  cosmask(:,nfftx/2+1:end)]; %neu
elseif mina<-90 && maxa>-90
    filtermask(mina<=atan2d(kx2g,ky2g) & maxa>=atan2d(kx2g,ky2g))=0;
    filtermask=[filtermask(:,1:nfftx/2) fliplr(flipud(filtermask(:,1:nfftx/2)))];
    cosmask=0.5+0.5*cosd(360/(maxa-mina)*(atan2d(kx2g,ky2g)-mina));
    cosmask=[cosmask(:,1:nfftx/2) fliplr(flipud(cosmask(:,1:nfftx/2)))]; %neu
else
    filtermask(mina<=atand(kx2g./ky2g) & maxa>=atand(kx2g./ky2g))=0;
    cosmask=0.5+0.5*cosd(360/(maxa-mina)*(atand(kx2g./ky2g)-mina));
end 
filtermask(filtermask==0)=cosmask(filtermask==0);

% apply filtermask
data_f_out=data_f.*filtermask;

% plot filtered spectrum
if plotflag==1
    gcf=fh_spectrum;
    subplot(1,2,2)
    imagesc(ky2,-kx2,abs(data_f_out'))
    set(gca,'FontSize',20,'Clim',[0 cmaxs])
    title('Spectrum after filtering')
    xlabel('ky [1/m]')
    ylabel('kx [1/m]')
    axis xy
    axis equal
    axis tight
    cb=colorbar;
    set(cb,'fontSize',20)

    sgtitle(['angle = ',num2str(angle),'° / opening angle = ',num2str(opening),'°'])

    % save figure
    saveas(fh_spectrum,sname,'png')
    close(fh_spectrum);
end


%%% 2D-IFFT
data_f_out=fftshift(data_f_out);  % shift quadrants back

data_new=ifft2(data_f_out,nffty,nfftx,'symmetric'); % Inverse 2D-fft
data_new=data_new(1:length(data(:,1)),1:length(data(1,:)));  % cut to original size
