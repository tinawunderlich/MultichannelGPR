function [data_out]=semblanceSmoothing(data,dt,sembwin,dtmin,dtmax,dtinc,sembexp,sembflag,qbal,qclip)

%%% Semblance smoothing of 3d-data-cube
% [data_out]=semblanceSmoothing(data,dt,sembwin,dtmin,dtmax,dtinc,sembexp,sembflag,qbal,qclip)
%
% code based on Wilken et al. 2019: Imaging a medieval shipwreck with the
% new PingPong 3D marine reflection seismic system. Archaeological
% Prospection.
% Originally written by Dr. Dennis Wilken, modified by Dr. Tina
% Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3d-datacube wiht traces along third dimension
% dt: time sampling interval
% sembwin: window for semblance in bins, has to be odd!
% dtmin: minimum time difference in sembwin [ns]
% dtmax: maximum time difference in sembwin [ns]
% dtinc: increment between dtmin and dtmax [ns]
% sembexp: filter exponent to calculate coherence from semblance^sembexp
% sembflag: =1: use double cones for semblance calculation, =2: use round rotating slice of diameter sembwin (is rotated in steps of 10°)
% qbal: normalize datacube after filtering? (yes=1)
% qclip: probability for normalization using quantiles in interval [0,1] (e.g. 0.98)
% 
% Output: data_out: semblance smoothed 3d-datacube (same size as data)
%


% get matlab version (to shut down parfor_wait)
rel=version('-release');
if str2num(rel(1:4))<2017
    pfw=0; % no parfor_wait
else
    pfw=1; % use parfor_wait
end


data(isnan(data))=0;    % set Nans in data to zero

if round(sembwin/2)==sembwin/2   % test if sembwin is odd
    sembwin=sembwin+1;
    disp(['sembwin is even. Using sembwin=',int2str(sembwin),' instead.']);
end

% create semblance operator in apex point (0,0,0) for all slopes
if sembflag==1 % using double cones
    diffdt=[dtmin:dtinc:dtmax]./dt;   % different dts in sembwin in time-bins
    h=round(diffdt./2);   % height of conus in time-bins
    R=((sembwin-1)/2)./h;    % radius in conus height 1
    xd=cell(1,length(diffdt));
    yd=cell(1,length(diffdt));
    td=cell(1,length(diffdt));
    for i=1:length(diffdt)  % for every slope

        % possible points in sembwin (bins):
        [ix,iy,it]=meshgrid([-(sembwin-1)/2:(sembwin-1)/2],[-(sembwin-1)/2:(sembwin-1)/2],[-h(i):h(i)]);

        % find points approx. on surface of conus:
        onsurface=(fix(ix.^2+iy.^2)==fix(R(i)^2*it.^2));

        xd{i}=[0 ix(onsurface)'];  % 0 is for apex point
        yd{i}=[0 iy(onsurface)'];
        td{i}=[0 it(onsurface)'];
    end
elseif sembflag==2 % use rotating round slice
    diffdt=[dtmin:dtinc:dtmax]./dt;   % different dts in sembwin in time-bins
    rota=0:10:170;  % rotation angle of slice in °
    xd=cell(1,length(diffdt)*length(rota));
    yd=cell(1,length(diffdt)*length(rota));
    td=cell(1,length(diffdt)*length(rota));
    
    % possible points in sembwin (bins):
    [ix,iy]=meshgrid([-(sembwin-1)/2:(sembwin-1)/2],[-(sembwin-1)/2:(sembwin-1)/2]);
    % find points in circle
    incircle=(sqrt(ix.^2+iy.^2)<=(sembwin-1)/2);

    for i=1:length(diffdt)   % for all slopes
        for j=1:length(rota)
            ij=(i-1)*length(rota)+j;
            % find x and y bins in circle
            xd{ij}=[ix(incircle)];
            yd{ij}=[iy(incircle)];
            
            iz=repmat(linspace(-round(diffdt(i)/2),round(diffdt(i)/2),sembwin),[sembwin,1]);    % t-values corresponding to ix and iy
            ztemp=iz(incircle); % t values (bins) for xd/yd
            Ra=[cosd(rota(j)) -sind(rota(j)) 0; sind(rota(j)) cosd(rota(j)) 0; 0 0 1];  % rotational matrix
            for k=1:length(xd{ij})    % rotate every point
                xytemp=[xd{ij}(k); yd{ij}(k); ztemp(k)];
                temp(:,k)=Ra*xytemp;    % one column for every point (1. line is x, 2. line is y and 3. line is t)
            end

            F=scatteredInterpolant(temp(1,:)',temp(2,:)',temp(3,:)');
            td{ij}=round(F(xd{ij},yd{ij}));
        end
    end
end



Nx=length(data(1,:,1));
Ny=length(data(:,1,1));
ns=length(data(1,1,:));
[ix,iy,it]=meshgrid([1:Nx],[1:Ny],[1:ns]);  % grids with indices

% initialize arrays:  
S=cell(1,numel(data));  % semblance
maxS=zeros(size(data)); % maximum semblance
whichS=zeros(size(data));   % which slope has maximum semblance for each point?
datanew=zeros(size(data));  % new data cube

n=numel(data);
if pfw==1
    WaitMessage = parfor_wait(n, 'Waitbar', false,'ReportInterval',100);
end
parfor i=1:n  % for every sample and every trace

    for j=1:length(xd)  % for every slope/rotation
        % find points in range of dataset
        inn=(iy(i)+yd{j}>=1 & iy(i)+yd{j}<=Ny & ix(i)+xd{j}>=1 & ix(i)+xd{j}<=Nx & it(i)+td{j}>=1 & it(i)+td{j}<=ns);
        
        top=sum(sum(sum(data(iy(i)+yd{j}(inn),ix(i)+xd{j}(inn),it(i)+td{j}(inn)))))^2;
        bottom=sum(sum(sum(data(iy(i)+yd{j}(inn),ix(i)+xd{j}(inn),it(i)+td{j}(inn)).^2)));
        
        % Semblance for current point:
        S{i}(j)=top./bottom; % (Wilken et al 2019, eq. 2)
        if isnan(S{i}(j))
            S{i}(j)=0;
        end
    end
    % maximum semblance
    [maxS(i),whichS(i)]=max(S{i}); % (eq. 3)
 
    % find points on surface of conus in range of dataset for maximum
    % semblance
    inn=(iy(i)+yd{whichS(i)}>=1 & iy(i)+yd{whichS(i)}<=Ny & ix(i)+xd{whichS(i)}>=1 & ix(i)+xd{whichS(i)}<=Nx & it(i)+td{whichS(i)}>=1 & it(i)+td{whichS(i)}<=ns);
      
    % smooth data (eq. 4)
    temp=data(iy(i)+yd{whichS(i)}(inn),ix(i)+xd{whichS(i)}(inn),it(i)+td{whichS(i)}(inn)); % Amplitude at points
    datanew(i)=mean(temp(:));
      
    if pfw==1
        WaitMessage.Send;
    end
end
if pfw==1
    WaitMessage.Destroy;
end

% coherence filter (eq. 5)
data_out=datanew.*maxS.^sembexp;


% NORMALIZE TRACES:
if(qbal==1)
    data_out=normalize3d(data_out,qclip);
end
