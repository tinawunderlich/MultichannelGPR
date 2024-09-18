function wigglesc(data,t,x,scal,ax)

%% scaled wiggle plot function
% Dr. Tina Wunderlich, CAU Kiel, July 2022
% tina.wunderlich@ifg.uni-kiel.de
%
% data: radargram with traces along columns
% x: vector of profile coordinates with same length as data(1,:)
% t: time vector of length(data(:,1))
% scal: scaling factor, if=1, max amplitude is scaled as dx
% ax: axes to plot in (optional)

if nargin==4
    ax=gca;
end

zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % anonymous function for finding the indices of zero crossings

dx=mean(diff(x));
if length(t(:,1))>length(t(1,:))
    t=t';
end

data=data-mean(data,1,'omitnan'); % subtract mean amplitude for each trace
data(1,:)=0;
data(end,:)=0;
n=length(data(:,1)); % number of samples

for i=1:length(data(1,:)) % for every trace
    % Vertices:
    V=[data(:,i)/max(data(:,i))*dx*scal t']; % scaled amplitude and time for this trace
    V=V(~isnan(V(:,1)),:); % delete rows with nan 
    if size(V,1)<n
        V=[V; 0 max(t)]; % add last point if it is missing
    end
    ind=zci(V); % indices of zero crossings
    ind=ind(ind<=n); % ind has to be <= number of samples
    V(ind,1)=0; % set zero crossing to 0
    V(end,1)=0; % set last point to zero
    % Faces:
    F1=1:length(V(:,1));
    F1=F1(V(:,1)<=0);  % white polygon = negative amplitudes (indices of points with negative amplitudes)
    F1=[1 F1 n]; % add first and last sample
    F2=1:length(V(:,1)); 
    F2=F2(V(:,1)>=0); % black polygon = positive amplitudes
    F2=[1 F2 n]; % add first and last sample
    V(:,1)=V(:,1)+x(i); % add x-value of trace to amplitudes
    patch(ax,'Faces',F1,'Vertices',V,'FaceColor','none');
    patch(ax,'Faces',F2,'Vertices',V,'FaceColor','k');
end
