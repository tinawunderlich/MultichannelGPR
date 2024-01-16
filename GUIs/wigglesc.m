function wigglesc(data,t,x,scal)

%% scaled wiggle plot function
% Dr. Tina Wunderlich, CAU Kiel, July 2022
% tina.wunderlich@ifg.uni-kiel.de
%
% data: radargram with traces along columns
% x: vector of profile coordinates with same length as data(1,:)
% t: time vector of length(data(:,1))
% scal: scaling factor, if=1, max amplitude is scaled as dx

zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % anonymous function for finding the indices of zero crossings

dx=mean(diff(x));
if length(t(:,1))>length(t(1,:))
    t=t';
end

hold on
for i=1:length(data(1,:)) % for every trace
    data(:,i)=data(:,i)-mean(data(:,i),'omitnan'); % subtract mean amplitude
    V=[data(:,i)/max(data(:,i))*dx*scal t'];
    V(isnan(V(:,1)),:)=[]; % delete rows with nan
    ind=zci(V); % indices of zero crossings
    V(ind,1)=0; % set zero crossings to zero
    F1=1:length(V(:,1)); % white polygon = negative amplitudes
    F1(V(:,1)>0)=[];
    F2=1:length(V(:,1)); % black polygon = positive amplitudes
    F2(V(:,1)<0)=[];
    F=NaN(max([length(F1),length(F2)]),2);
    F(1:length(F1),1)=F1;
    F(1:length(F2),2)=F2;
    V(:,1)=V(:,1)+x(i);
    patch('Faces',F(:,1)','Vertices',V,'FaceColor','none');
    patch('Faces',F(:,2)','Vertices',V,'FaceColor','k');
end
