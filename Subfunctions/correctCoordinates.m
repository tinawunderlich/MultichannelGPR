function [data_out,trh_out,deltraces]=correctCoordinates(data,trh,h_GNSS,offsetGNSS_X,offsetGNSS_Y,removeStartEnd,smooth_coords,numsamp_smooth)

% Function for correcting measured GNSS coordinates with offsets between
% GNSS antenna and GPR antenna and GPS height also taking the topography
% into account
%
% Dr. Tina Wunderlich, CAU Kiel, 2020-2026
%
% Input:
% data: traces as matrix (or if no traces shall be removed =[])
% trh: trace header with fields x, y, z and optionally other fields (if
% channum: correction is done for each channel individually)
% h_GNSS: heigth of GNSS antenna above ground [m]
% offsetGNSS_X,offsetGNSS_Y: offsets between GNSS antenna and GPR antenna
% inline (X) and crossline (Y) [m]
% (in profile direction GNSS left of antenna -> positive X)
% (if GNSS behind antenna midpoint -> positive Y)
% removeStartEnd: if =1: remove standstill traces at beginning and end of profile
% smooth_coords: if =1: smooth coordinates with movmean using
% numsamp_smooth samples
% numsamp_smooth: number of samples for smoothing
%
%
% Output:
% data: corrected traces as matrix (opt. standstill traces removed)
% trh: corrected coordinates in trace header fields x, y, z
% deltraces: logical vector containing those traces of original data that need to be deleted

if nargin==5
    removeStartEnd=0;
    smooth_coords=0;
    numsamp_smooth=1;
end

if isfield(trh,'channum')
    channels=unique(trh.channum);
else
    channels=1;
    trh.channum=ones(size(trh.x));
end

% get all field names in struct trh
fn=fieldnames(trh);

savetrh=trh; % save old trh for splitting
savedata=data; % save original data for splitting

% initialize new output:
trh_out=trh;
for ii=1:length(fn) % for all fields in struct take only one channel
    trh_out.(fn{ii})=[]; % set all fields empty
end
data_out=[];

for ch=channels
    %%% SPLIT in different channels first
    for ii=1:length(fn) % for all fields in struct take only one channel
        trh.(fn{ii})=savetrh.(fn{ii})(ch==savetrh.channum);
    end
    data=savedata(:,ch==savetrh.channum);

    % if there are traces at the same position -> delete
    if ~isempty(data)
        d1=diff(trh.x);
        d2=diff(trh.y);
        d1=[d1(:); d1(end)];
        d2=[d2(:); d2(end)];
        deltraces=(d1==0 & d2==0);
        if removeStartEnd==1
            diffe = sqrt(diff(trh.x).^2 + diff(trh.y).^2);
            standstill = diffe <= max(diffe(1),diffe(end));
            % check how many points would be affected:
            if sum(standstill)<length(standstill)/6 % only delete if <1/6 of profile length
                first = find(~standstill,1,'first')-1;
                last = length(trh.x) - find(fliplr(~standstill),1,'first')-1;
                deltraces(1:first) = 1;
                deltraces(last:end) = 1;
            end
        end

        % correct all fields in trh:
        for ii=1:length(fn)
            trh.(fn{ii})=trh.(fn{ii})(~deltraces);
        end
        % also delete traces in data:
        data=data(:,~deltraces);
    end

    if smooth_coords==1
        trh.x=movmean(trh.x,numsamp_smooth);
        trh.y=movmean(trh.y,numsamp_smooth);
        trh.z=movmean(trh.z,numsamp_smooth);
    end

    % initialize temporary coordinates:
    xn=zeros(size(trh.x));
    yn=zeros(size(trh.x));
    zn=zeros(size(trh.x));

    % now correct for offsets:
    anz2=round(0.5/mean(sqrt(diff(trh.x).^2+diff(trh.y).^2))); % number of points for direction determination (using mean trace spacing for 0.5 m distance)
    if anz2/2==round(anz2/2)
        anz2=anz2+1; % make odd
    end
    for ii=1:length(trh.x)-anz2
        dist=sqrt((trh.x(ii)-trh.x(ii+anz2))^2+(trh.y(ii)-trh.y(ii+anz2))^2); % distance between two points anz2-traces away
        dist_xy=[trh.x(ii+anz2)-trh.x(ii) trh.y(ii+anz2)-trh.y(ii)]; % differences in x and y direction
        dist_all=sqrt((trh.x(ii)-trh.x(ii:ii+anz2)).^2+(trh.y(ii)-trh.y(ii:ii+anz2)).^2); % distance to following points until anz2-traces away
        % correct for GPS height and associated xy-shift if
        % topography is present
        if any(trh.z~=0)
            % fit line to height of following points:
            p=polyfit(dist_all,trh.z(ii:ii+anz2),1);
            alpha=atand(abs(polyval(p,0)-polyval(p,dist))/dist); % mean slope angle in degree
            hnew=h_GNSS/cosd(alpha); % new height of GPS above ground at measured point
            b=hnew*sind(alpha); % shift along surface to correct point
            d_x=b*cosd(alpha); % horizontal shift of point in profile direction
            d_z=b*sind(alpha); % vertical shift of point in profile direction

            % angle of profile direction towards north
            if trh.x(ii)==trh.x(ii+anz2) % north/south
                beta=0;
            elseif trh.y(ii)==trh.y(ii+anz2) % west/east
                beta=90;
            else
                beta=atand(abs(trh.x(ii+anz2)-trh.x(ii))/abs(trh.y(ii+anz2)-trh.y(ii)));
            end
            % split horizontal shift into x and y component:
            d=[d_x*cosd(beta) d_x*sind(beta)];
            % differentiate between different cases
            updownflag=trh.z(ii)<trh.z(ii+anz2); % up=1, down=0
            if all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards NNE
                d_xx=min(d);
                d_xy=max(d);
            elseif all(dist_xy>0) && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards NNE
                d_xx=-min(d);
                d_xy=-max(d);
            elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards NEE
                d_xx=max(d);
                d_xy=min(d);
            elseif all(dist_xy>0) && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards NEE
                d_xx=-max(d);
                d_xy=-min(d);
            elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards SEE
                d_xx=max(d);
                d_xy=-min(d);
            elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards SEE
                d_xx=-max(d);
                d_xy=min(d);
            elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards SSE
                d_xx=min(d);
                d_xy=-max(d);
            elseif dist_xy(1)>0 && dist_xy(2)<0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards SSE
                d_xx=-min(d);
                d_xy=max(d);
            elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==1 % up, towards SSW
                d_xx=-min(d);
                d_xy=-max(d);
            elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))<abs(dist_xy(2)) && updownflag==0 % down, towards SSW
                d_xx=min(d);
                d_xy=max(d);
            elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==1 % up, towards SWW
                d_xx=-max(d);
                d_xy=-min(d);
            elseif dist_xy(1)<0 && dist_xy(2)<0 && abs(dist_xy(1))>=abs(dist_xy(2)) && updownflag==0 % down, towards SWW
                d_xx=max(d);
                d_xy=min(d);
            elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==1 % up, towards NWW
                d_xx=-max(d);
                d_xy=min(d);
            elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))>abs(dist_xy(2)) && updownflag==0 % down, towards NWW
                d_xx=max(d);
                d_xy=-min(d);
            elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==1 % up, towards NNW
                d_xx=-min(d);
                d_xy=max(d);
            elseif dist_xy(1)<0 && dist_xy(2)>0 && abs(dist_xy(1))<=abs(dist_xy(2)) && updownflag==0 % down, towards NNW
                d_xx=min(d);
                d_xy=-max(d);
            elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==1 % up, towards N
                d_xx=0;
                d_xy=max(d);
            elseif dist_xy(1)==0 && dist_xy(2)>0 && updownflag==0 % down, towards N
                d_xx=0;
                d_xy=-max(d);
            elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==1 % up, towards S
                d_xx=0;
                d_xy=-max(d);
            elseif dist_xy(1)==0 && dist_xy(2)<0 && updownflag==0 % down, towards S
                d_xx=0;
                d_xy=max(d);
            elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==1 % up, towards E
                d_xx=max(d);
                d_xy=0;
            elseif dist_xy(2)==0 && dist_xy(1)>0 && updownflag==0 % down, towards E
                d_xx=-max(d);
                d_xy=0;
            elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==1 % up, towards W
                d_xx=-max(d);
                d_xy=0;
            elseif dist_xy(2)==0 && dist_xy(1)<0 && updownflag==0 % down, towards W
                d_xx=max(d);
                d_xy=0;
            else
                d_xx=NaN;
                d_xy=NaN;
            end

            % coordinate transformation
            localC=[zeros(numel(dist_all),1) dist_all(:)];
            globalC=[trh.x(ii:ii+anz2)+d_xx trh.y(ii:ii+anz2)+d_xy];
            temp=helmert([offsetGNSS_X offsetGNSS_Y],localC,globalC);

            xn(ii)=temp(1);
            yn(ii)=temp(2);
            zn(ii)=trh.z(ii)-hnew+d_z; % correct height
        else % no topography present
            localC=[zeros(numel(dist_all),1) dist_all(:)];
            globalC=[trh.x(ii:ii+anz2) trh.y(ii:ii+anz2)];
            temp=helmert([offsetGNSS_X offsetGNSS_Y],localC,globalC);
            xn(ii)=temp(1);
            yn(ii)=temp(2);
        end
    end
    % extrapolate for the last few traces:
    ende=length(trh.x);
    dist=sqrt((trh.x(length(trh.x)-anz2+1)-trh.x(ende))^2+(trh.y(length(trh.x)-anz2+1)-trh.y(ende))^2);
    for ii=length(trh.x)-anz2+1:ende-1
        interdist(ii,1)=sqrt((trh.x(ii)-trh.x(ii+1))^2+(trh.y(ii)-trh.y(ii+1))^2);
    end
    interdist=interdist(length(trh.x)-anz2+1:ende-1);
    transdata=[repmat(offsetGNSS_X,[numel(interdist)+1,1]) repmat(offsetGNSS_Y,[numel(interdist)+1,1])+[0; interdist]];

    % coordinate transformation
    temp=helmert(transdata,[0 0; 0 dist],[trh.x(length(trh.x)-anz2+1)+d_xx trh.y(length(trh.x)-anz2+1)+d_xy; trh.x(ende)+d_xx trh.y(ende)+d_xy]);

    xn(length(trh.x)-anz2+1:ende)=temp(:,1);
    yn(length(trh.x)-anz2+1:ende)=temp(:,2);
    if any(trh.z~=0)
        zn(length(trh.x)-anz2+1:ende)=trh.z(length(trh.x)-anz2+1:ende)-hnew+d_z;
    end



    % Output:
    data_out=[data_out data];

    % set new/corrected coordinates in original fields in trace
    % header
    trh_out.x=[trh_out.x xn];
    trh_out.y=[trh_out.y yn];
    trh_out.z=[trh_out.z zn];
    trh_out.channum=[trh_out.channum zeros(size(zn))+ch];
    if isfield(trh_out,'tracenum')
        trh_out.tracenum=[trh_out.tracenum 1:length(zn)];
    else
        trh_out.tracenum=[1:length(zn)]';
    end
    % also add other fields:
    for i=1:length(fn)
        if ~strcmp(fn{i},'x') && ~strcmp(fn{i},'y') && ~strcmp(fn{i},'z') && ~strcmp(fn{i},'channum') && ~strcmp(fn{i},'tracenum')
            trh_out.(fn{i})=[trh_out.(fn{i}) trh.(fn{i})];
        end
    end

    clear xn yn zn tempx_flip tempy_flip tempz_flip;
end

