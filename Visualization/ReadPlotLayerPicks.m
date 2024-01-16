clear all
close all
clc

% Reads the picks from LayerPicking and plots them. Please feel free to
% change according to your needs.

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!

% get file name
if ispc
    if exist('picks.temp') % read last opened folder from temp.temp
        fid=fopen('picks.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
        fileattrib('picks.temp','-h');
        fid=fopen('picks.temp','wt');
        fprintf(fid,'%s',foldername);
        fclose(fid);
        fileattrib('picks.temp','+h');
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks'); % path to radargram-folder

        fid=fopen('picks.temp','wt');
        fprintf(fid,'%s',fullfile(folder,file));
        fclose(fid);
        fileattrib('picks.temp','+h');
    end
else
    if exist('.picks.temp') % read last opened folder from temp.temp
        fid=fopen('.picks.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks'); % path to radargram-folder
    end

    fid=fopen('.picks.temp','wt');
    fprintf(fid,'%s',fullfile(folder,file));
    fclose(fid);
end


if ispc
    if exist('temp.temp') % read last opened folder from temp.temp
        fid=fopen('temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
        fileattrib('temp.temp','-h');
        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder

        fid=fopen('temp.temp','wt');
        fprintf(fid,'%s',pfad);
        fclose(fid);
        fileattrib('temp.temp','+h');
    end
else
    if exist('.temp.temp') % read last opened folder from temp.temp
        fid=fopen('.temp.temp','r');
        fn=textscan(fid,'%s');
        fclose(fid);
        if ~isempty(fn{1})
            pfad=uigetdir(fn{1}{1},'Choose folder with radargrams');
        else
            pfad=uigetdir([],'Choose folder with radargrams');
        end
    else
        pfad=uigetdir([],'Choose folder with radargrams'); % path to radargram-folder
    end

    fid=fopen('.temp.temp','wt');
    fprintf(fid,'%s',pfad);
    fclose(fid);
end

%---------------------------------------------------------
% read pick file
fid=fopen(fullfile(folder,file),'r');
anz=fscanf(fid,'%d',1);
for i=1:anz
    temp{i}=textscan(fid,'%s',3);
end
tempcol=textscan(fid,'%f%f%f',anz);
col=[tempcol{1} tempcol{2} tempcol{3}];
% make layer strings 
for i=1:anz
    layerlist{i}=temp{i}{1}{3};
    layerID(i)=str2num(temp{i}{1}{1});
end
temp=textscan(fid,'%f%f%f%f%d%d%d','Headerlines',2);
% picks=[x z E N ID profnum linenum];
picks=[temp{1} temp{2} temp{3} temp{4} double(temp{5}) double(temp{6}) double(temp{7})];
fclose(fid);

%%%
% Now you have picks=[x z E N ID profnum linenum], where x and z (or t) are
% local profile coordinates and E/N global_coordinates. col is a matrix of
% colors for the different layers.
%%%

clear temp, clear tempcol;

 
% Plot maps for every ID
ids=unique(picks(:,5));
for i=1:length(ids)  % for every ID
    figure
    scatter(picks(picks(:,5)==ids(i),3),picks(picks(:,5)==ids(i),4),20,picks(picks(:,5)==ids(i),2),'fill')
    set(gca,'DataAspectratio',[1 1 1])
    colormap(jet)
    title(['LayerID ',int2str(ids(i))])
    grid on
    colorbar
end

%%% 
% If you find out that several IDs are in reality the same layer, you can
% merge them:
% join several IDs to one layer
% picks2=picks;
% picks2(picks2(:,5)==13,7)=picks2(picks2(:,5)==13,7)+10; % increase line numbers
% picks2(picks2(:,5)==13,5)=4; % ID 13 = ID 4
% picks2(picks2(:,5)==2,7)=picks2(picks2(:,5)==2,7)+10;
% picks2(picks2(:,5)==2,5)=7; % 2=7
% picks2(picks2(:,5)==8,7)=picks2(picks2(:,5)==8,7)+10;
% picks2(picks2(:,5)==8,5)=9; % 8=9
% picks2(picks2(:,5)==11,7)=picks2(picks2(:,5)==11,7)+10;
% picks2(picks2(:,5)==11,5)=12; % 11=12
% 
% % Plot maps for every ID with joined layers
% ids=unique(picks2(:,5));
% for i=1:length(ids)  % for every ID
%     figure
%     scatter(picks2(picks2(:,5)==ids(i),3),picks2(picks2(:,5)==ids(i),4),20,picks2(picks2(:,5)==ids(i),2),'fill')
%     set(gca,'DataAspectratio',[1 1 1],'Clim',[-2 0])
%     colormap(jet)
%     title(['LayerID ',int2str(ids(i))])
%     grid on
% end
% 
% picks=picks2; % overwrite picks with joined layers



%%%
% Plot radargrams and picks
rad=[1 12 14 22 30 33 42 52 64]; % only choose some radargrams: number of radargrams
load(fullfile(pfad,'radargrams.mat'));
load(fullfile(pfad,'global_coords.mat'));
load(fullfile(pfad,'t.mat')); % zmig in m!

for i=1:length(rad)
    % calculate profile coordinates:
    xprof{rad(i)}=sqrt((global_coords{rad(i)}(:,1)-global_coords{rad(i)}(1,1)).^2+(global_coords{rad(i)}(:,2)-global_coords{rad(i)}(1,2)).^2);
    
    p=unique(picks(picks(:,6)==rad(i),5)); % IDs in this profile
    pi=picks(picks(:,6)==rad(i),:); % all picks in this profile
    
    % colorlims
    coldata=sort(unique(radargrams{rad(i)}(~isnan(radargrams{rad(i)}(:)))));
    cmin=coldata(round(length(coldata)/100*2));  % 2% colorscale
    cmax=coldata(end-round(length(coldata)/100*2));
    
    figure
    imagesc(xprof{rad(i)},t,radargrams{rad(i)})
    colormap(flipud(gray))
    hold on
    for j=1:length(p)   % for all IDs
        lines=unique(pi(pi(:,5)==p(j),7));  % all line numbers for this ID
        for k=1:length(lines)
            plot(pi(pi(:,5)==p(j) & pi(:,7)==lines(k),1),pi(pi(:,5)==p(j) & pi(:,7)==lines(k),2),'linewidth',2,'Color',col(p(j),:))
        end
    end
    axis xy
    set(gca,'CLim',[cmin cmax],'Fontsize',20)
    title(['Profile ',num2str(rad(i))])
    set(gca,'Dataaspectratio',[4 1 1])
end


%%% 
% Plot picks in bird view:
dz=0.1; % depth interval for cutting
start=-2:dz:0; % depths for cutting of layers
dx=0.2;
[xgrid,ygrid]=meshgrid([min(picks(:,3)):dx:max(picks(:,3))],[min(picks(:,4)):dx:max(picks(:,4))]);
% picks vorbereiten
for j=1:length(ids)
    test=[picks(picks(:,5)==ids(j),3),picks(picks(:,5)==ids(j),4),picks(picks(:,5)==ids(j),2)];
    F=scatteredInterpolant(test(~isnan(test(:,1)),1),test(~isnan(test(:,1)),2),test(~isnan(test(:,1)),3),'linear','none');
    line{j}=F(xgrid,ygrid); % interpolate every layer
end
% plotten
for i=1:length(start)-1
    
    figure
    set(gca,'Dataaspectratio',[1 1 1])
    title(['z = ',num2str(start(i)),' - ',num2str(start(i+1)),' m'])
    axis xy
    hold on
    
    % picks plotten
    for j=1:length(ids)
        contour(xgrid,ygrid,line{j},[start(i)+dz/2 start(i)+dz/2],'Color',col(ids(j),:),'Linewidth',2);
    end
end

