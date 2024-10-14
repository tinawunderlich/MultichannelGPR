clear all
close all
clc

% Script for creating a 1D velocity function from picks (*.txt-file made with
% Velocity_2D.m)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de

fit=1;  % Fitting funtion: =1 polynom, =2 exp
ord=1;  % Order of polynom for fitting of 1D velocity function

%--------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!

% get file name
if ispc
    if exist('v1d.temp') % read last opened folder from temp.temp
        fid=fopen('v1d.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
    end
    % save last selected folder in file
    fid=fopen('v1d.temp','wt');
    fprintf(fid,fullfile(folder,file));
    fclose(fid);
else
    if exist('.v1d.temp') % read last opened folder from temp.temp
        fid=fopen('.v1d.temp','r');
        if fid~=-1
            fn=textscan(fid,'%s');
        else
            fn{1}=[];
        end
        fclose(fid);
        if ~isempty(fn{1})
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks',fn{1}{1});
        else
            [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
        end
    else
        [file,folder]=uigetfile('*.txt','Choose *.txt file with picks');
    end
    % save last selected folder in file
    fid=fopen('.v1d.temp','wt');
    fprintf(fid,fullfile(folder,file));
    fclose(fid);
end



fid=fopen(fullfile(folder,file),'r');
temp=textscan(fid,'%f%f%f%f%f%f','Headerlines',1);
fclose(fid);
data=[temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}]; % Profnum channel x y t v


% plot 3d
figure
scatter3(data(:,3),data(:,4),data(:,5),30,data(:,6),'fill')
xlabel('x [m]')
ylabel('y [m]')
zlabel('t [ns]')
set(gca,'FontSize',20)
grid on
cb=colorbar;
ylabel(cb,'vrms [cm/ns]')
set(cb,'FontSize',20)


% plot 2d
figure
scatter(data(:,3),data(:,4),data(:,5).*4,data(:,6),'fill')
grid on
xlabel('x [m]')
ylabel('y [m]')
set(gca,'FontSize',20,'Dataaspectratio',[1 1 1])
cb=colorbar;
ylabel(cb,'vrms [cm/ns]')
set(cb,'FontSize',20)

figure
plot(data(:,5),data(:,6),'*','Linewidth',2)
hold on
tv=0:ceil(max(data(:,5)));
if fit==1
    p=polyfit(data(:,5),data(:,6),ord);
    plot(tv,polyval(p,tv),'r')
    vrms=polyval(p,tv);
elseif fit==2
    fun=@(a,x) a(1)*exp(a(2)*x)+a(3);
    f=lsqcurvefit(fun,[100 -1 1],data(:,5),data(:,6));
    plot(tv,fun(f,tv),'r')
    vrms=fun(f,tv);
end
ylabel('vrms [cm/ns]')
xlabel('t [ns]')
set(gca,'Fontsize',20)
grid on
axis ij

% save for migration
save(fullfile(folder,'vrms.mat'),'vrms','-v7.3');
save(fullfile(folder,'tv.mat'),'tv','-v7.3');
disp('1D-function saved as vrms.mat and tv.mat!')