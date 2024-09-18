clear all
close all
clc

%------------------------ CMP_Processing ----------------------------------
%
% Dr. Tina Wunderlich, September 2024, tina.wunderlich@ifg.uni-kiel.de, CAU Kiel
%
% For using this script you need your data in the MultichannelGPR-format
% (radargrams.mat, t.mat, global_coords.mat, x.mat), but with no
% coordinates in x.mat and glbal_coords.mat. For GSSI-data use
% DZT_Convert.m first.
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% DO NOT CHANGE FROM HERE ON!
warning('off');


% set path temporarily:
oldpath=path;
currentFile = pwd;
curFold=fileparts(currentFile);
addpath(fullfile(curFold,'GUIs'));

% Start GUI
CMPproc_exported;

waitfor(findall(0,'type','figure'));

% set original path
path(oldpath);