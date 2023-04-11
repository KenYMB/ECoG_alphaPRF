%% s7_makeFigures
% Reproducing figures from "Spatial Tuning of Alpha Oscillations in Human Visual Cortex"

% 20220223 Yuasa

%% Initialize
run_checkPath;
clearvars;

%% Directory Setup
% %-- Output path 
% datPth      = fullfile(analysisRootPath,'Data',filesep);
% figPth      = fullfile(analysisRootPath, 'Figures',filesep);
%-- Subject list
subjectList_fname = 'subjectlist.tsv';

%%
makeFigure1;
%%
makeFigure2;
%%
makeFigure3;
%%
makeFigure4;
%%
makeFigure5;
%%
makeFigure6;
%%
makeFigure7;
%%
makeFigure9;
%%
makeFigure10;

%% Finish session
close all;
