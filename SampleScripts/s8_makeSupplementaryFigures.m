%% s8_makeSupplementaryFigures
% Reproducing supplementary figures from "Spatial Tuning of Alpha Oscillations in Human Visual Cortex"

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
makeFigureS1;
%%
makeFigureS2;
%%
makeFigureS3_S4;
%%
makeFigureS5;
%%
makeFigureS6;
%%
makeFigureS7;
%%
makeFigureS8;

%% Finish session
close all;
