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
makeFigureS2_S3;
%%
makeFigureS4;
%%
makeFigureS5;
%%
makeFigureS6;

%% Finish session
close all;
