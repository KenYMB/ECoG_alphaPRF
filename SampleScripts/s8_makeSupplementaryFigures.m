%% s8_makeSupplementaryFigures
% Reproducing supplementary figures from "Precise Spatial Tuning of Visually Driven Alpha Oscillations in Human Visual Cortex"

% 20220223 Yuasa
% 20230712 Yuasa - update figure numbers

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
%%
% FigureS9 is the schematic.
%%
makeFigureS10;

%% Finish session
close all;
