%% s7_makeFigures
% Reproducing figures from "Precise Spatial Tuning of Visually Driven Alpha Oscillations in Human Visual Cortex"

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
makeFigure8;
%%
makeFigure9;

%% Finish session
close all;
