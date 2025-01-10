%% s7_makeFigures
% Reproducing figures from "Precise Spatial Tuning of Visually Driven Alpha Oscillations in Human Visual Cortex"

% 20220223 Yuasa
% 20230712 Yuasa - update figure numbers
% 20240830 Yuasa - update figure numbers for revision

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
makeFigure4_4S1_4S2;
%%
makeFigure5_5S1;
%%
makeFigure6_6S2;
%%
makeFigure7_7S1;
%%
makeFigure8;
%%
% Figure9 is not produced by MATLAB.
%%
makeFigure10;
%%
makeFigure11;

%% Finish session
close all;
