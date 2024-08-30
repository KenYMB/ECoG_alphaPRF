%% s8_makeSupplementaryFigures
% Reproducing supplementary figures from "Precise Spatial Tuning of Visually Driven Alpha Oscillations in Human Visual Cortex"

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
makeFigure1S1;
%%
makeFigure3S1;
%%
% makeFigure4_4S1;      % This figure is produced in s7_makeFigures.
%%
makeFigure6S1;
%%
% makeFigure6_6S2;      % This figure is produced in s7_makeFigures.
%%
% makeFigure7_7S1;      % This figure is produced in s7_makeFigures.
%%
makeFigure8S1_8S2;
%%
makeFigure8S3;
%%
makeFigure10S1;

%% Finish session
close all;
