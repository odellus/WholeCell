close all; clear; clc
set(0,'DefaultFigureWindowStyle','docked')

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'
             'rxnFluxes'};

% Calculate principal components for each datatype and save them into the
% "compiled_data" directory
for a = 1:length(datatypes)
  disp(a)
  f000_predictPerturbs_pairwiseCorr(datatypes{a})
end

% Run the script to build estimates for perturbVec using random forest. 
s006_improved_random_forest_est