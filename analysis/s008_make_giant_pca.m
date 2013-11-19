clear; close all; clc
addpath('C:\Users\Alex\Documents\GitHub\WholeCell')
datatypes = {'metConcs'
  'dnaSeq'
  'rnaSeq'
  'chipSeq'
  'rnaArray'
  'protArray'
  'rxnFluxes'};

tic

all_data_stacked = [];
for a = 1:length(datatypes)
  S = load(['./compiled_data/f000_' datatypes{a} '.mat']);
  all_data_stacked = cat(2,all_data_stacked,S.data_stacked);
  pV_stacked = S.pV_stacked;
end

minCols = min(all_data_stacked);
rngCols = range(all_data_stacked);
nr = size(all_data_stacked,1);
data_normlzd = zscore(all_data_stacked);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);
clear all_data_stacked % reduce junk in memory
clear data_normlzd

save('./compiled_data/s008_giant_pca','Coeff','Score','Explained','pV_stacked')