
for zzz = 1:100
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
    %sPCA.(datatypes{a}).Score = S.Score(:,1:25);
    %sPCA.(datatypes{a}).Score = S.Explained(1:25);
    pV_stacked = S.pV_stacked;
  end
  
  minCols = min(all_data_stacked);
  rngCols = range(all_data_stacked);
  nr = size(all_data_stacked,1);
  data_normlzd = zscore(all_data_stacked);
  [Coeff,Score,~,~,Explained] = pca(data_normlzd);
  clear all_data_stacked % reduce junk in memory
  clear data_normlzd
  
  x1 = Score(1:(end-1),1:50); % first 25 principal components for all simulations.
  x2 = pV_stacked';
  x1_pred = Score(end,1:50); % The target 'mutant' phenotype
  x2_pred = nan(1,30);        % List of estimated params
  
  y_est = nan(30,1);
  pv_ind = 1:30;
  a_hist = [];
  toc
  for a = randperm(30)
    disp(a)
    tic
    y = pV_stacked(a,:)';
    x2 = pV_stacked(a_hist,:)';
    B = TreeBagger(5e2,[x1 x2],y,'Method','regression');
    assert(~any(isnan(x2_pred(a_hist))))
    y_est(a) = B.predict([x1_pred x2_pred(a_hist)]);
    x2_pred(a) = y_est(a);
    a_hist = [a_hist; a];
    toc
  end
  
  load('./compiled_data/pv_est_randForest_v3.mat')
  pv_est_randForest_v3(:,end+1) = y_est;
  save('./compiled_data/pv_est_randForest_v3.mat','pv_est_randForest_v3')
end