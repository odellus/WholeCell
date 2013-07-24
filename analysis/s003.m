% Need all the f000_*.mat files to be made. Can be done by running s002
% script.

clear; close all; clc

datatypes = {'metConcs'
             'dnaSeq'
             'rnaSeq'
             'chipSeq'
             'rnaArray'
             'protArray'}; % 'rxnFluxes' left out, zscore is bad..

all_data_stacked = [];
for a = 1:length(datatypes)
  S = load(['f000_' datatypes{a} '.mat']);
  all_data_stacked = cat(2,all_data_stacked,S.data_stacked);
  %sPCA.(datatypes{a}).Score = S.Score(:,1:25);
  %sPCA.(datatypes{a}).Score = S.Explained(1:25);
  pV_stacked = S.pV_stacked;
end

data_normlzd = zscore(all_data_stacked,0,2);
[Coeff,Score,~,~,Explained] = pca(data_normlzd);


% % Figure suggesting that 5 principal components and ~450 trees is enough
% figure; hold all
% for nComp = [5 10 15 25 50 100]
%   X = Score(1:(end-1),1:nComp);
%   y = pV_stacked(3,:)';
%   
%   opts = statset('UseParallel','always');
%   B = TreeBagger(1e3,X,y,'Method','regression','oobvarimp','on','Options',opts);
%   plot(oobError(B))
%   xlabel('number of grown trees')
%   ylabel('out-of-bag classification error')
% end
% legend('5','10','15','25','50','100')

% Figure showing our relative ability to predict each variable:
figure; hold on
nComp = 10;
nTrees = 300;
perturbVec_est = nan(30,1);
opts = statset('UseParallel','always');
for param = 1:30
  y = pV_stacked(param,:)';
  y_est = nan(15,1);
  for trial = 1:15
    PCs = datasample(1:size(Score,2),nComp,'Replace',false,'Weights',Explained);
    X = Score(1:(end-1),PCs);
    x_pred = Score(end,PCs);
  
    B = TreeBagger(nTrees,X,y,'Method','regression','Options',opts);
    y_est(trial) = B.predict(x_pred);
  end
  plot(param,y_est,'ob')
  perturbVec_est(param) = mean(y_est);
end
xlabel('Different Parameters')
ylabel('Estimate of Perturbation')
save('s003_perturbVec_est','perturbVec_est')
