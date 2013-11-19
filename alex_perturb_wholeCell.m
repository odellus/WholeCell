function [fns,avg_fn] = alex_perturb_wholeCell(perturbVec, nTrials)
% alex_perturb_wholeCell: Runs a simulation of the whole cell model. Each
%                         of the 30 free parameters are perturbed by the
%                         entries given in "perturbVec"
%
% INPUT:  "perturbVec" - 30x1 vector, which specifies how each of the 30
%                        free parameters are perturbed. A value of 1.0
%                        means that the entry is unperturbed. To meet the
%                        criteria listed in the competition description,
%                        all entries should be on the interval [0.0660,
%                        1.906].
%         "nTrials"    - Integer to specify the number of trials to run.
%                        If nTrials > 1, then each simulation is run with a
%                        different seed.
%
% OUTPUT: "fns"    - Cell array holding the filenames of each simulation.
%                    These are saved in the directory "./output". If
%                    nTrials == 1, they are saved in the sub-directory
%                    "single"; if nTrials > 1, they are saved in the
%                    sub-directory "multiple".
%         "avg_fn" - String holding the file holding the averaged
%                    simulation data (across nTrials). If nTrials == 1,
%                    then avg_fn = nan.
%
% Author: Alex Williams

assert(length(perturbVec) == 30);
assert(all(perturbVec > 0));
assert(nTrials >= 1);

rng('shuffle')

genesTusRxns = {
    'MG_006' 'TU_003' 'Tmk'
    'MG_023' 'TU_011' 'Fba'
    'MG_047' 'TU_027' 'MetK'
    'MG_111' 'TU_069' 'Pgi'
    'MG_272' 'TU_180' 'AceE'
    'MG_299' 'TU_203' 'Pta'
    'MG_330' 'TU_233' 'CmkA2'
    'MG_357' 'TU_260' 'AckA'
    'MG_407' 'TU_294' 'Eno'
    'MG_431' 'TU_307' 'TpiA'
    };

setPath();
setWarnings();
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

% The index "a" keeps track of the parameter (within a gene) that is being
% perturbed {order = 'promoter affinity','rna half life','kcat'}. This
% index is incremented by three every iteration of the loop.
a = 1;
for i_gene = 1:10
  tuId = genesTusRxns{i_gene, 2};
  rxnId = genesTusRxns{i_gene, 3};
  % Promotor Affinity Parameter
  sim.applyRnaPolTuBindingProbs(struct(tuId,...
                               perturbVec(a)*rnaPolTuBindingProbs.(tuId)));
  % Half Life Parameter
  sim.applyRnaHalfLives(struct(tuId,...
                               perturbVec(a+1)*rnaHalfLives.(tuId)));
  % Reaction kCat Parameter
  sim.applyMetabolicReactionKinetics(struct(rxnId,...
                  struct('for', perturbVec(a+2)*rxnKinetics.(rxnId).for)));
  a = a+3;
end

% Run the whole-cell simulation three times and save in "output"
if nTrials == 1
  outDir = './output/single'; % save single trials here
else
  outDir = './output/multiple'; % save replicate trials here
end

dirLen = length(dir(outDir))-2;
tag = randi(1e7);
fns = cell(nTrials,1);
for trial = 1:nTrials
  seed = randi(220516568);
  sim.applyOptions('lengthSec', 65000, 'seed', seed);
  
  parameterVals = sim.getAllParameters();
  
  fns{trial} = sprintf('whole-cell-sim-%06i-%08i-%09i-%i.mat',dirLen,tag,seed,trial);
  simulateHighthroughputExperiments(...
    'seed', seed, ...
    'parameterVals', parameterVals, ...
    'simPath', [outDir '/' fns{trial}] ...
    );
end

if nTrials > 1
  % Average the three simulations and save in "averaged_output"
  averageHighthroughputExperiments(...
    'simPathPattern', sprintf('output/multiple/whole-cell-sim-%06i-%08i-*.mat',dirLen,tag), ...
    'avgValsPath', sprintf('averaged_output/averaged_sim-%06i-%08i.mat',dirLen,tag) ...
    );
  
  % Append perturb vector to file
  save(sprintf('./averaged_output/averaged_sim-%06i-%08i.mat',dirLen,tag),'perturbVec','-append')
  avg_fn = sprintf('averaged_sim-%06i-%08i.mat',dirLen,tag);
else
  avg_fn = nan;
end