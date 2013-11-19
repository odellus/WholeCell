% alex_update_single_table: This script updates the table in
%                        "single-whole-cell-table.mat" which contains
%                        all of the previous simulation results with
%                        nTrials==1. A different table holds the results of
%                        simulations with nTrials>1;
%                        "averaged-whole-cell-table.mat"
%
% REMARKS:
%  -- This script calculates the "prediction distance" for each simulation.
%  -- The whole-cell-table is sorted by the "prediction distance"
%  -- This script is called by alex_runWholeCell_inf after a simulation
%       finishes.

load single-whole-cell-table.mat
w = what('./output/single/');

for a = 1:length(w.mat)
  if ~any(strcmp(sing_cell_table.filename,w.mat{a}))
    % Calculate prediction distance
    est = load(['./output/single/' w.mat{a}]);
    if isfield(est,'perturbVec')
      tru = load('./mutant_WT_avgData');
      predDist = alex_calcPredDist(tru,est);
      est.predDist = predDist;
      save(['./output/single/' w.mat{a}],'predDist','-append');
      % Incorporate new information into the table
      sing_cell_table.perturbVecs = [sing_cell_table.perturbVecs est.perturbVec];
      sing_cell_table.predDist = [sing_cell_table.predDist; est.predDist];
      sing_cell_table.number = [sing_cell_table.number; length(sing_cell_table.number)+1];
      sing_cell_table.filename{end+1,1} = w.mat{a};
      % Sort everything based on whole_cell_table.predDist
      [~,I] = sort(sing_cell_table.predDist);
      sing_cell_table.predDist = sing_cell_table.predDist(I);
      sing_cell_table.perturbVecs = sing_cell_table.perturbVecs(:,I);
      sing_cell_table.number = sing_cell_table.number(I);
      sing_cell_table.filename = sing_cell_table.filename(I);
    end
  end
end
save('single-whole-cell-table.mat','sing_cell_table');
