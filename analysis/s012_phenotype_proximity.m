T = load('C:\Users\Alex\Documents\GitHub\WholeCell\Mutant_WT_avgData.mat');
w = what('C:\Users\Alex\Documents\GitHub\WholeCell\output\');

T_vol = T.singleCell.volume';
T_vol(T_vol==0) = nan;

% min_norm = inf;
% for a = 1:length(w.mat)
%   fn = w.mat{a};
%   E = load(['C:\Users\Alex\Documents\GitHub\WholeCell\output\' fn]);
%   vol_diff = (repmat(E.singleCell.volume',1,32) - T.singleCell.volume').^2;
%   vol_norm = sum(vol_diff(:));
%   if vol_norm < min_norm
%     vol_best = E.singleCell.volume;
%     best_sim = fn;
%     min_norm = vol_norm;
%   end
% end
% 
% figure; hold on
% plot(T_vol,'-b')
% plot(vol_best,'-r','LineWidth',2)
% title(fn)

figure; hold on

for a = 1:length(w.mat)
  fn = w.mat{a};
  E = load(['C:\Users\Alex\Documents\GitHub\WholeCell\output\' fn]);
  clf; hold on
  plot(T_vol,'-b')
  plot(E.singleCell.volume,'-r','LineWidth',2)
  title(fn)
  drawnow
end