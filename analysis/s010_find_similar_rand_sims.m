%load('./compiled_data/s008_giant_pca.mat')
close all

min_ss_all = inf;
for ind = 1:1127
  not_ind = true(1128,1);
  not_ind(end) = false; % last one is the target params
  not_ind(ind) = false;
  [min_ss_ind,best_sim] = min(sum((Score(not_ind,1:end) - repmat(Score(ind,1:end),1126,1)).^2,2));
  if ind <= best_sim
    best_sim = best_sim+1;
  end
  if min_ss_all > min_ss_ind
    closest_pair = [ind, best_sim];
    min_ss_all = min_ss_ind;
  end
end

figure; hold on
plot(Score(closest_pair(1),1:50),'ob','MarkerFaceColor','b')
plot(Score(closest_pair(2),1:50),'xr','LineWidth',1.5,'MarkerSize',8)
plot(quantile((Score(:,1:50)),0.95),'-m')
plot(quantile((Score(:,1:50)),0.05),'-m')
title([num2str(closest_pair(1)) '//' num2str(closest_pair(2))])
legend(['Sim ' num2str(closest_pair(1))],['Sim ' num2str(closest_pair(2))])
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')

figure; hold on
plot(Score(closest_pair(1),1:end),'ob','MarkerFaceColor','b')
plot(Score(closest_pair(2),1:end),'xr','LineWidth',1.5,'MarkerSize',8)
plot(quantile((Score(:,1:end)),0.95),'-m')
plot(quantile((Score(:,1:end)),0.05),'-m')
title([num2str(closest_pair(1)) '//' num2str(closest_pair(2))])
legend(['Sim ' num2str(closest_pair(1))],['Sim ' num2str(closest_pair(2))])
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')


figure; hold on
plot(pV_stacked(:,closest_pair(1)),'ob','MarkerFaceColor','b')
plot(pV_stacked(:,closest_pair(2)),'xr','LineWidth',1.5,'MarkerSize',8)
plot([0 30],[1 1],':k')
title([num2str(closest_pair(1)) '//' num2str(closest_pair(2))])
legend(['Sim ' num2str(closest_pair(1))],['Sim ' num2str(closest_pair(2))])
xlabel('Model Parameters')
ylabel('Proportional Change')

% yy_all = inf;
% for ind = 1:1127
% %   if ind == 927 || ind == 125
% %     continue
% %   end
%   not_ind = true(1128,1);
%   not_ind(end) = false; % last one is the target params
%   not_ind(ind) = false;
%   [min_ss_ind,best_sim] = min(sum((Score(not_ind,1:50) - repmat(Score(ind,1:50),1126,1)).^2,2));
%   ss_ind = sum((pV_stacked(:,best_sim) - pV_stacked(:,ind)).^2);
%   yy = (min_ss_ind - 0*ss_ind);
%   if yy_all > yy
%     closest_pair = [ind, best_sim];
%     yy_all = yy;
%   end
% end
% 
% figure; hold on
% plot(Score(closest_pair(1),1:50),'ob','MarkerFaceColor','b')
% plot(Score(closest_pair(2),1:50),'xr','LineWidth',1.5,'MarkerSize',8)
% plot(quantile((Score(:,1:50)),0.95),'-m')
% plot(quantile((Score(:,1:50)),0.05),'-m')
% title([num2str(closest_pair(1)) '//' num2str(closest_pair(2))])
% legend('Sim1','Sim2')
% 
% figure; hold on
% plot(pV_stacked(:,closest_pair(1)),'ob','MarkerFaceColor','b')
% plot(pV_stacked(:,closest_pair(2)),'xr','LineWidth',1.5,'MarkerSize',8)
% title([num2str(closest_pair(1)) '//' num2str(closest_pair(2))])
% legend('Sim1','Sim2')


% yy_all = inf;
% for ind = 1:1127
%   not_ind = true(1128,1);
%   not_ind(end) = false; % last one is the target params
%   not_ind(ind) = false;
%   [max_ss_ind,best_sim] = max(sum((Score(not_ind,1:50) - repmat(Score(ind,1:50),1126,1)).^2,2));
%   ss_ind = sum((pV_stacked(:,closest_pair(1)) - pV_stacked(:,closest_pair(2))).^2);
%   yy = (-max_ss_ind + 0.5*ss_ind);
%   if yy_all > yy
%     closest_pair = [ind, best_sim];
%     yy_all = yy;
%   end
% end
% 
% figure; hold on
% plot(Score(closest_pair(1),1:50),'ob','MarkerFaceColor','b')
% plot(Score(closest_pair(2),1:50),'xr','LineWidth',1.5,'MarkerSize',8)
% plot(quantile((Score(:,1:50)),0.95),'-m')
% plot(quantile((Score(:,1:50)),0.05),'-m')
% title([num2str(ind) '//' num2str(best_sim)])
% legend('Sim1','Sim2')
% 
% figure; hold on
% plot(pV_stacked(:,closest_pair(1)),'ob','MarkerFaceColor','b')
% plot(pV_stacked(:,closest_pair(2)),'xr','LineWidth',1.5,'MarkerSize',8)
% title([num2str(ind) '//' num2str(best_sim)])
% legend('Sim1','Sim2')