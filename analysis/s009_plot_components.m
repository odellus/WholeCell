load('./compiled_data/s008_giant_pca.mat')
close all

[~,best_sim] = min(sum((Score(1:(end-1),1:10) - repmat(Score(end,1:10),1127,1)).^2,2));

figure; hold on
plot(Score(end,1:50),'ob','MarkerFaceColor','b')
plot(Score(best_sim,1:50),'xr','LineWidth',1.5,'MarkerSize',8)
plot(quantile((Score(:,1:50)),0.95),'-m')
plot(quantile((Score(:,1:50)),0.05),'-m')
legend('Target','Closest Simulation','95% range')
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')

figure; hold on
plot(Score(end,:),'.b')
plot(quantile(Score,0.95),'-m')
plot(quantile(Score,0.05),'-m')
xlim([-20 size(Score,2)])
legend('Target','95% range')
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')

figure; hold on
pc_range = 400:425;
[~,best_sim] = min(sum((Score(1:(end-1),pc_range) - repmat(Score(end,pc_range),1127,1)).^2,2));
plot(pc_range,Score(end,pc_range),'ob','MarkerFaceColor','b')
plot(pc_range,Score(best_sim,pc_range),'xr','LineWidth',1.5,'MarkerSize',8)
plot(pc_range,quantile((Score(:,pc_range)),0.95),'-m')
plot(pc_range,quantile((Score(:,pc_range)),0.05),'-m')
legend('Target','Closest Simulation','95% range')
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')

figure; hold on
pc_range = 400:425;
[~,best_sim] = max(sum(abs(Score(1:(end-1),pc_range)),2));
plot(pc_range,Score(end,pc_range),'ob','MarkerFaceColor','b')
plot(pc_range,Score(best_sim,pc_range),'xr','LineWidth',1.5,'MarkerSize',8)
plot(pc_range,quantile((Score(:,pc_range)),0.95),'-m')
plot(pc_range,quantile((Score(:,pc_range)),0.05),'-m')
legend('Target','Closest Simulation','95% range')
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')

tidy_allFigs

figure; hold on
plot(Score(88,:),'.b')
plot(Score(214,:),'xr')
plot(quantile(Score,0.95),'-m')
plot(quantile(Score,0.05),'-m')
xlim([-20 size(Score,2)])
legend('Target','95% range')
xlabel('Principal Components (ranked by variance explained)')
ylabel('PC Score')
