% Code to produce manuscript figures 

clearvars; clc; 

monkeyT = (1/255)*[216,216,216]; 
monkeyB = (1/255)*[129,213,218];

load('./processed_data/fr1_wm_NSP0-1_neuralboundary_delay.mat')
load('./processed_data/fr1_wm_NSP0-1_cell_selectivity.mat')

%% F3b,f,g : example centroids , distance matrices 

ff = 15; 

figure; scatter3( session_info(ff).all_cells.correct.centroids(:,1), session_info(ff).all_cells.correct.centroids(:,2), ...
    session_info(ff).all_cells.correct.centroids(:,3), 150, 1:9, 'filled' )
title( 'Figure 3b')

figure; imagesc( session_info(ff).all_cells.correct.dm); colormap bone; 
c = colorbar; c.Label.String = 'Normalized Euclidean Distance';
xlabel('Condition'); ylabel('Condition')
title( 'Figure 3f')

load('./processed_data/fr1_wm_trajectories.mat')
map = [12:17, 1:11]; 
figure; imagesc( session_data(map(ff)).traj_dm); colormap bone; 
c = colorbar; c.Label.String = 'Normalized Frechet Distance';
xlabel('Condition'); ylabel('Condition')
title( 'Figure 3g') 

%% F3e : neural boundaries 

store = struct; 
for ff = 1: 17 
    % observed values 
    times = S(ff).times(S(ff).time_selective,:); 
    times = reshape( times, 1, []);
    store(ff).time = times;
end

% shuffle times
shuff = 7000; % max number of ms to shuffle spike times within trial 
iterations = 100; % number of shuffles 
bins = 0:250:7000 ;

counts = zeros( 100, length(bins)-1); 
for it = 1:iterations
    shuffle_times = [store.time] + (rand(size([store.time]))-0.5).*shuff; 
    shuffle_times = mod(shuffle_times, 7000); 
    counts(it, :) = histcounts( shuffle_times, bins)./length(shuffle_times);
end

figure; hold on 
B = histcounts( [store(7:17).time], bins); B = B./length([store(7:17).time]); 
T = histcounts( [store(1:6).time], bins); T = T./length([store(1:6).time]); 
plot(bins(1:end-1),mean(counts,1), 'color', monkeyT, 'linewidth', 2)
bar( bins(1:end-1),B, 'FaceColor', monkeyB, 'FaceAlpha', 1); 
bar(bins(1:end-1),T, 'FaceColor', monkeyT, 'FaceAlpha', 0.7)
plot([3000,3000], [0,0.15], '-k'); plot([5000,5000], [0,0.15], '-k')
legend( 'Shuffle','NHP B', 'NHP T' ); 
xlabel( 'spike time (ms)')
ylabel( 'fraction of trials')
title( 'Figure 3e') 

%% F3h : correct vs incorrect correlation 

trls = [1,2,3,4,5,6,7,8,9,11,12] ; 
% all of these have all 9 targets represented in both correct and incorrect. 
% 6 sessions for NHP T, 5 for NHP B
correct = zeros(1, length(trls)); 
incorrect = zeros(1, length(trls)); 

figure; hold on
for ff = 1:6 
    plot([1,3], [session_info(trls(ff)).all_cells.correct.cc, session_info(trls(ff)).all_cells.incorrect.cc], 'color', [0.5,0.5,0.5])
    p1 = plot(1, session_info(trls(ff)).all_cells.correct.cc, '.', 'Markersize', 20, 'color', monkeyT);
    plot(3, session_info(trls(ff)).all_cells.incorrect.cc, '.', 'Markersize', 20, 'color', monkeyT)
    correct(ff) = session_info(trls(ff)).all_cells.correct.cc;
    incorrect(ff) = session_info(trls(ff)).all_cells.incorrect.cc; 
end
for ff = 7:length(trls) 
    plot([1,3], [session_info(trls(ff)).all_cells.correct.cc, session_info(trls(ff)).all_cells.incorrect.cc], 'color', [0.5,0.5,0.5])
    p2 = plot(1, session_info(trls(ff)).all_cells.correct.cc, '.', 'Markersize', 20, 'color', monkeyB);
    plot(3, session_info(trls(ff)).all_cells.incorrect.cc, '.', 'Markersize', 20, 'color', monkeyB)
    correct(ff) = session_info(trls(ff)).all_cells.correct.cc;
    incorrect(ff) = session_info(trls(ff)).all_cells.incorrect.cc; 
end
p3 = plot( [0.8,1.2], [mean(correct), mean(correct)], 'k', 'Linewidth', 1);
plot([2.8,3.2], [mean(incorrect),mean(incorrect)], 'k', 'Linewidth', 1); 
ylim([0,1]); ylabel('Spearman Correlation');
xlim([0,4]); title('Centroid-Target Distance Matrix Correlation by Session')
xticks([1,3]); xticklabels({'Correct', 'Incorrect'}); 
legend( [p1,p2, p3], {'NHP T','NHP B', 'Mean'})
title( 'Figure 3h') 


%% 4h VR vs ODR

% load ODR struct
load('./processed_data/ODRCorrelationData_2023.mat') 
ODR = zeros( 1, length(session_info));
for ff = 1:length(session_info)
    ODR(ff) = session_info(ff).all_cells.correct.cc; 
end

% load VR struct 
load('./processed_data/fr1_wm_NSP0-1_taskboundary_delay.mat')
VR = zeros( 1, length(session_info)); 
for ff = 1:length(session_info)
    VR(ff) = session_info(ff).all_cells.correct.cc;
end

% correlation plot
figure; hold on; 
scatter( 0.2*(rand(1,length(VR))-0.5)+1, VR(1,:), 50, 'blue', 'filled');
plot( [0.8, 1.2], [median(VR(1,:)), median(VR(1,:))], 'k', 'Linewidth', 1);
scatter( 0.2*(rand(1,length(ODR))-0.5)+2, ODR(1,:), 50, 'red', 'filled');
plot( [1.8, 2.2], [median(ODR(1,:)), median(ODR(1,:))], 'k', 'Linewidth', 1);
set(gca, 'XTick', [1,2], 'XTickLabels', {'VR Task', 'ODR Task'})
ylabel( 'Spearman Correlation')
xlim([0,3])
ylim([0,1])
title( 'Figure 4h') 

%% Figure 5: Ketamine 

% fig 5 e) 
clearvars; clc
load('./processed_data/fr1_preKet_NSP0-1_neuralboundary_delay.mat')
pre = session_info;
load('./processed_data/fr1_postKet_NSP0-1_neuralboundary_delay.mat')
post = session_info;
load('./processed_data/fr1_post30Ket_NSP0-1_neuralboundary_delay.mat')
latepost = session_info;

preKet = zeros(1,17); post30Ket = zeros(1,17); postKet = zeros(1,17);
figure; hold on
for ff = 1:17
    preKet(ff)=pre(ff).all_cells.correct.cc;
    postKet(ff)=post(ff).all_cells.correct.cc;
    post30Ket(ff)=latepost(ff).all_cells.correct.cc;

    plot( [1,2,3], [preKet(ff), postKet(ff), post30Ket(ff)], 'color', (1/255)*[216,216,216])
end
% EXCLUDE SESSION 2 --- NO CORRECT TRIALS FOR TARGETS 1,2 
preKet(2) = []; postKet(2) = []; post30Ket(2) = []; 
o = ones( 1,16); 
plot( o(1:16), preKet, 'g.', 'Markersize', 20);
plot( 2*o, postKet, 'b.', 'Markersize', 20);
plot( 3*o, post30Ket, 'r.', 'Markersize', 20);
plot( [0.9, 1.1], [mean(preKet), mean(preKet)], '-k', 'Linewidth', 2)
plot( [1.9, 2.1], [mean(postKet), mean(postKet)], '-k', 'Linewidth', 2)
plot( [2.9, 3.1], [mean(post30Ket), mean(post30Ket)], '-k', 'Linewidth', 2)
ylim([0,1]); ylabel( 'Spearman Correlation')
set(gca, 'XTick', 1:3, 'XTickLabels', {'Pre', 'Post', 'Late-Post'})
title( 'Figure 5e') 

% 5d projections

ff=14;
% pre ketamine centroids
centroids = pre(ff).all_cells.correct.centroids; 
figure; scatter3( centroids(:,1), centroids(:,2), centroids(:,3), 150, 1:9, 'filled')
title( 'Figure 5d top') 
% post ketamine centroids 
centroids = post(ff).all_cells.correct.centroids; 
figure; scatter3( centroids(:,1), centroids(:,2), centroids(:,3), 150, 1:9, 'filled')
title( 'Figure 5d middle') 
% late-post ketamine centroids 
centroids = latepost(ff).all_cells.correct.centroids; 
figure; scatter3( centroids(:,1), centroids(:,2), centroids(:,3), 150, 1:9, 'filled')
title( 'Figure 5d bottom') 
