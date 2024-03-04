function [centroids, dm, cc, cc_shuff ] = get_correlations( prj, targets, goal_loc, iterations, type)

% INPUTS
% prj = projection matrix
% targets = trial labels for each condition
% goal_loc = distance matrix of targets
% iterations = number of shuffles for null model 
% type = 'Spearman' or 'Pearson' correlation

% OUTPUTS
% centroids = ntarg x ndim matrix of centroids 
% dm = distance matrix (NaN entries for targets that did not appear) 
% cc = observed correlation between centroids distance matrix and targets
% distance matrix 

ndim = size(prj,2); 
u = unique(targets); 
ntarg = length(u); 

centroids = zeros(ntarg, ndim); % get centroids 
for ii = 1:ntarg 
    centroids(ii,:) = mean( prj( targets == u(ii), :), 1, 'omitnan');
end

DM = zeros( ntarg ); % get distance matrix 
for ii = 1:ntarg 
    for jj = 1: ii-1
        DM(ii,jj) = sqrt( sum( (centroids(ii,:) - centroids(jj,:)).^2 ));
    end
end
DM = DM + DM'; 
DM = abs(DM); DM = DM./ max(max(DM)); % normalize 


% need to only consider the targets in goal_loc which actually appear in
% this session (some sessions don't have all incorrect trials represented) 
GL = goal_loc( u, u); 
GL = abs(GL); GL = GL./ max(max(GL)); % normalize 

cc = corr(DM(:), GL(:), 'Type', type); % get correlation

cc_shuff = zeros(1,iterations); 

% SHUFFLE TARGETS 
for it = 1: iterations
    targets = targets( randperm(length(targets))); 
    cent = zeros(ntarg, ndim); % get centroids 
    for ii = 1:ntarg 
        cent(ii,:) = mean( prj( targets == u(ii), :), 1, 'omitnan');
    end
    dist = zeros( ntarg ); % get distance matrix
    for ii = 1:ntarg
        for jj = 1: ii-1
            dist(ii,jj) = sqrt( sum( (cent(ii,:) - cent(jj,:)).^2 ));
        end
    end
    dist = dist + dist'; dist = abs(dist); dist = dist./ max(max(dist));

    cc_shuff(it) = corr(dist(:), GL(:), 'Type', type);
end

dm = NaN*ones( ntarg,ntarg ); dm( u, u) = DM; 