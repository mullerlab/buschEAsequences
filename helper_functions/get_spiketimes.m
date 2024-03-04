function spikes = get_spiketimes(raster, neuronCond, timeOnset, timeOffset, s)

% INPUTS
% raster = raster data (neurons X time X trials) 
% neuronCond = vector of which cells to include
% timeOnset = start of trial
% timeOffset = end of trial
% s = std for sdf maker 

% OUTPUTS
% spikes = neurons X trials matrix of spike times
% ( entry (i,j) gives time at which cell i reaches max firing rate in trial j) 

trlNum = size(raster, 3); neuronNum = sum(neuronCond); 

% get SDFs
[SDFs] = sdfMaker(raster(neuronCond,timeOnset:timeOffset,:), s);

% get spike times = times of max firing rate 
spikes = zeros( neuronNum, trlNum ); 
for trl = 1: trlNum
    t = SDFs(:,:,trl); % get sdf for that trial 
    [~, maxT] = max(t, [], 2); % get time at which max firing rate occurs 
    spikes(:,trl) = maxT; % save spike time
end

