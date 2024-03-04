function [raster_correct, raster_incorrect, neuronCond] = get_raster_info(data, array, tuned, pval, timeOnset, timeOffset )

% INPUTS
% data = data struct 
% array = 0 for NPS0, array = 1 for NPS1, array = 2 for both
% tuned = 0 for untuned cells, tuned = 1 for tuned cells, tuned = 2 for both
% pval = threshold for considering tuned neurons
% timeOnset = start of trial
% timeOffset = send of trial

% OUTPUTS
% correct = raster of correct trials 
% incorrect = raster of incorrect trials
% neuronCond = vector of which cells are included

% CORRECT TRIALS
correctTrl=data.WM.trialOutcome;
cond=data.WM.cond(correctTrl);
% create raster
if array == 0
    raster_correct=data.WM.rasters.NSP0(:,:,correctTrl); %neuronXTimeXTrl
elseif array == 1
    raster_correct=data.WM.rasters.NSP1(:,:,correctTrl); %neuronXTimeXTrl
elseif array == 2
    raster_correct = cat(1, data.WM.rasters.NSP0, data.WM.rasters.NSP1);
    raster_correct = raster_correct(:,:, correctTrl);
end
trlNum=size(raster_correct,3);

% get neuron tuning
[tuning] = tuner(raster_correct,cond,repmat(timeOnset,trlNum,1),repmat(timeOffset,trlNum,1));
if tuned == 1 % select tuned neurons
    tuned_neurons = find(tuning.pVal <= pval) ;
    raster_correct = raster_correct(tuned_neurons, :,:);
elseif tuned == 0 % select untuned neurons 
    tuned_neurons = find(tuning.pVal > pval) ;
    raster_correct = raster_correct(tuned_neurons, :,:);
end

% select neurons with firing rate > 1 on all trial conditions (correct trials) 
neuronCond=sum(tuning.meanFR'>1)== length(unique(cond));
if tuned ~= 2
    neuronCond = neuronCond(tuned_neurons); % when selecting for tuned/untuned neurons
end

% INCORRECT TRIALS

correctTrl = ~correctTrl; 
cond=data.WM.cond(correctTrl);

% create raster
if array == 0
    raster_incorrect=data.WM.rasters.NSP0(:,:,correctTrl); %neuronXTimeXTrl
elseif array == 1
    raster_incorrect=data.WM.rasters.NSP1(:,:,correctTrl); %neuronXTimeXTrl
elseif array == 2
    raster_incorrect = cat(1, data.WM.rasters.NSP0, data.WM.rasters.NSP1);
    raster_incorrect = raster_incorrect(:,:, correctTrl);
end
trlNum=size(raster_incorrect,3);

% get neuron tuning
%[tuning] = tuner(raster_incorrect,cond,repmat(timeOnset,trlNum,1),repmat(timeOffset,trlNum,1));
if tuned == 1 % select tuned neurons
    %tuned_neurons = find(tuning.pVal <= pval) ;
    raster_incorrect = raster_incorrect(tuned_neurons, :,:);
elseif tuned == 0 % select untuned neurons 
    %tuned_neurons = find(tuning.pVal > pval) ;
    raster_incorrect = raster_incorrect(tuned_neurons, :,:);
end




end