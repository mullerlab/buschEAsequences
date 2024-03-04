function phases = get_epoch_sequences(spiketimes, epoch, delayOnset, delayOffset)

% INPUTS
% spiketimes = cells x trials matrix of spike sequences for corresponding trials 
% epoch = 'cue' , 'delay', 'nav'
% [delayOnset, delayOffset] = cue-delay and delay-nav boundary times 

% OUTPUT
% updated spiketimes struct
phases = spiketimes; 
if strcmp(epoch, 'cue')
    phases(spiketimes >= delayOnset) = NaN;
elseif strcmp(epoch, 'delay')
    phases( (spiketimes < delayOnset) | (spiketimes > delayOffset)) = NaN;
elseif strcmp(epoch, 'nav')
    phases(spiketimes<=delayOffset) = NaN;
elseif strcmp(epoch, 'delnav')
    phases( spiketimes < delayOnset) = NaN;
end

