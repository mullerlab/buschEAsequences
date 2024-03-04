function [ SDFs ] = sdfMaker( rasters, stdev )
%% use convolution theorem to get continuous estimates of firing rates
% code modified from Analyzing Neural Time Series by Mike X Cohen
% input: rasters should be size nUnits x Time x numTrials
% output: SDFs

% Create Gaussian kernel
time = -1:1/1000:1; %each time point is separated by 1ms
ss = stdev; %15ms standard deviation for gaussian window
gaussianWin = 1/(sqrt(2*ss^2*pi)) * exp(-time.^2./(2*ss^2));
% half of the wavelet size, useful for chopping off edges after convolution.
halfWinSize = (length(gaussianWin)-1)/2;

numUnits = size(rasters,1);

numPts = size(rasters,2); %number of time points
numTrials = size(rasters,3); % number of trials for this session
nData = numPts * numTrials; %concatenate all trials to make it run faster

%preallocate array
SDFs = zeros([numUnits,numPts,numTrials]);

%take the fft of the window
n_conv = length(gaussianWin) + nData - 1;
n_conv_pow2 = pow2(nextpow2(n_conv)); %for fft purposes, calculate next power of 2
fftWin = fft(gaussianWin,n_conv_pow2);

for unit = 1:numUnits
    
    % do convolution take the fft of the signal
    fftSig = fft(reshape(squeeze(rasters(unit,:,:)),1,nData),n_conv_pow2);
    
    %Convolution theorem
    unitSDF = ifft(fftWin.*fftSig);
    unitSDF = unitSDF(1:n_conv);
    unitSDF = unitSDF(halfWinSize+1:end-halfWinSize);
    convData = reshape(unitSDF,numPts,numTrials);
    
    SDFs(unit,:,:) = convData;
    
    
end

end