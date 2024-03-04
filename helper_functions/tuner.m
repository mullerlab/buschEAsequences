function [ tuning ] = tuner( rasterData, conds, win,win2)
%TUNER returns tuning for each unit in rasterData
% rasterData is numUnits x numTime x numTrials
%   uses Kruskal-wallis to assess difference between mdn spike counts
%counts spikes in window win, then throws into kruskalwallis to determine
%if sign. different and returns max meanrank

[numUnits,numTime,numTrials]=size(rasterData);

spikeCounts = zeros(numUnits,numTrials);

for ii = 1:numTrials
    spikeCounts(:,ii) = squeeze(sum(rasterData(:,win(ii):win2(ii),ii),2));
end

% spikeCounts = squeeze(sum(rasterData(:,win,:),2));
duration=mean((win2-win)/1000);
pVal = zeros(numUnits,1);
maxCond = zeros(numUnits,1);
allConds = unique(conds);

numConds = numel(allConds);
medFR = zeros(numUnits,numConds);
meanFR = zeros(numUnits,numConds);

for unit = 1:numUnits
    
    [pVal2(unit),~,stats]=kruskalwallis(spikeCounts(unit,:),conds,'off');
    [pVal(unit),terms,stats2]=anovan(spikeCounts(unit,:),{conds},'display','off');
    SSgroup(unit,1)=terms{2,2};
    SStotal(unit,1)=terms{4,2};
    df(unit,1)=terms{2,3};
    mse(unit,1)=stats2.mse;
  
    
    for cond = 1:numConds  
        medFR(unit,cond) = median(spikeCounts(unit,conds == allConds(cond)));
        meanFR(unit,cond) = mean(spikeCounts(unit,conds == allConds(cond)));
    end
    
    [~,maxInd] = max(stats.meanranks);
    maxCond(unit) = allConds(maxInd);
    [~,minInd] = min(stats.meanranks);
    minCond(unit) = allConds(minInd);
    
end

 tuning.pVal = pVal;
 tuning.SSgroup=SSgroup;
 tuning.SStotal=SStotal;
 tuning.mse=mse;
 tuning.df=df;
%  tuning.maxCond = maxCond;
%  tuning.minCond = minCond';
  tuning.gNames = stats.gnames;
%  tuning.nTrials = stats.n;
 % tuning.nTrials=numTrials;
tuning.medFR = medFR/duration;
tuning.meanFR = meanFR/duration;
end