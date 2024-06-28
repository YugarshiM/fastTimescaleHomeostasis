function statStruct = measureTrace(V,tme,tStart,tEnd)
	% just assess if something bursts

	statStruct.depBlock = []; statStruct.tonic = []; statStruct.burst = []; statStruct.silence = []; statStruct.irregBurster = [];
	statStruct.classification.number = NaN;
	statStruct.classification.name = '';


	[pkTops,pkLoc,widthAtHalfMaxPk] = findpeaks(V(tStart:tEnd),tme(tStart:tEnd),'MinPeakProminence',2);

	% we assume this means there's no spiking, so silence or dep block
	if isempty(pkLoc)
		return
	end

	[trsTops,trsLoc,widthAtHalfMaxTrs] = findpeaks(-V(tStart:tEnd),tme(tStart:tEnd),'MinPeakProminence',2);
	trsBots = -trsTops;

	cnt = 0;
	for i = -85:1:-49
		cnt = cnt + 1;
		hypLocsIx = find(trsBots < i);
		storeSeq(cnt,1) = i;
		storeSeq(cnt,2) = length(hypLocsIx);
	end
	% nothing is hyperpolarizing correctly OR neuron is in silence or dep block
	if sum(storeSeq(:,2)) == 0 || length(trsBots) == 0
		return
	end

	botsDetected = find(storeSeq(:,2));
	hypLocsIx = find(trsBots < storeSeq(botsDetected(1),1)); % set the level at the first sign of bottoms!, assume that they're all within a mV
	hypBots = trsBots(hypLocsIx);

	% there's only one min, maybe window needs to be bigger, but doubt it
	if length(hypBots) == 1 || length(hypBots) == 0
		return
	end

	%check number of minima between hypBots.
	period = []; burstInterval = []; IBI = []; spikeHeights = []; hyperPolarizationValues = []; slowWaveAmp = [];
	pksBwBots = [];
	for j = 1:length(hypBots)-1
		tS = trsLoc(hypLocsIx(j)); tE = trsLoc(hypLocsIx(j+1)); 
		period = [period, tE-tS];
		tSeries = V(find(tS==tme):find(tE==tme));
		[pkTops,pkLoc,widthAtHalfMaxPk] = findpeaks(tSeries,'MinPeakProminence',2);
		%findpeaks(tSeries,tme(find(tS==tme):find(tE==tme))-tS,'MinPeakProminence',2);
		%keyboard
		burstTime = pkLoc(end) - pkLoc(1);
		burstInterval = [burstInterval, burstTime];
		IBI = [IBI, (tE-tS)-burstTime];
		[trsBots2,trsLoc2] = findpeaks(-tSeries,'MinPeakProminence',2);
		spikeHeights = [spikeHeights, mean(pkTops) - (-max(trsBots2))];

		hyperPolarizationValues = [hyperPolarizationValues, hypBots(j)];
		slowWaveAmp = [slowWaveAmp, (-max(trsBots2)) - hypBots(j)];
		pksBwBots = [pksBwBots,length(pkTops)];
		
	end


	% first condition checks if the hyperpolarization bottoms are within 3 mV of each other
	% the second condition examines checks if there is more than 1 spike between bots
	% the third condition check periodicity 
	if (range(hypBots) < 3) && prod(pksBwBots > 1) && range(pksBwBots) < 1.5
		statStruct.classification.number = 4;
		statStruct.classification.name = 'Burster';

		statStruct.burst.vals = [mean(period), mean(burstInterval), mean(IBI), mean(spikeHeights), mean(hyperPolarizationValues), mean(slowWaveAmp)];
		%statStruct.burst.rng = [range(period), range(burstInterval), range(IBI), range(spikeHeights), range(hyperPolarizationValues), range(slowWaveAmp)];
		statStruct.burst.cv = [std(period)/mean(period), std(burstInterval)/mean(burstInterval), std(IBI)/mean(IBI), std(spikeHeights)/mean(spikeHeights), std(hyperPolarizationValues)/mean(-1*hyperPolarizationValues), std(slowWaveAmp)/mean(slowWaveAmp)];

		if max(statStruct.burst.cv) > .03
			'there are large variances, pausing'
			statStruct.burst.vals
			statStruct.burst.cv
		end
	end


end