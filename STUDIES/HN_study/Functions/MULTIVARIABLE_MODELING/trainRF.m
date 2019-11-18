function [RF] = trainRF(X,Y,nBoot)

[trainSets,~] = buildBootSet(Y,nBoot,'NoAdjust');

RF = applyEnsembleRF(X(trainSets(:,1),:),Y(trainSets(:,1),1));
for n = 2:nBoot
   RFtemp = applyEnsembleRF(X(trainSets(:,n),:),Y(trainSets(:,n),1));
   RF = append(RF,RFtemp);
end

end