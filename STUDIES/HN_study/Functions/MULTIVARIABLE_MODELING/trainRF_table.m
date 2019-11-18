function [RF] = trainRF_table(table,Y,categories,nBoot,cost)

[trainSets,~] = buildBootSet(Y,nBoot,'NoAdjust');

RF = applyEnsembleRF_table(table(trainSets(:,1),:),Y(trainSets(:,1)),categories,cost);
for n = 2:nBoot
   RFtemp = applyEnsembleRF_table(table(trainSets(:,n),:),Y(trainSets(:,n)),categories,cost);
   RF = append(RF,RFtemp);
end

end