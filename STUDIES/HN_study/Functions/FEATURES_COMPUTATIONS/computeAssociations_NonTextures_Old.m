function computeAssociations_NonTextures_Old(pathUnivariate,pathFeatures,cohorts,outcomes,timeToEvent,featType,nonTextName)

startpath = pwd;

nameOutcomes = fieldnames(outcomes.(cohorts{1})); nOutcomes = numel(nameOutcomes); nCohort = numel(cohorts);
nFeatType = numel(featType); nNonText = numel(nonTextName);

cd(pathFeatures)
load('nonText_HGJ_GTVp'), nonTextVar = cell(nNonText*nFeatType,1);
rsMat_Binary = zeros(nNonText*nFeatType,nOutcomes); pMat_Binary = zeros(nNonText*nFeatType,nOutcomes); stringCell_Binary = cell(nNonText*nFeatType,nOutcomes);
rsMat_Time = zeros(nNonText*nFeatType,nOutcomes); pMat_Time = zeros(nNonText*nFeatType,nOutcomes); stringCell_Time = cell(nNonText*nFeatType,nOutcomes);
count = 0;
for i = 1:nNonText
    for j = 1:nFeatType
        count = count + 1;
        nonTextVar{count} = [nonTextName{i},'_',featType{j}];
        for o = 1:nOutcomes
            var = []; outcome = []; time = [];
            for c = 1:nCohort
                cohort = cohorts{c};
                load(['nonText_',cohort,'_',featType{j}])
                var = [var;nonText.(nonTextName{i}).Data];
                outcome = [outcome;outcomes.(cohort).(nameOutcomes{o})];
                time = [time;timeToEvent.(cohort).(nameOutcomes{o})];
            end
            [rsMat_Binary(count,o),pMat_Binary(count,o)] = corr(var,outcome,'type','Spearman','rows','pairwise');
            [rsMat_Time(count,o),pMat_Time(count,o)] = corr(var,time,'type','Spearman','rows','pairwise');
            if pMat_Binary(count,o) < 0.01
                stringCell_Binary{count,o} = ['rs = ',num2str(rsMat_Binary(count,o),'%.2f'),', p = ',num2str(pMat_Binary(count,o),'%.2i')];
            else
                stringCell_Binary{count,o} = ['rs = ',num2str(rsMat_Binary(count,o),'%.2f'),', p = ',num2str(pMat_Binary(count,o),'%.2f')];
            end
            if pMat_Time(count,o) < 0.01
                stringCell_Time{count,o} = ['rs = ',num2str(rsMat_Time(count,o),'%.2f'),', p = ',num2str(pMat_Time(count,o),'%.2i')];
            else
                stringCell_Time{count,o} = ['rs = ',num2str(rsMat_Time(count,o),'%.2f'),', p = ',num2str(pMat_Time(count,o),'%.2f')];
            end
        end
    end
end
for j = 1:nOutcomes
    pVal = pMat_Binary(:,j);
    [significance] = benjamini_hochberg(pVal,0.10);
    for i = 1:nNonText*nFeatType
        if significance(i)
            stringCell_Binary{i,j} = ['*',stringCell_Binary{i,j},'*']; % Showing significance;
        end
    end
    pVal = pMat_Time(:,j);
    [significance] = benjamini_hochberg(pVal,0.10);
    for i = 1:nNonText*nFeatType
        if significance(i)
            stringCell_Time{i,j} = ['*',stringCell_Time{i,j},'*']; % Showing significance;
        end
    end
end
cd(pathUnivariate)
Locoregional = stringCell_Binary(:,1); Distant = stringCell_Binary(:,2); Death = stringCell_Binary(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',nonTextVar);
results.rsMat = rsMat_Binary; results.pMat = pMat_Binary;
save('nonText_UniV_Binary','results'), clear results
Locoregional = stringCell_Time(:,1); Distant = stringCell_Time(:,2); Death = stringCell_Time(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',nonTextVar);
results.rsMat = rsMat_Time; results.pMat = pMat_Time;
save('nonText_UniV_Time','results'), clear results

cd(startpath)
end