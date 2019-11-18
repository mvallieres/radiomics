function computeAssociations_Clinical_Old(pathUnivariate,cohorts,clinical,outcomes,timeToEvent)

startpath = pwd;
cd(pathUnivariate)

nameOutcomes = fieldnames(outcomes.(cohorts{1})); nOutcomes = numel(nameOutcomes); nCohort = numel(cohorts);
clinicVar = fieldnames(clinical.(cohorts{1})); nClinic = numel(clinicVar);
rsMat_Binary = zeros(nClinic,nOutcomes); pMat_Binary = zeros(nClinic,nOutcomes); stringCell_Binary = cell(nClinic,nOutcomes);
rsMat_Time = zeros(nClinic,nOutcomes); pMat_Time = zeros(nClinic,nOutcomes); stringCell_Time = cell(nClinic,nOutcomes);
for i = 1:nClinic
    for j = 1:nOutcomes
        var = []; outcome = []; time = [];
        for c = 1:nCohort
            cohort = cohorts{c};
            var = [var;clinical.(cohort).(clinicVar{i})];
            outcome = [outcome;outcomes.(cohort).(nameOutcomes{j})];
            time = [time;timeToEvent.(cohort).(nameOutcomes{j})];
        end
        [rsMat_Binary(i,j),pMat_Binary(i,j)] = corr(var,outcome,'type','Spearman','rows','pairwise');
        [rsMat_Time(i,j),pMat_Time(i,j)] = corr(var,time,'type','Spearman','rows','pairwise');
        if pMat_Binary(i,j) < 0.01
            stringCell_Binary{i,j} = ['rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2i')];
        else
            stringCell_Binary{i,j} = ['rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2f')];
        end
        if pMat_Time(i,j) < 0.01
            stringCell_Time{i,j} = ['rs = ',num2str(rsMat_Time(i,j),'%.2f'),', p = ',num2str(pMat_Time(i,j),'%.2i')];
        else
            stringCell_Time{i,j} = ['rs = ',num2str(rsMat_Time(i,j),'%.2f'),', p = ',num2str(pMat_Time(i,j),'%.2f')];
        end
    end
end
for j = 1:nOutcomes
    pVal = pMat_Binary(:,j);
    [significance] = benjamini_hochberg(pVal,0.10);
    for i = 1:nClinic
        if significance(i)
            stringCell_Binary{i,j} = ['*',stringCell_Binary{i,j},'*']; % Showing significance;
        end
    end
    pVal = pMat_Time(:,j);
    [significance] = benjamini_hochberg(pVal,0.10);
    for i = 1:nClinic
        if significance(i)
            stringCell_Time{i,j} = ['*',stringCell_Time{i,j},'*']; % Showing significance;
        end
    end
end
Locoregional = stringCell_Binary(:,1); Distant = stringCell_Binary(:,2); Death = stringCell_Binary(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',clinicVar);
results.rsMat = rsMat_Binary; results.pMat = pMat_Binary;
save('clinical_UniV_Binary','results'), clear results
Locoregional = stringCell_Time(:,1); Distant = stringCell_Time(:,2); Death = stringCell_Time(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',clinicVar);
results.rsMat = rsMat_Time; results.pMat = pMat_Time;
save('clinical_UniV_Time','results'), clear results

cd(startpath)
end