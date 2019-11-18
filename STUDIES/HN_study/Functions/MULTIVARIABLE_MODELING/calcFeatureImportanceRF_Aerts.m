function calcFeatureImportanceRF_Aerts(pathRF,nPerms,seed)

startpath = pwd;
cd(pathRF), load('testingVariableImportance') % 'variableImportance' gets out of there

[percentAUCdecrease,varNames] = calcFeatureImportanceRF_HN(pathRF,nPerms,'DeathSign','CT',seed);
variableImportance.DeathSign.CT.percentAUCdecrease = percentAUCdecrease;
variableImportance.DeathSign.CT.varNames = varNames;
cd(pathRF), save('testingVariableImportance','variableImportance')

cd(startpath)
end