function [percentAUCchange,varNames,pVal] = calcFeatureImportanceRF_HN(pathRF,nPerms,nameOutcome,fSetName,seed)

startpath = pwd;

% INITIALIZATION
cd(pathRF), load('training'), load('testing')
if strcmp(nameOutcome,'DeathSign')
    outcome = testing.outcomes.Death;
else
    outcome = testing.outcomes.(nameOutcome);
end
nInst = numel(outcome); % Number of testing patients
results = load(['testResultsRF_',fSetName,'_',nameOutcome]); results = struct2cell(results); results = results{1};
RF = load(['RF_',fSetName,'_',nameOutcome]); RF = struct2cell(RF); RF = RF{1};
if ~strcmp(fSetName,'clinic')
    if ~isempty(strfind(fSetName,'clinic'))
        ind = strfind(fSetName,'clinic');
        fSetName = fSetName(1:ind-1);
    end
    indClinic = training.clinical.bestAdd.(nameOutcome).(fSetName);
    text = testing.textures.(nameOutcome).(fSetName); nText = size(text,2);
else
    load('clinicalPARAM') % parameters gets out of there
    indClinic = parameters.clinical.(nameOutcome);
    text = table;
end


% GETTING PERM AND BOOT SETS (there will be one permutation per bootstrap sample)
perms = zeros(nInst,nPerms); rng(seed)
nBoot = nPerms;
for p = 1:nPerms
    perms(:,p) = datasample(1:nInst,nInst,'Replace',false)';
end
[bootSets,~] = buildBootSet(outcome,nBoot); % No imbalance-adjustment here ('IALR' option)



% PRODUCING RESULTS
tableTest = [text,testing.clinical.table(:,indClinic)];
varNames = tableTest.Properties.VariableNames'; nVar = numel(varNames);
percentAUCchange = zeros(nVar,1);
pVal = zeros(nVar,1);

% Getting all the true bootstrap AUCs
AUC_true = zeros(nPerms,1);
for b = 1:nBoot
    tableBoot = tableTest(bootSets(:,p),:);
    outcomeBoot = outcome(bootSets(:,p));
    [prob] = predictRF(tableBoot,RF);
    [AUC_true(b),~,~,~] = calcPerformMetrics(prob,outcomeBoot,0.5);
end

% Getting all the perm AUCs
for v = 1:nVar
    percentAUCchange_var = 0;
    AUC_perm = zeros(nPerms,1);
    for p = 1:nPerms
        tableBoot = tableTest(bootSets(:,p),:);
        outcomeBoot = outcome(bootSets(:,p));
        tablePerm = tableBoot;
        tablePerm(:,v) = tableBoot(perms(:,p),v);
        [prob] = predictRF(tablePerm,RF);
        [AUC_perm(p),~,~,~] = calcPerformMetrics(prob,outcomeBoot,0.5);
        percentAUCchange_var = percentAUCchange_var + (AUC_perm(p) - AUC_true(p))/AUC_true(p); 
    end
    % if percentAUCdecrease_var > 0, percentAUCdecrease_var = percentAUCdecrease_var * -1; end
    percentAUCchange_var = percentAUCchange_var/nPerms;
    percentAUCchange(v) = percentAUCchange_var;
    pVal(v) = ranksum(AUC_true,AUC_perm,'tail','right');
end
[percentAUCchange,ind] = sort(percentAUCchange,'descend');
varNames = varNames(ind);
pVal = pVal(ind);

cd(startpath)
end