function computeAllFinalRF_HN(pathRF,nBoot,seed)

startpath = pwd;

cd(pathRF), load('training')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);

for o = 1:nOutcomes
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = training.outcomes.Death;
    else
        outcome = training.outcomes.(nameOutcomes{o});
    end
    for f = 1:nFset
        text = training.textures.(nameOutcomes{o}).(fSetNames{f}); nText = size(text,2);
        indClinic = training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f});
        cost = training.cost.(nameOutcomes{o}).(fSetNames{f});
        tableTrain = [text,training.clinical.table(:,indClinic)];
        cat = logical([zeros(1,nText),training.clinical.categories(indClinic)]);
        rng(seed), [RF] = trainRF_table(tableTrain,outcome,cat,nBoot,cost);
        RF = compact(RF); % Compact version
        save(['RF_',[fSetNames{f},'clinic'],'_',nameOutcomes{o}],'RF')
    end
end

cd(startpath)
end