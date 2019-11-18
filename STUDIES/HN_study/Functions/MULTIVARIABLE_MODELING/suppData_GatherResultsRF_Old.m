function suppData_GatherResultsRF_Old(pathRF,pathSUPP)

startpath = pwd;

% INITIALIZATION
tables = struct;
nameOutcomes = {'Locoregional','Distant','Death','Death','DeathSign'}; nExp = numel(nameOutcomes);
fSetNames = {'PET','CT','PETCT','CT','CT'};
rowNames = cell(nExp,1);
for o = 1:nExp
    rowNames{o} = [nameOutcomes{o},'_',fSetNames{o}];
end
metrics = {'AUC','Sensitivity','Specificity','Accuracy','CI'}; nMetrics = numel(metrics);
strMetrics = [];
for m = 1:nMetrics
    strMetrics = [strMetrics,metrics{m},','];
end
strMetrics(end) = [];
paths = {pathSUPP,pathSUPP,pathRF}; nComp = numel(paths);
prefixRead = {'radiomics','clinical',''};
fields = {'radiomics','clinical','combined'};


% PROCESSING 
for c = 1:nComp
    cd(paths{c})
    for m = 1:nMetrics
        eval([metrics{m},' = zeros(nExp,1);']);
    end
    for o = 1:nExp
        load(['testResultsRF',prefixRead{c},'_',fSetNames{o},'_',nameOutcomes{o}]) % results gets out of there
        for m = 1:nMetrics
            eval([metrics{m},'(o,1) = results.',metrics{m},';']);
        end
    end
    eval(['tableExp = table(',strMetrics,',''RowNames'',rowNames);']);
    tables.(fields{c}) = tableExp;
end

cd(pathSUPP), save('summaryTables','tables')

cd(startpath)
end