function suppData_GatherResultsRF(pathRF,pathSUPP)

startpath = pwd;

% INITIALIZATION
tables = struct;
% nameOutcomes = {'Locoregional','Distant','Death','Death','DeathSign'};
% fSetNames = {'PET','CT','PETCT','CT','CT'};
nameOutcomes = {'Locoregional','Distant','Death'};
fSetNames = {'PET','CT','CT'};
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
    if strcmp(fields{c},'clinical')
        nExp = 3; % 3 clinical outcomes
        rowNames = cell(nExp,1); 
        for o = 1:nExp
            rowNames{o} = nameOutcomes{o};
        end
    else
        nExp = numel(nameOutcomes);
        rowNames = cell(nExp,1);
        for o = 1:nExp
            rowNames{o} = [nameOutcomes{o},'_',fSetNames{o}];
        end
    end
    cd(paths{c})
    for m = 1:nMetrics
        eval([metrics{m},' = zeros(nExp,1);']);
    end
    for o = 1:nExp
        if strcmp(fields{c},'clinical') || strcmp(fields{c},'volClinical')
            load(['testResultsRF',prefixRead{c},'_',nameOutcomes{o}]) % results gets out of there
        else
            load(['testResultsRF',prefixRead{c},'_',fSetNames{o},'_',nameOutcomes{o}]) % results gets out of there
        end
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