function organizeExperiments_Radiomics(pathFeatures,pathLR,pathCR,outcomes,timeToEvent,nPatient,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)

startpath = pwd;

% INITIALIZATION
roi = 'GTVtot'; % Using "GTVtot" ROI for all experiments
permTrain = randperm(nPatient.HGJ + nPatient.CHUS); permTest = randperm(nPatient.HMR + nPatient.CHUM); % Random permutation of instances.
nNonText = numel(nonTextName);
nFset = numel(fSetNames);
nameOutcomes = fieldnames(outcomes.HGJ); nOutcomes = numel(nameOutcomes);
nTextType = numel(textType);


% DATA ORGANIZATION
training = struct; testing = struct;
for o = 1:nOutcomes
    
    % Organizing non-texture features
    cd(pathFeatures)
    train1 = load(['nonText_HGJ_',roi]); train1 = struct2cell(train1); train1 = train1{1};
    train2 = load(['nonText_CHUS_',roi]); train2 = struct2cell(train2); train2 = train2{1};
    test1 = load(['nonText_HMR_',roi]); test1 = struct2cell(test1); test1 = test1{1};
    test2 = load(['nonText_CHUM_',roi]); test2 = struct2cell(test2); test2 = test2{1};
    for i = 1:nNonText
        training.(nameOutcomes{o}).nonText.(nonTextName{i}).Data = [train1.(nonTextName{i}).Data;train2.(nonTextName{i}).Data];
        testing.(nameOutcomes{o}).nonText.(nonTextName{i}).Data = [test1.(nonTextName{i}).Data;test2.(nonTextName{i}).Data];
        training.(nameOutcomes{o}).nonText.(nonTextName{i}).Data = training.(nameOutcomes{o}).nonText.(nonTextName{i}).Data(permTrain);
        testing.(nameOutcomes{o}).nonText.(nonTextName{i}).Data = testing.(nameOutcomes{o}).nonText.(nonTextName{i}).Data(permTest);
    end

    % Organizing texture features
    for f = 1:nFset
        if strcmp(fSetNames{f},'PET') % Using only PET
            train1 = load(['text_HGJ_PT_',roi]); train1 = struct2cell(train1); train1 = train1{1};
            train2 = load(['text_CHUS_PT_',roi]); train2 = struct2cell(train2); train2 = train2{1};
            test1 = load(['text_HMR_PT_',roi]); test1 = struct2cell(test1); test1 = test1{1};
            test2 = load(['text_CHUM_PT_',roi]); test2 = struct2cell(test2); test2 = test2{1};
            train1 = {train1}; train2 = {train2};
            test1 = {test1}; test2 = {test2};
            training.(nameOutcomes{o}).text.(fSetNames{f}).param = paramSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).baseline = baselineSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).freedom = freedomSEP;
            training.(nameOutcomes{o}).text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
            textCellName = {'PT'};
        elseif strcmp(fSetNames{f},'CT') % Using only CT
            train1 = load(['text_HGJ_CT_',roi]); train1 = struct2cell(train1); train1 = train1{1};
            train2 = load(['text_CHUS_CT_',roi]); train2 = struct2cell(train2); train2 = train2{1};
            test1 = load(['text_HMR_CT_',roi]); test1 = struct2cell(test1); test1 = test1{1};
            test2 = load(['text_CHUM_CT_',roi]); test2 = struct2cell(test2); test2 = test2{1};
            train1 = {train1}; train2 = {train2};
            test1 = {test1}; test2 = {test2};
            training.(nameOutcomes{o}).text.(fSetNames{f}).param = paramSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).baseline = baselineSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).freedom = freedomSEP;
            training.(nameOutcomes{o}).text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
            textCellName = {'CT'};
        elseif strcmp(fSetNames{f},'PETCT') % Using both PET and CT
            train1_1 = load(['text_HGJ_PT_',roi]); train1_1 = struct2cell(train1_1); train1_1 = train1_1{1};
            train1_2 = load(['text_HGJ_CT_',roi]); train1_2 = struct2cell(train1_2); train1_2 = train1_2{1};
            train2_1 = load(['text_CHUS_PT_',roi]); train2_1 = struct2cell(train2_1); train2_1 = train2_1{1};
            train2_2 = load(['text_CHUS_CT_',roi]); train2_2 = struct2cell(train2_2); train2_2 = train2_2{1};
            test1_1 = load(['text_HMR_PT_',roi]); test1_1 = struct2cell(test1_1); test1_1 = test1_1{1};
            test1_2 = load(['text_HMR_CT_',roi]); test1_2 = struct2cell(test1_2); test1_2 = test1_2{1};
            test2_1 = load(['text_CHUM_PT_',roi]); test2_1 = struct2cell(test2_1); test2_1 = test2_1{1};
            test2_2 = load(['text_CHUM_CT_',roi]); test2_2 = struct2cell(test2_2); test2_2 = test2_2{1};
            train1 = {train1_1,train1_2}; train2 = {train2_1,train2_2};
            test1 = {test1_1,test1_2}; test2 = {test2_1,test2_2};
            training.(nameOutcomes{o}).text.(fSetNames{f}).param = paramSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).baseline = baselineSEP; 
            training.(nameOutcomes{o}).text.(fSetNames{f}).freedom = freedomSEP;
            training.(nameOutcomes{o}).text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
            textCellName = {'PT','CT'};
        end

        for i = 1:numel(train1)
            sizeParam = size(train1{1});
            training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}) = cell(sizeParam);
            testing.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}) = cell(sizeParam);
            for n = 1:numel(train1{1})
                for type = 1:nTextType
                    for text = 1:numel(textName{type})
                        training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = [train1{i}{n}.(textType{type}).(textName{type}{text}).Data;train2{i}{n}.(textType{type}).(textName{type}{text}).Data];
                        testing.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = [test1{i}{n}.(textType{type}).(textName{type}{text}).Data;test2{i}{n}.(textType{type}).(textName{type}{text}).Data];
                        training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = training.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data(permTrain);
                        testing.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = testing.(nameOutcomes{o}).text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data(permTest);
                    end
                end
            end
        end
    end
end

% Organizing outcomes
for o = 1:nOutcomes
    training.(nameOutcomes{o}).outcome = [outcomes.HGJ.(nameOutcomes{o});outcomes.CHUS.(nameOutcomes{o})];
    testing.(nameOutcomes{o}).outcome = [outcomes.HMR.(nameOutcomes{o});outcomes.CHUM.(nameOutcomes{o})];
    training.(nameOutcomes{o}).timeToEvent = [timeToEvent.HGJ.(nameOutcomes{o});timeToEvent.CHUS.(nameOutcomes{o})];
    testing.(nameOutcomes{o}).timeToEvent = [timeToEvent.HMR.(nameOutcomes{o});timeToEvent.CHUM.(nameOutcomes{o})];
    training.(nameOutcomes{o}).outcome = training.(nameOutcomes{o}).outcome(permTrain);
    testing.(nameOutcomes{o}).outcome = testing.(nameOutcomes{o}).outcome(permTest);
    training.(nameOutcomes{o}).timeToEvent = training.(nameOutcomes{o}).timeToEvent(permTrain);
    testing.(nameOutcomes{o}).timeToEvent = testing.(nameOutcomes{o}).timeToEvent(permTest);
end


cd(pathLR), save('training','training'), save('testing','testing'), save('permTrain','permTrain'), save('permTest','permTest')
cd(pathCR), save('training','training'), save('testing','testing'), save('permTrain','permTrain'), save('permTest','permTest')

cd(startpath)
end