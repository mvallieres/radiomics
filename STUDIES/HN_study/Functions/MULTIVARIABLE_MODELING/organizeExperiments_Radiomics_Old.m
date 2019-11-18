function organizeExperiments_Radiomics_Old(pathFeatures,pathLR,pathCR,outcomes,timeToEvent,nPatient,nonTextName,fSetNames,paramSEP,baselineSEP,freedomSEP,textType,textName)

startpath = pwd;

permTrain = randperm(nPatient.HGJ + nPatient.CHUS); permTest = randperm(nPatient.HMR + nPatient.CHUM); % Random permutation of instances.

nNonText = numel(nonTextName);
nFset = numel(fSetNames);
nameOutcomes = fieldnames(outcomes.HGJ); nOutcomes = numel(nameOutcomes);
nTextType = numel(textType);

training.nonText = struct; training.text = struct; training.outcomes = struct; training.timeToEvent = struct; 
testing.nonText = struct; testing.text = struct; testing.outcomes = struct; testing.timeToEvent = struct;

% Organizing non-texture features
cd(pathFeatures)
train1 = load('nonText_HGJ_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
train2 = load('nonText_CHUS_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
test1 = load('nonText_HMR_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
test2 = load('nonText_CHUM_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
for i = 1:nNonText
    training.nonText.(nonTextName{i}).Data = [train1.(nonTextName{i}).Data;train2.(nonTextName{i}).Data];
    testing.nonText.(nonTextName{i}).Data = [test1.(nonTextName{i}).Data;test2.(nonTextName{i}).Data];
    training.nonText.(nonTextName{i}).Data = training.nonText.(nonTextName{i}).Data(permTrain);
    testing.nonText.(nonTextName{i}).Data = testing.nonText.(nonTextName{i}).Data(permTest);
end

% Organizing texture features
for f = 1:nFset
    training.text.(fSetNames{f}) = struct; testing.text.(fSetNames{f}) = struct;
    if strcmp(fSetNames{f},'PET') % Using only PET
        train1 = load('text_HGJ_PT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
        train2 = load('text_CHUS_PT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
        test1 = load('text_HMR_PT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
        test2 = load('text_CHUM_PT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
        train1 = {train1}; train2 = {train2};
        test1 = {test1}; test2 = {test2};
        training.text.(fSetNames{f}).param = paramSEP; 
        training.text.(fSetNames{f}).baseline = baselineSEP; 
        training.text.(fSetNames{f}).freedom = freedomSEP;
        training.text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
        textCellName = {'PT'};
    elseif strcmp(fSetNames{f},'CT') % Using only CT
        train1 = load('text_HGJ_CT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
        train2 = load('text_CHUS_CT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
        test1 = load('text_HMR_CT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
        test2 = load('text_CHUM_CT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
        train1 = {train1}; train2 = {train2};
        test1 = {test1}; test2 = {test2};
        training.text.(fSetNames{f}).param = paramSEP; 
        training.text.(fSetNames{f}).baseline = baselineSEP; 
        training.text.(fSetNames{f}).freedom = freedomSEP;
        training.text.(fSetNames{f}).paramName = {'Scale','Quant.algo','Ng'};
        textCellName = {'CT'};
    end
    
    for i = 1:numel(train1)
        sizeParam = size(train1{1});
        training.text.(fSetNames{f}).(textCellName{i}) = cell(sizeParam);
        testing.text.(fSetNames{f}).(textCellName{i}) = cell(sizeParam);
        for n = 1:numel(train1{1})
            for type = 1:nTextType
                for text = 1:numel(textName{type})
                    training.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = [train1{i}{n}.(textType{type}).(textName{type}{text}).Data;train2{i}{n}.(textType{type}).(textName{type}{text}).Data];
                    testing.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = [test1{i}{n}.(textType{type}).(textName{type}{text}).Data;test2{i}{n}.(textType{type}).(textName{type}{text}).Data];
                    training.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = training.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data(permTrain);
                    testing.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data = testing.text.(fSetNames{f}).(textCellName{i}){n}.(textType{type}).(textName{type}{text}).Data(permTest);
                end
            end
        end
    end
end


% Organizing outcomes
for o = 1:nOutcomes
    training.outcomes.(nameOutcomes{o}) = [outcomes.HGJ.(nameOutcomes{o});outcomes.CHUS.(nameOutcomes{o})];
    testing.outcomes.(nameOutcomes{o}) = [outcomes.HMR.(nameOutcomes{o});outcomes.CHUM.(nameOutcomes{o})];
    training.timeToEvent.(nameOutcomes{o}) = [timeToEvent.HGJ.(nameOutcomes{o});timeToEvent.CHUS.(nameOutcomes{o})];
    testing.timeToEvent.(nameOutcomes{o}) = [timeToEvent.HMR.(nameOutcomes{o});timeToEvent.CHUM.(nameOutcomes{o})];
    training.outcomes.(nameOutcomes{o}) = training.outcomes.(nameOutcomes{o})(permTrain);
    testing.outcomes.(nameOutcomes{o}) = testing.outcomes.(nameOutcomes{o})(permTest);
    training.timeToEvent.(nameOutcomes{o}) = training.timeToEvent.(nameOutcomes{o})(permTrain);
    testing.timeToEvent.(nameOutcomes{o}) = testing.timeToEvent.(nameOutcomes{o})(permTest);
end

% Re-arranging by outcomes (same features for all outcomes here)
tempTraining = training; tempTesting = testing; clear training testing
for o = 1:nOutcomes
    training.(nameOutcomes{o}).nonText = tempTraining.nonText;
    training.(nameOutcomes{o}).text = tempTraining.text;
    training.(nameOutcomes{o}).outcome = tempTraining.outcomes.(nameOutcomes{o});
    training.(nameOutcomes{o}).timeToEvent = tempTraining.timeToEvent.(nameOutcomes{o});
    testing.(nameOutcomes{o}).nonText = tempTesting.nonText;
    testing.(nameOutcomes{o}).text = tempTesting.text;
    testing.(nameOutcomes{o}).outcome = tempTesting.outcomes.(nameOutcomes{o});
    testing.(nameOutcomes{o}).timeToEvent = tempTesting.timeToEvent.(nameOutcomes{o});
end

cd(pathLR), save('training','training'), save('testing','testing'), save('permTrain','permTrain'), save('permTest','permTest')
cd(pathCR), save('training','training'), save('testing','testing'), save('permTrain','permTrain'), save('permTest','permTest')

cd(startpath)
end