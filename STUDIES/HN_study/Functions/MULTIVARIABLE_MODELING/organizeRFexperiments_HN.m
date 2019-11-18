function organizeRFexperiments_HN(pathWORK,pathExperimentTextures,pathRF,nameOutcomes,fSetNames)

startpath = pwd;

% INITIALIZATION
nOutcomes = numel(nameOutcomes);
nFset = numel(fSetNames);
cd(pathWORK), load('clinical'), load('subTypes')
cd(pathExperimentTextures), load('permTrain'), load('permTest')
trainText = load('training'); trainText = struct2cell(trainText); trainText = trainText{1};
testText = load('testing'); testText = struct2cell(testText); testText = testText{1};
training = struct; testing = struct;


% ORGANIZING CLINICAL DATA
age_train = [clinical.HGJ.Age;clinical.CHUS.Age]; age_test = [clinical.HMR.Age;clinical.CHUM.Age];
subType_train = [subTypes.HGJ;subTypes.CHUS]; subType_test = [subTypes.HMR;subTypes.CHUM];
tStage_train = [clinical.HGJ.T_stage;clinical.CHUS.T_stage]; tStage_test = [clinical.HMR.T_stage;clinical.CHUM.T_stage];
nStage_train = [clinical.HGJ.N_stage;clinical.CHUS.N_stage]; nStage_test = [clinical.HMR.N_stage;clinical.CHUM.N_stage];
tnmStage_train = [clinical.HGJ.TNM_stage;clinical.CHUS.TNM_stage]; tnmStage_test = [clinical.HMR.TNM_stage;clinical.CHUM.TNM_stage];
age_train = age_train(permTrain); age_test = age_test(permTest);
subType_train = subType_train(permTrain); subType_test = subType_test(permTest);
tStage_train = tStage_train(permTrain); tStage_test = tStage_test(permTest);
nStage_train = nStage_train(permTrain); nStage_test = nStage_test(permTest);
tnmStage_train = tnmStage_train(permTrain); tnmStage_test = tnmStage_test(permTest); 
Age = age_train; SubType = subType_train; T_Stage = tStage_train; N_Stage = nStage_train; TNM_Stage = tnmStage_train;
training.clinical.table = table(Age,SubType,T_Stage,N_Stage,TNM_Stage);
training.clinical.categories = [0,1,1,1,1];
Age = age_test; SubType = subType_test; T_Stage = tStage_test; N_Stage = nStage_test; TNM_Stage = tnmStage_test;
testing.clinical.table = table(Age,SubType,T_Stage,N_Stage,TNM_Stage);
testing.clinical.categories = [0,1,1,1,1];


% ORGANIZING TEXTURES: PET, CT and PETCT (sep) feature sets
cd(pathExperimentTextures)
for o = 1:nOutcomes
    training.outcomes.(nameOutcomes{o}) = trainText.(nameOutcomes{o}).outcome;
    training.timeToEvents.(nameOutcomes{o}) = trainText.(nameOutcomes{o}).timeToEvent;
    testing.outcomes.(nameOutcomes{o}) = testText.(nameOutcomes{o}).outcome;
    testing.timeToEvents.(nameOutcomes{o}) = testText.(nameOutcomes{o}).timeToEvent;
    for f = 1:nFset
        results = load(['testResultsLR_',fSetNames{f},'_',nameOutcomes{o}]); results = struct2cell(results); results = results{1};
        order = results.model.order;
        
        % Training
        strVar = [];
        for r = 1:order
            if results.model.(['Feature',num2str(r)]).isTexture
                name = [results.model.(['Feature',num2str(r)]).scan,'_',results.model.(['Feature',num2str(r)]).textName];
            else
                name = results.model.(['Feature',num2str(r)]).nonText;
            end
            trainData = results.trainData.data(:,r);
            eval([name,' = trainData;']);
            strVar = [strVar,name,','];
        end
        strVar(end) = [];
        eval(['training.textures.(nameOutcomes{o}).(fSetNames{f}) = table(',strVar,');']);
        
        % Testing
        strVar = [];
        for r = 1:order
            if results.model.(['Feature',num2str(r)]).isTexture
                name = [results.model.(['Feature',num2str(r)]).scan,'_',results.model.(['Feature',num2str(r)]).textName];
            else
                name = results.model.(['Feature',num2str(r)]).nonText;
            end
            testData = results.testData.data(:,r);
            eval([name,' = testData;']);
            strVar = [strVar,name,','];
        end
        strVar(end) = [];
        eval(['testing.textures.(nameOutcomes{o}).(fSetNames{f}) = table(',strVar,');']);
    end
end


% SAVING ORGANIZED DATA
cd(pathRF)
save('training','training'), save('testing','testing')
save('permTrain','permTrain'), save('permTest','permTest')

cd(startpath)
end