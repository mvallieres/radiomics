function organizeRFexp_Binary(pathWORK,exp,pathTextureResults,pathExperimentsBinary)

startpath = pwd;

cd(pathWORK)
load('clinical')
load('subTypes')

% LOADING VARIABLES
cd(fullfile(pathTextureResults,['Experiment',num2str(exp)]))
load('permTrain'), load('permTest')
trainText = load('training'); trainText = struct2cell(trainText); trainText = trainText{1};
testText = load('testing'); testText = struct2cell(testText); testText = testText{1};


% ORGANIZING DATA
cd(pathExperimentsBinary), mkdir(['Experiment',num2str(exp)]), cd(['Experiment',num2str(exp)]), pathExperiment = pwd;
age_train = [clinical.HGJ.Age;clinical.CHUS.Age]; age_test = [clinical.HMR.Age;clinical.CHUM.Age];
subType_train = [subTypes.HGJ;subTypes.CHUS]; subType_test = [subTypes.HMR;subTypes.CHUM];
tStage_train = [clinical.HGJ.T_stage;clinical.CHUS.T_stage]; tStage_test = [clinical.HMR.T_stage;clinical.CHUM.T_stage];
nStage_train = [clinical.HGJ.N_stage;clinical.CHUS.N_stage]; nStage_test = [clinical.HMR.N_stage;clinical.CHUM.N_stage];
age_train = age_train(permTrain); age_test = age_test(permTest);
subType_train = subType_train(permTrain); subType_test = subType_test(permTest);
tStage_train = tStage_train(permTrain); tStage_test = tStage_test(permTest);
nStage_train = nStage_train(permTrain); nStage_test = nStage_test(permTest);

% Locoregional
cd(fullfile(pathTextureResults,['Experiment',num2str(exp)]))
load('testResults_CT_Locoregional'), resultsCT = results;
load('testResults_PET_Locoregional'), resultsPT = results;
training.Locoregional.variables.PET.Data = {resultsPT.trainData.data(:,1),resultsPT.trainData.data(:,2),resultsPT.trainData.data(:,3),...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Locoregional.variables.PET.Categories = logical([0,0,0,0,1,1,1]);
training.Locoregional.variables.PET.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Locoregional.variables.PET.Data = {resultsPT.testData.data(:,1),resultsPT.testData.data(:,2),resultsPT.testData.data(:,3),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Locoregional.variables.PET.Categories = logical([0,0,0,0,1,1,1]);
testing.Locoregional.variables.PET.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Locoregional.variables.CT.Data = {resultsCT.trainData.data(:,1),resultsCT.trainData.data(:,2),resultsCT.trainData.data(:,3),resultsCT.trainData.data(:,4)...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Locoregional.variables.CT.Categories = logical([0,0,0,0,0,1,1,1]);
training.Locoregional.variables.CT.Info = {resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.nonText;resultsCT.model.Feature4.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Locoregional.variables.CT.Data = {resultsCT.testData.data(:,1),resultsCT.testData.data(:,2),resultsCT.testData.data(:,3),resultsCT.testData.data(:,4),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Locoregional.variables.CT.Categories = logical([0,0,0,0,0,1,1,1]);
testing.Locoregional.variables.CT.Info = {resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.nonText;resultsCT.model.Feature4.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Locoregional.variables.PET_CT.Data = {resultsPT.trainData.data(:,1),resultsPT.trainData.data(:,2),resultsPT.trainData.data(:,3),resultsCT.trainData.data(:,1),resultsCT.trainData.data(:,2),resultsCT.trainData.data(:,3),resultsCT.trainData.data(:,4),...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Locoregional.variables.PET_CT.Categories = logical([0,0,0,0,0,0,0,0,1,1,1]);
training.Locoregional.variables.PET_CT.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.FullName;resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.nonText;resultsCT.model.Feature4.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Locoregional.variables.PET_CT.Data = {resultsPT.testData.data(:,1),resultsPT.testData.data(:,2),resultsPT.testData.data(:,3),resultsCT.testData.data(:,1),resultsCT.testData.data(:,2),resultsCT.testData.data(:,3),resultsCT.testData.data(:,4),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Locoregional.variables.PET_CT.Categories = logical([0,0,0,0,0,0,0,0,1,1,1]);
testing.Locoregional.variables.PET_CT.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.FullName;resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.nonText;resultsCT.model.Feature4.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Locoregional.outcome = trainText.Locoregional.outcome;
training.Locoregional.timeToEvent = trainText.Locoregional.timeToEvent;
testing.Locoregional.outcome = testText.Locoregional.outcome;
testing.Locoregional.timeToEvent = testText.Locoregional.timeToEvent;

% Distant
cd(fullfile(pathTextureResults,['Experiment',num2str(exp)]))
load('testResults_CT_Distant'), resultsCT = results;
load('testResults_PET_Distant'), resultsPT = results;
training.Distant.variables.PET.Data = {resultsPT.trainData.data(:,1),resultsPT.trainData.data(:,2),resultsPT.trainData.data(:,3),...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Distant.variables.PET.Categories = logical([0,0,0,0,1,1,1]);
training.Distant.variables.PET.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.nonText;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Distant.variables.PET.Data = {resultsPT.testData.data(:,1),resultsPT.testData.data(:,2),resultsPT.testData.data(:,3),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Distant.variables.PET.Categories = logical([0,0,0,0,1,1,1]);
testing.Distant.variables.PET.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.nonText;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Distant.variables.CT.Data = {resultsCT.trainData.data(:,1),resultsCT.trainData.data(:,2),resultsCT.trainData.data(:,3),...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Distant.variables.CT.Categories = logical([0,0,0,0,1,1,1]);
training.Distant.variables.CT.Info = {resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Distant.variables.CT.Data = {resultsCT.testData.data(:,1),resultsCT.testData.data(:,2),resultsCT.testData.data(:,3),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Distant.variables.CT.Categories = logical([0,0,0,0,1,1,1]);
testing.Distant.variables.CT.Info = {resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Distant.variables.PET_CT.Data = {resultsPT.trainData.data(:,1),resultsPT.trainData.data(:,2),resultsPT.trainData.data(:,3),resultsCT.trainData.data(:,1),resultsCT.trainData.data(:,2),resultsCT.trainData.data(:,3),...
                                            age_train,subType_train,tStage_train,nStage_train};
training.Distant.variables.PET_CT.Categories = logical([0,0,0,0,0,0,0,1,1,1]);
training.Distant.variables.PET_CT.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.nonText;resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
testing.Distant.variables.PET_CT.Data = {resultsPT.testData.data(:,1),resultsPT.testData.data(:,2),resultsPT.testData.data(:,3),resultsCT.testData.data(:,1),resultsCT.testData.data(:,2),resultsCT.testData.data(:,3),...
                                            age_test,subType_test,tStage_test,nStage_test};
testing.Distant.variables.PET_CT.Categories = logical([0,0,0,0,0,0,0,1,1,1]);
testing.Distant.variables.PET_CT.Info = {resultsPT.model.Feature1.nonText;resultsPT.model.Feature2.FullName;resultsPT.model.Feature3.nonText;resultsCT.model.Feature1.FullName;resultsCT.model.Feature2.FullName;resultsCT.model.Feature3.FullName;...
                                            'Age';'SubType';'T_Stage';'N_Stage'};
training.Distant.outcome = trainText.Distant.outcome;
training.Distant.timeToEvent = trainText.Distant.timeToEvent;
testing.Distant.outcome = testText.Distant.outcome;
testing.Distant.timeToEvent = testText.Distant.timeToEvent;

cd(pathExperiment)
save('permTrain','permTrain'), save('permTest','permTest')
save('training','training'), save('testing','testing')


% CREATING FEAURE SETS
mkdir('FSET'), cd('FSET')
fSetNames = {'PET','CT','PET_CT'};
nameOutcomes = {'Locoregional','Distant'};
for f = 1:numel(fSetNames)
    for o = 1:numel(nameOutcomes)
        fSet.Data = training.(nameOutcomes{o}).variables.(fSetNames{f}).Data;
        fSet.Categories = training.(nameOutcomes{o}).variables.(fSetNames{f}).Categories;
        fSet.Info = training.(nameOutcomes{o}).variables.(fSetNames{f}).Info;
        save(['FSET_',fSetNames{f},'_',nameOutcomes{o}],'fSet')
        clear fSet
    end
end

cd(startpath)
end