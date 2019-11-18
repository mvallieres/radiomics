function organizeRFexperiments_OldMitaineWithAerts(pathWORK,pathExperimentTextures,pathRF)

% IMPORTANT! CETTE FONCTION DOIT ABSOLUMENT ETRE AUTOMATISEE!!!

startpath = pwd;

cd(pathWORK)
load('clinical')
load('subTypes')

% LOADING VARIABLES
cd(pathExperimentTextures)
load('permTrain'), load('permTest')
trainText = load('training'); trainText = struct2cell(trainText); trainText = trainText{1};
testText = load('testing'); testText = struct2cell(testText); testText = testText{1};
PT_loco = load('testResultsModels_PET_Locoregional'); PT_loco = struct2cell(PT_loco); PT_loco = PT_loco{1};
CT_loco = load('testResultsModels_CT_Locoregional'); CT_loco = struct2cell(CT_loco); CT_loco = CT_loco{1};
PT_distant = load('testResultsModels_PET_Distant'); PT_distant = struct2cell(PT_distant); PT_distant = PT_distant{1};
CT_distant = load('testResultsModels_CT_Distant'); CT_distant = struct2cell(CT_distant); CT_distant = CT_distant{1};
PT_death = load('testResultsModels_PET_Death'); PT_death = struct2cell(PT_death); PT_death = PT_death{1};
CT_death = load('testResultsModels_CT_Death'); CT_death = struct2cell(CT_death); CT_death = CT_death{1};


% ORGANIZING OUTCOMES
training.outcomes.Locoregional = trainText.Locoregional.outcome;
training.outcomes.Distant = trainText.Distant.outcome;
training.outcomes.Death = trainText.Death.outcome;
testing.outcomes.Locoregional = testText.Locoregional.outcome;
testing.outcomes.Distant = testText.Distant.outcome;
testing.outcomes.Death = testText.Death.outcome;
training.timeToEvents.Locoregional = trainText.Locoregional.timeToEvent;
training.timeToEvents.Distant = trainText.Distant.timeToEvent;
training.timeToEvents.Death = trainText.Death.timeToEvent;
testing.timeToEvents.Locoregional = testText.Locoregional.timeToEvent;
testing.timeToEvents.Distant = testText.Distant.timeToEvent;
testing.timeToEvents.Death = testText.Death.timeToEvent;


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



% ORGANIZING TEXTURES (change names if necessary)

% For Locoregional
PET_GLV = PT_loco.trainData.data(:,1); PET_GLN = PT_loco.trainData.data(:,2); PET_LGRE = PT_loco.trainData.data(:,3);
training.textures.Locoregional.PET = table(PET_GLV,PET_GLN,PET_LGRE);
PET_GLV = PT_loco.testData.data(:,1); PET_GLN = PT_loco.testData.data(:,2); PET_LGRE = PT_loco.testData.data(:,3);
testing.textures.Locoregional.PET = table(PET_GLV,PET_GLN,PET_LGRE);

CT_LGZE = CT_loco.trainData.data(:,1); CT_Correlation = CT_loco.trainData.data(:,2); CT_Busyness = CT_loco.trainData.data(:,3);
training.textures.Locoregional.CT = table(CT_LGZE,CT_Correlation,CT_Busyness);
CT_LGZE = CT_loco.testData.data(:,1); CT_Correlation= CT_loco.testData.data(:,2); CT_Busyness = CT_loco.testData.data(:,3);
testing.textures.Locoregional.CT = table(CT_LGZE,CT_Correlation,CT_Busyness);

PET_GLV = PT_loco.trainData.data(:,1); PET_GLN = PT_loco.trainData.data(:,2); PET_LGRE = PT_loco.trainData.data(:,3);
CT_LGZE = CT_loco.trainData.data(:,1); CT_Correlation = CT_loco.trainData.data(:,2); CT_Busyness = CT_loco.trainData.data(:,3);
training.textures.Locoregional.PETCT = table(PET_GLV,PET_GLN,PET_LGRE,CT_LGZE,CT_Correlation,CT_Busyness);
PET_GLV = PT_loco.testData.data(:,1); PET_GLN = PT_loco.testData.data(:,2); PET_LGRE = PT_loco.testData.data(:,3);
CT_LGZE = CT_loco.testData.data(:,1); CT_Correlation= CT_loco.testData.data(:,2); CT_Busyness = CT_loco.testData.data(:,3);
testing.textures.Locoregional.PETCT = table(PET_GLV,PET_GLN,PET_LGRE,CT_LGZE,CT_Correlation,CT_Busyness);


% For Distant
PT_SRHGE = PT_distant.trainData.data(:,1); PT_SUVmax = PT_distant.trainData.data(:,2); PT_ZSV = PT_distant.trainData.data(:,3);
training.textures.Distant.PET = table(PT_SRHGE,PT_SUVmax,PT_ZSV);
PT_SRHGE = PT_distant.testData.data(:,1); PT_SUVmax = PT_distant.testData.data(:,2); PT_ZSV = PT_distant.testData.data(:,3);
testing.textures.Distant.PET = table(PT_SRHGE,PT_SUVmax,PT_ZSV);

CT_LRHGE = CT_distant.trainData.data(:,1); CT_ZSV = CT_distant.trainData.data(:,2); CT_ZSN = CT_distant.trainData.data(:,3);
training.textures.Distant.CT = table(CT_LRHGE,CT_ZSV,CT_ZSN);
CT_LRHGE = CT_distant.testData.data(:,1); CT_ZSV = CT_distant.testData.data(:,2); CT_ZSN = CT_distant.testData.data(:,3);
testing.textures.Distant.CT = table(CT_LRHGE,CT_ZSV,CT_ZSN);

PT_SRHGE = PT_distant.trainData.data(:,1); PT_SUVmax = PT_distant.trainData.data(:,2); PT_ZSV = PT_distant.trainData.data(:,3);
CT_LRHGE = CT_distant.trainData.data(:,1); CT_ZSV = CT_distant.trainData.data(:,2); CT_ZSN = CT_distant.trainData.data(:,3);
training.textures.Distant.PETCT = table(PT_SRHGE,PT_SUVmax,PT_ZSV,CT_LRHGE,CT_ZSV,CT_ZSN);
PT_SRHGE = PT_distant.testData.data(:,1); PT_SUVmax = PT_distant.testData.data(:,2); PT_ZSV = PT_distant.testData.data(:,3);
CT_LRHGE = CT_distant.testData.data(:,1); CT_ZSV = CT_distant.testData.data(:,2); CT_ZSN = CT_distant.testData.data(:,3);
testing.textures.Distant.PETCT = table(PT_SRHGE,PT_SUVmax,PT_ZSV,CT_LRHGE,CT_ZSV,CT_ZSN);


% For Death
PT_SZLGE = PT_death.trainData.data(:,1); PT_LGZE = PT_death.trainData.data(:,2); PT_Contrast = PT_death.trainData.data(:,3);
training.textures.Death.PET = table(PT_SZLGE,PT_LGZE,PT_Contrast);
PT_SZLGE = PT_death.testData.data(:,1); PT_LGZE = PT_death.testData.data(:,2); PT_Contrast = PT_death.testData.data(:,3);
testing.textures.Death.PET = table(PT_SZLGE,PT_LGZE,PT_Contrast);

CT_GLN = CT_death.trainData.data(:,1); CT_SZE = CT_death.trainData.data(:,2); CT_Energy = CT_death.trainData.data(:,3);
training.textures.Death.CT = table(CT_GLN,CT_SZE,CT_Energy);
CT_GLN = CT_death.testData.data(:,1); CT_SZE = CT_death.testData.data(:,2); CT_Energy = CT_death.testData.data(:,3);
testing.textures.Death.CT = table(CT_GLN,CT_SZE,CT_Energy);

PT_SZLGE = PT_death.trainData.data(:,1); PT_LGZE = PT_death.trainData.data(:,2); PT_Contrast = PT_death.trainData.data(:,3);
CT_GLN = CT_death.trainData.data(:,1); CT_SZE = CT_death.trainData.data(:,2); CT_Energy = CT_death.trainData.data(:,3);
training.textures.Death.PETCT = table(PT_SZLGE,PT_LGZE,PT_Contrast,CT_GLN,CT_SZE,CT_Energy);
PT_SZLGE = PT_death.testData.data(:,1); PT_LGZE = PT_death.testData.data(:,2); PT_Contrast = PT_death.testData.data(:,3);
CT_GLN = CT_death.testData.data(:,1); CT_SZE = CT_death.testData.data(:,2); CT_Energy = CT_death.testData.data(:,3);
testing.textures.Death.PETCT = table(PT_SZLGE,PT_LGZE,PT_Contrast,CT_GLN,CT_SZE,CT_Energy);


% For Death (Aerts signature)
PT_death = load('testResultsAerts_PET_Death'); PT_death = struct2cell(PT_death); PT_death = PT_death{1};
CT_death = load('testResultsAerts_CT_Death'); CT_death = struct2cell(CT_death); CT_death = CT_death{1};

PT_Energy = PT_death.trainData.data(:,1); PT_Compactness = PT_death.trainData.data(:,2); PT_GLN = PT_death.trainData.data(:,3); PT_GLN_HLH = PT_death.trainData.data(:,4);
training.textures.DeathSign.PET = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH);
PT_Energy = PT_death.testData.data(:,1); PT_Compactness = PT_death.testData.data(:,2); PT_GLN = PT_death.testData.data(:,3); PT_GLN_HLH = PT_death.testData.data(:,4);
testing.textures.DeathSign.PET = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH);

CT_Energy = CT_death.trainData.data(:,1); CT_Compactness = CT_death.trainData.data(:,2); CT_GLN = CT_death.trainData.data(:,3); CT_GLN_HLH = CT_death.trainData.data(:,4);
training.textures.DeathSign.CT = table(CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);
CT_Energy = CT_death.testData.data(:,1); CT_Compactness = CT_death.testData.data(:,2); CT_GLN = CT_death.testData.data(:,3); CT_GLN_HLH = CT_death.testData.data(:,4);
testing.textures.DeathSign.CT = table(CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);

PT_Energy = PT_death.trainData.data(:,1); PT_Compactness = PT_death.trainData.data(:,2); PT_GLN = PT_death.trainData.data(:,3); PT_GLN_HLH = PT_death.trainData.data(:,4);
CT_Energy = CT_death.trainData.data(:,1); CT_Compactness = CT_death.trainData.data(:,2); CT_GLN = CT_death.trainData.data(:,3); CT_GLN_HLH = CT_death.trainData.data(:,4);
training.textures.DeathSign.PETCT = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH,CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);
PT_Energy = PT_death.testData.data(:,1); PT_Compactness = PT_death.testData.data(:,2); PT_GLN = PT_death.testData.data(:,3); PT_GLN_HLH = PT_death.testData.data(:,4);
CT_Energy = CT_death.testData.data(:,1); CT_Compactness = CT_death.testData.data(:,2); CT_GLN = CT_death.testData.data(:,3); CT_GLN_HLH = CT_death.testData.data(:,4);
testing.textures.DeathSign.PETCT = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH,CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);

% SAVING ORGANIZED DATA
cd(pathRF)
save('training','training'), save('testing','testing')
save('permTrain','permTrain'), save('permTest','permTest')

cd(startpath)
end