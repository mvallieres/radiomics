function testFinalModelsTimeAertsComplete_HN(pathExperiment)

% RADIOMIC SIGNATURE COEFFICIENTS and MEDIAN HAZARD RATIO
coeff = [2.42e-11,-5.38e-03,-1.47e-04,9.39e-06]'; % From original work of (Aerts & Velazquez et al., 2014), courtesy of Hugo Aerts <Hugo_Aerts@dfci.harvard.edu>
medianHR = -0.1191567; % From original work of (Aerts & Velazquez et al., 2014), courtesy of Hugo Aerts <Hugo_Aerts@dfci.harvard.edu>
coeff = -coeff; medianHR = -medianHR; % Due to cox regression implementation differences in this work with original work of (Aerts & Velazquez et al., 2014)
order = numel(coeff);


startpath = pwd;
cd(pathExperiment), load('testing')

% INITIALIZATION
results.model.Name = {'CT_Energy','CT_Compactness','CT_GLN','CT_GLN_HLH'};
results.model.coeff = coeff;
results.model.medianHR = medianHR;
results.model.order = order;

% TESTING THE MODEL
time = testing.Death.timeToEvent;
censoring = 1 - testing.Death.outcome;
resp = zeros(numel(time),1);
data = zeros(numel(time),order);
for n = 1:order
    data(:,n) = testing.Death.sign.CTorig(:,n);
    resp = resp + coeff(n)*data(:,n);
end
results.testData.data = data; results.testData.response = resp; results.testData.outcome = 1 - censoring; results.testData.timeToEvent = time; results.testData.censoring = censoring;
CI = calcCI(resp,time,censoring);
results.CI = CI;

cd(pathExperiment), save(['testResultsCR_','CTorigsignComplete','_','Death'],'results')

cd(startpath)
end