function CI = applyCoxfit(Xtrain,Ytrain,Xtest,Ytest,censTrain,censTest)

coeff = coxphfit(Xtrain,Ytrain,'censoring',censTrain,'baseline',0);
response = responseCox(Xtest,coeff);
CI = calcCI(response,Ytest,censTest);

end