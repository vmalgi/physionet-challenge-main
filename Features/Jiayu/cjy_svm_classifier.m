function predict = cjy_svm_classifier(features,model)
 predict = svmclassify(model,features);
 if predict == 0
     predict = -1;
 end

