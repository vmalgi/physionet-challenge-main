function classifyResult=ClassifyUSB()
%%oldMaxHeapSize = com.mathworks.services.Prefs.getIntegerPref('JavaMemHeapMax');
%com.mathworks.services.Prefs.setIntegerPref('JavaMemHeapMax', 950);


% Specify Weka path

javaaddpath('weka.jar');
import java.io.File;
%## read file
loader =  weka.core.converters.ArffLoader();
loader.setFile(File('Test.arff'));
TestDataset = loader.getDataSet();
TestDataset.setClassIndex(TestDataset.numAttributes -1);
%% Read model
MyModel = weka.core.SerializationHelper.read ('CM-0180-3I-43F.model');

%% Read the first Instance
NewPCG= TestDataset.instance(0);

%% Predict the class
pred=MyModel.classifyInstance(NewPCG);

% [dis1]=MyModel1.distributionForInstance(NewPCG);

if pred>0
    classifyResult=1;
else
    classifyResult=-1;
end;