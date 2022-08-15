classdef Network
    properties
        sizes
        nLayers
        biases
        weights
        velocities
        layerTypes
        costFunction
    end
    methods
        function obj = initialize_network(obj,varargin)
            networkSizes = varargin{1};
            obj.sizes = networkSizes;
            obj.nLayers = length(networkSizes);
            obj.biases = cell(1,obj.nLayers-1);
            obj.weights = cell(1,obj.nLayers-1);
            obj.velocities = cell(1,obj.nLayers-1);
            if nargin == 2
                obj.layerTypes = repmat({'sigmoid'},[1, obj.nLayers-1]);
                obj.costFunction = 'crossEntropy';
            elseif nargin == 3
                obj.layerTypes = varargin{2};
                if obj.nLayers ~= length(obj.layerTypes) + 1
                    error('Wrong number of layers described')
                end
                obj.costFunction = 'crossEntropy';
            elseif nargin == 4
                obj.layerTypes = varargin{2};
                if obj.nLayers ~= length(obj.layerTypes) + 1
                    error('Wrong number of layers described')
                end
                obj.costFunction = varargin{3};
            else
                error('Calling Network with wrong number of input arguments')
            end
            if obj.sizes(end) == 1 && strcmpi(obj.layerTypes{end},'softmax')
               error('Cant use softmax if only 1 output') 
            end
        end
        
        function obj = initialize_weights_old(obj)
            for i = 1:obj.nLayers-1
                obj.biases{i} = normrnd(0,1,obj.sizes(i+1),1);
                obj.weights{i} = normrnd(0,1,obj.sizes(i+1),obj.sizes(i));
                obj.velocities{i} = zeros(size(obj.weights{i}));
            end
        end
        
        function obj = initialize_weights_new(obj)
            for i = 1:obj.nLayers-1
                obj.biases{i} = normrnd(0,1,obj.sizes(i+1),1);
                obj.weights{i} = (normrnd(0,1,obj.sizes(i+1),obj.sizes(i)))./(sqrt(obj.sizes(i)));
                obj.velocities{i} = zeros(size(obj.weights{i}));
            end
            
        end
        
        function obj = initialize_weights_He(obj)
            for i = 1:obj.nLayers-1
                obj.biases{i} = zeros(obj.sizes(i+1),1);
                obj.weights{i} = (normrnd(0,1,obj.sizes(i+1),obj.sizes(i))).*(2./(sqrt(obj.sizes(i))));
                obj.velocities{i} = zeros(size(obj.weights{i}));
            end
        end
        
        function obj = initialize_weights_Glorot(obj)
            for i = 1:obj.nLayers-1
                obj.biases{i} = normrnd(0,1,obj.sizes(i+1),1);
                obj.weights{i} = -sqrt(6./(obj.sizes(i)+obj.sizes(i+1))) + (2.*sqrt(6./(obj.sizes(i)+obj.sizes(i+1)))).*rand(obj.sizes(i+1),obj.sizes(i));
                obj.velocities{i} = zeros(size(obj.weights{i}));
            end
        end
        function a = feed_forward(obj,a)
            for i = 1:(obj.nLayers-1)
                switch obj.layerTypes{i}
                    case 'sigmoid'
                        a = sigmoid(obj.weights{i}*a + obj.biases{i});
                    case 'tanh'
                        a = tanh(obj.weights{i}*a + obj.biases{i});
                    case 'rectifiedLinear'
                        a = rectLinear(obj.weights{i}*a + obj.biases{i});
                    case 'leakyRectifiedLinear'
                        a = leakyRectLinear(obj.weights{i}*a + obj.biases{i});
                    case 'softmax'
                        a = softmax(obj.weights{i}*a + obj.biases{i});
                    case 'sigmoidPlusLinear'
                        a = sigmoidPlusLinear(obj.weights{i}*a + obj.biases{i});
                    otherwise
                        error('Layer Types are wrong')
                end
            end
        end
        
        function obj = stochastic_gradient_descent(obj,trainingFeatures,trainingOutcomes,epochs,miniBatchSize,learningRate,lmdba,mu,dropoutMethod,percentToDrop)
            if length(percentToDrop) ~= length(obj.biases)
                error('Percentage of nodes/weights to dropout/dropconnect should be equal to number of hidden layers + 1 (output layer)')
            end
            if length(learningRate) ~= length(obj.biases)
                if length(learningRate) == 1
                    learningRate = repmat(learningRate,size(obj.biases));
                else
                    error('Learning Rate vector should be equal to number of hidden layers + 1 (output layer)')
                end
            end
            for i = 1:epochs
                p = randperm(size(trainingFeatures,2));
                trainingFeatures = trainingFeatures(:,p);
                trainingOutcomes = trainingOutcomes(:,p);
                miniBatchFeatures = cell(size(trainingFeatures,1),floor(size(trainingFeatures,2)./miniBatchSize));
                miniBatchOutcomes = cell(size(trainingOutcomes,1),floor(size(trainingOutcomes,2)./miniBatchSize));
                
                for j  = 1:size(miniBatchFeatures,2)
                    if strcmpi(dropoutMethod,'Dropout')                                    
                        removedWeights = cell(size(obj.weights));
                        indicesToRemoveCell = cell(size(obj.weights));
                        for k = 1:length(obj.weights)
                            removedWeights{k} = zeros(size(obj.weights{k}));
                            totalNumberOfNodes = obj.sizes(k);
                            indicesToRemove = randsample([1:totalNumberOfNodes],ceil(totalNumberOfNodes.*((percentToDrop(k))/100)));
                            indicesToRemoveCell{k} = indicesToRemove;
                            removedWeights{k}(:,indicesToRemove) = obj.weights{k}(:,indicesToRemove);
                            if k~=1
                                removedWeights{k-1}(indicesToRemove,:) = obj.weights{k-1}(indicesToRemove,:);
                            end
                        end
                        for k = 1:length(obj.weights)
                            obj.weights{k}(:,indicesToRemoveCell{k}) = 0;
                            if k~=1
                                obj.weights{k-1}(indicesToRemoveCell{k},:) = 0;
                            end
                            
                        end
                        
                        
%                         for k = 1:length(obj.weights)
%                             figure('Name',strcat(j,k,'Weights'))
%                             temp = obj.weights{k};
%                             temp(temp~=0) = 1;
%                             imagesc(temp);
%                             colormap gray
%                             colorbar
%                             
%                             
%                             figure('Name',strcat(j,k,'Removed Weights'))
%                             temp = removedWeights{k};
%                             temp(temp~=0) = 1;
%                             imagesc(temp);
%                             colormap gray
%                             colorbar
%                         end
                        
                    elseif strcmpi(dropoutMethod,'Dropconnect')
                        
                        removedWeights = cell(size(obj.weights));
                        indicesToRemoveCell = cell(size(obj.weights));
                        for k = 1:length(obj.weights)
                            removedWeights{k} = zeros(size(obj.weights{k}));
                            totalNumberOfWeights = (size(obj.weights{k},1)*size(obj.weights{k},2));
                            indicesToRemove = randsample([1:totalNumberOfWeights],ceil(totalNumberOfWeights.*((percentToDrop(k))/100)));
                            indicesToRemoveCell{k} = indicesToRemove;
                            removedWeights{k}(indicesToRemove) = obj.weights{k}(indicesToRemove);
                            obj.weights{k}(indicesToRemove) = 0;
                        end
                    end
                    
                    
                    miniBatchFeatures{j} = trainingFeatures(:,(1+((j-1).*miniBatchSize):(j.*miniBatchSize)));
                    miniBatchOutcomes{j} = trainingOutcomes(:,(1+((j-1).*miniBatchSize):(j.*miniBatchSize)));
                    obj = updateWeightsAndBiases(obj,miniBatchFeatures{j},miniBatchOutcomes{j},learningRate,lmdba,miniBatchSize.*floor(size(trainingFeatures,2)./miniBatchSize),mu);

                     %Restore weights
                     if strcmpi(dropoutMethod,'Dropout')
                         for k = 1:length(obj.weights)
                             obj.weights{k}(:,indicesToRemoveCell{k}) = 0;
                             if k~=1
                                 obj.weights{k-1}(indicesToRemoveCell{k},:) = 0;
                                 obj.weights{k-1} = obj.weights{k-1} + removedWeights{k-1};
                             end
                             if k == length(obj.weights)
                                 obj.weights{k} = obj.weights{k} + removedWeights{k};
                                 
                             end
                         end
                     elseif strcmpi(dropoutMethod,'Dropconnect')
                         for k = 1:length(obj.weights)
                             obj.weights{k}(indicesToRemoveCell{k}) = 0;
                             obj.weights{k} = obj.weights{k} + removedWeights{k};
                         end
                     end
                     
%                      for k = 1:length(obj.weights)
%                          figure('Name',strcat(j,k,'Weights'))
%                          temp = obj.weights{k};
%                          temp(temp~=0) = 1;
%                          imagesc(temp);
%                          colormap gray
%                          colorbar
%                      end
                     
                end
                %strcat('Epoch ',num2str(i),' is complete')
                if (strcmpi(dropoutMethod,'Dropout') || strcmpi(dropoutMethod,'Dropconnect')) && i == epochs
                    %Networks has been trained using fewer nodes/weights than the fully connected to need to multiply all weights by a factor
                    for j = 1:length(obj.weights)
                        totalNumberOfWeights = (size(obj.weights{k},1)*size(obj.weights{k},2));
                        fractionToDropToGiveInteger = ceil(totalNumberOfWeights.*((percentToDrop(j))/100))./(totalNumberOfWeights);
                        obj.weights{j} = obj.weights{j}.*(1-fractionToDropToGiveInteger);
                    end
                end
            end
            
        end
        
        function obj = updateWeightsAndBiases(obj,miniBatchFeatures,miniBatchOutcomes,learningRate,lmdba,n,mu)
            nablaB = cell(size(obj.biases));
            nablaW = cell(size(obj.weights));
            for i = 1:length(nablaB)
                nablaB{i} = zeros(size(obj.biases{i}));
                nablaW{i} = zeros(size(obj.weights{i}));
            end
            for j = 1:size(miniBatchFeatures,2)
                [deltaNablaB, deltaNablaW] = backpropagate(obj,miniBatchFeatures(:,j),miniBatchOutcomes(:,j));
                for i = 1:length(nablaB)
                    nablaB{i} = nablaB{i}+deltaNablaB{i};
                    nablaW{i} = nablaW{i}+deltaNablaW{i};
                end
            end
            for i = 1:length(nablaB)
                obj.velocities{i} = ((mu).*(obj.velocities{i})) - ((learningRate(i)./size(miniBatchFeatures,2)).*nablaW{i});
                obj.biases{i} = obj.biases{i}  - ((learningRate(i)./size(miniBatchFeatures,2)).*nablaB{i});
                obj.weights{i} = obj.weights{i} + obj.velocities{i} - (learningRate(i).*(lmdba./n).*(obj.weights{i})) ;
            end
            
        end
        
        function [deltaNablaB, deltaNablaW] = backpropagate(obj,miniBatchFeatures,miniBatchOutcome)
            deltaNablaB = cell(size(obj.biases));
            deltaNablaW = cell(size(obj.weights));
            for i = 1:length((obj.weights))
                deltaNablaB{i} = zeros(size(obj.biases{i}));
                deltaNablaW{i} = zeros(size(obj.weights{i}));
            end
            %Feedforward
            activation = miniBatchFeatures;
            activations = cell(1,length(obj.biases)+1);
            activations{1} = activation;
            layerInputs = cell(size(obj.weights));
            for i = 1:length(obj.weights)
                layerInput = obj.weights{i}*activation + obj.biases{i};
                layerInputs{i} = layerInput;
                switch obj.layerTypes{i}
                    case 'sigmoid'
                        activation = sigmoid(layerInput);
                    case 'tanh'
                        activation = tanh(layerInput);
                    case 'rectifiedLinear'
                        activation = rectLinear(layerInput);
                    case 'leakyRectifiedLinear'
                        activation = leakyRectLinear(layerInput);
                    case 'softmax'
                        activation = softmax(layerInput);
                    case 'sigmoidPlusLinear'
                        activation = sigmoidPlusLinear(layerInput);
                    otherwise
                        error('Layer Types are wrong')
                end
                activations{i+1} = activation;
            end
            %backpass
            %NB THIS NEEDS EDITING TO ALLOW SOFTMAX IN THE OUTPUT LAYER AND
            %DIFFERENT COST FUNCTIONS
            [obj delta] = deltaFun(obj,activations{end},miniBatchOutcome,layerInputs{end});
            %delta = costDev.*sigmoidPrime(layerInputs{end});
            deltaNablaB{end} = delta;
            deltaNablaW{end} = delta*((activations{end-1})') ;
            for i = 1:length(obj.weights)-1
                switch obj.layerTypes{end-i}
                    case 'sigmoid'
                        activationFunPrime = sigmoidPrime(layerInputs{end-i});
                    case 'tanh'
                        activationFunPrime = 1-(tanh(layerInputs{end-i}).^2);
                    case 'rectifiedLinear'
                        activationFunPrime = rectLinearPrime(layerInputs{end-i});
                    case 'leakyRectifiedLinear'
                        activationFunPrime = leakyRectLinearPrime(layerInputs{end-i});
                    case 'sigmoidPlusLinear'
                        activationFunPrime = sigmoidPlusLinearPrime(layerInputs{end-i});
                    otherwise
                        error('Layer Types are wrong')
                end
                delta = (((obj.weights{end-i+1})')*delta).*activationFunPrime;
                deltaNablaB{end-i} = delta;
                deltaNablaW{end-i} = delta*((activations{end-i-1})');
            end
        end
        
        
        function [obj delta] = deltaFun(obj,outputActivation,outcome,z)
            switch obj.costFunction
                case 'Quadratic'
                    switch obj.layerTypes{end-i}
                        case 'sigmoid'
                            delta = (outputActivation - outcome).*sigmoidPrime(z);
                        case 'tanh'
                            delta = (outputActivation - outcome).*(1-(tanh(z).^2));
                        case 'rectifiedLinear'
                            delta = (outputActivation - outcome).*rectLinearPrime(z);
                        case 'leakyRectifiedLinear'
                            delta = (outputActivation - outcome).*leakyRectLinearPrime(z);
                        case 'sigmoidPlusLinear'
                            delta = (outputActivation - outcome).*sigmoidPlusLinearPrime(z);
                        otherwise
                            error('Layer Types are wrong')
                    end
                    
                case 'crossEntropy'
                    
                    if ~strcmp(obj.layerTypes{end},'sigmoid')
                        error('Do not use cross entropy cost function if final layer is not sigmoid')
                    end
                    %Cross Entropy Cost
                    delta = (outputActivation - outcome);
                    
                case 'logLikelyhood'
                    %Cross Entropy Cost
                    delta = (outputActivation - outcome);
                otherwise
                    error('Trying to call unsupported Cost function')
            end
            
        end
    end
end