function meanOneBeat = select_correlated_beats_and_mean(ManyBeats,dimensionToAverageOver)
if length(size(ManyBeats)) ~= 3 && length(size(ManyBeats)) ~= 2
    error('This method can only average over 2 or 3 dimensional feature vectors')
end

if length(size(ManyBeats)) == 3
    total = size(ManyBeats,1).*size(ManyBeats,2).*size(ManyBeats,3);
    %for a 3d vector concatenate the two dimensions you are not averaging
    %over
    eachBeat = zeros(size(ManyBeats,dimensionToAverageOver),total./size(ManyBeats,dimensionToAverageOver));
    for i = 1:size(ManyBeats,dimensionToAverageOver)
        if dimensionToAverageOver == 1
            eachBeat(i,:) =  reshape(ManyBeats(i,:,:),[1, total./size(ManyBeats,dimensionToAverageOver)]);
        elseif dimensionToAverageOver == 2
            eachBeat(i,:) =  reshape(ManyBeats(:,i,:),[1, total./size(ManyBeats,dimensionToAverageOver)]);
        elseif dimensionToAverageOver == 3
            eachBeat(i,:) =  reshape(ManyBeats(:,:,i),[1, total./size(ManyBeats,dimensionToAverageOver)]);
        else
            error('Dimension to average over must be 1,2,3')
        end
        
    end
    
elseif length(size(ManyBeats)) == 2
    if dimensionToAverageOver == 2
        eachBeat = ManyBeats';
    else
        eachBeat = ManyBeats;
    end
end


Closest = pdist(eachBeat);
Z = squareform(Closest);
Z(Z == 0) = NaN;
[minCol, iCol] = min(Z);
[minValue, iRow] = min(minCol);
FirstPoint = iRow;
SecondPoint = iCol(iRow);
extraPoint = [];
for k = 1:length(Z)
    if k ~= FirstPoint || k~=SecondPoint
        if (Z(FirstPoint,k) < 1.5*minValue && Z(SecondPoint,k) < 1.5*minValue)
            extraPoint = [extraPoint k];
        end
    end
end
HeartBeatsSelected = [extraPoint FirstPoint SecondPoint];

eachBeat = eachBeat(HeartBeatsSelected,:);
eachBeatMean = mean(eachBeat,1);

if length(size(ManyBeats)) == 3
    if dimensionToAverageOver == 1
        meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,2),size(ManyBeats,3)]);
    elseif dimensionToAverageOver == 2
        meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,1),size(ManyBeats,3)]);
    elseif dimensionToAverageOver == 3
        meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,1),size(ManyBeats,2)]);
    else
        error('Dimension to average over must be 1,2,3')
    end
    
    
elseif length(size(ManyBeats)) == 2
    meanOneBeat = eachBeatMean;
end
end