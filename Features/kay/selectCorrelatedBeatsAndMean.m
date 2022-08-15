function meanOneBeat = selectCorrelatedBeatsAndMean(ManyBeats)
Closest = pdist(ManyBeats);
Z = squareform(Closest);
Z(Z == 0) = NaN;
[minCol, iCol] = min(Z);
[minValue, iRow] = min(minCol);
FirstPoint = iRow;
SecondPoint = iCol(iRow);
extraPoint = [];
for k = 1:length(Z)
    if k ~= FirstPoint || k~=SecondPoint
        if (Z(FirstPoint,k) < 1.5*minValue && Z(FirstPoint,k) < 1.5*minValue)
            extraPoint = [extraPoint k];
        end
    end
end
HeartBeatsSelected = [extraPoint FirstPoint SecondPoint];
ManyBeats = ManyBeats(HeartBeatsSelected,:);
meanOneBeat = mean(ManyBeats,1);
end