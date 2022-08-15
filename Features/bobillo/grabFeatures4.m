function features=grabFeatures4(T,U,ind)
%
%Classify a training set, tensor T, given basis matrices U{1},{2} and
%{3}, and vector "ind" of indices for features of vectorized core tensor.
%
s=size(T);
for i=1:s(4)
    Gi=tmprod(T(:,:,:,i),U,[1 2 3],'T');
    aux=tens2vec(Gi,[1 2 3])';
    if exist('ind')
       features(i,:)=aux(ind); %each row is an observation
    else
       features(i,:)=aux;
    end
end
