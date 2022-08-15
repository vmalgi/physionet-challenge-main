function Tf=fillTensor4(T)
%
% Given a Tensor of size #features x #states x #cycles x #records, it fills trailing
% zeros of short recordings by replicating series from the start.
%
st=size(T);
for i=1:st(4)
   L=max(find(T(1,4,:,i)~=0));
   fill=st(3)-L;
   folds=floor(fill/L);
   tail=fill-folds*L;
   for f=1:folds+1;
       Tf(:,:,(f-1)*L+1:f*L,i)=T(:,:,1:L,i);
   end
   Tf(:,:,(folds+1)*L+1:st(3),i)=T(:,:,1:tail,i);
end

