

function h=dualhess(x,lambda,sarmats)
% function h=dualhess(x,lambda,sarmats)

nchannels=size(sarmats,1);
nsarmats=size(sarmats,3);
nunk=size(x,1);

h=zeros(nunk,nunk);

for n=1:nsarmats
    h(1:2*nchannels,1:2*nchannels)=h(1:2*nchannels,1:2*nchannels) + lambda.ineqnonlin(n,1)*2*[real(sarmats(:,:,n)) imag(sarmats(:,:,n));-imag(sarmats(:,:,n)) real(sarmats(:,:,n))];
end
