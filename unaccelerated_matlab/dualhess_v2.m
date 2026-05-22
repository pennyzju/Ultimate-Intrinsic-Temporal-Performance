

function h=dualhess_v2(xrf,Q,lambda)
% function h=dualhess_v2(xrf,Q,lambda)

nchannels=size(Q,1);
rf=xrf(1:nchannels) + 1j*xrf(nchannels+1:2*nchannels);

subh=2.0*Q;
h=[real(subh) -imag(subh);imag(subh) real(subh)];


