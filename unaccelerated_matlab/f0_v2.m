

function [f df]=f0_v2(xrf,Q)
% function [f df]=f0_v2(xrf,Q)

nchannels=size(Q,1);
rf = xrf(1:nchannels,1) + 1j*xrf(nchannels+1:2*nchannels,1);

tmp=Q*rf;
f= real( rf'*tmp );

if nargout>1
    df = [ 2*real(tmp); 2*imag(tmp) ];
end



