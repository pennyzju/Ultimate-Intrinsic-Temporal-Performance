


function [c ceq dc dceq]=fn(A,x,b1target,sarmats)
% function [c ceq dc dceq]=fn(gamma,A,x,b1target,sarmats)

nchannels=size(A,2);
rf=x(1:nchannels,1) + 1j*x(nchannels+1:2*nchannels,1);
gamma=x(end,1);

% equality constraint
ceq=[];
dceq=[];

% inequality constraints
nsarmats=size(sarmats,3);
c=zeros(nsarmats,1);
dc=[];
for i=1:nsarmats
    tmp=sarmats(:,:,i)*conj(rf);
    c(i,1)=real( (rf.')*tmp ) - gamma;
    dc=[dc [2*real(tmp);-2*imag(tmp);-1.0]];
end



