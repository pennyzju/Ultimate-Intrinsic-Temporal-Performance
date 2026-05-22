

function x0=compute_starting_point(A,b1target,sarmats)
% function x0=compute_starting_point(A,b1target,sarmats)

[x0 flag]=pcg(A'*A,A'*b1target,1e-6,5000);
x0=[real(x0);imag(x0);0];
[c ceq dc dceq]=fn(A,x0,b1target,sarmats);
x0(end,1)=1.1*max(c);  % ensure feasibility






