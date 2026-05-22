

function [f df]=f0(x)
% function [f df]=f0(x)

f=x(end,1);
df=[zeros(size(x,1)-1,1);1.0];

