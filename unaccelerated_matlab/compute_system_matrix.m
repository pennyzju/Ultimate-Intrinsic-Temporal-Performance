

function A=compute_system_matrix(b1maps,xpix,ypix,zpix)
% function A=compute_system_matrix(b1maps,xpix,ypix,zpix)

A=reshape( b1maps(xpix,ypix,zpix,:),1,size(b1maps,4) );

