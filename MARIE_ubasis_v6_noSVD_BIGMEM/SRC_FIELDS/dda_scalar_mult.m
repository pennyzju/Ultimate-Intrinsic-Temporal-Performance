function out = dda_scalar_mult(scalar, vector)
out = vector;
out(:,:,:,1) = out(:,:,:,1).*scalar;
out(:,:,:,2) = out(:,:,:,2).*scalar;
out(:,:,:,3) = out(:,:,:,3).*scalar;