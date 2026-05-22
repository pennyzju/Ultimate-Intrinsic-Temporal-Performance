cd front-1
load BASIS_Einc;
nfields = size(Easis,2);
nchannels = nfields;
Einc = zeros(nx,ny,nz,nchannels);
for i = 1:1
    
    if mod(i,100) == 0
        fprintf('[%d/%d]\n',i,nfields);
    end
    
    tmp( idxS )=Ebasis(:,i);    
    b1maps(:,:,:,i) = tmp;
end
load('object_def.mat')
Einc=zeros( size(mask) );
Einc2=zeros( size(mask) );
Einc(idxS) = Ebasis(1:193542,1);
Einc2(idxS) = Ebasis(193543:193542+193542,1);
volumeViewer
load('H:\UISNR\20240801_UISNR_output\front_4\front_4_material_maps.mat')
E_squared_weighted = abs(Einc).^2;
E_squared_weighted2 = abs(Einc2).^2;
E_squared_weighted_new = cat(4,E_squared_weighted ,E_squared_weighted2);
    % Element-wise multiplication of conductivity and weighted electric field intensity
    %alpha = sqrt(rotated_sigma_e .* E_squared_weighted);
alpha = rotated_sigma_e .* E_squared_weighted;
k = 7e6;
usnr_map2 = k*usnr_map;
    tSNR_map = usnr_map2./sqrt(1+alpha+coefficient_of_variation.*usnr_map2);
    k2 = 35e6;
usnr_map3 = k2*usnr_map;
    tSNR_map2 = usnr_map3./sqrt(1+alpha+coefficient_of_variation.*usnr_map3);
 k3 = 63e6;
usnr_map4 = k3*usnr_map;
    tSNR_map3= usnr_map4./sqrt(1+alpha+coefficient_of_variation.*usnr_map4);