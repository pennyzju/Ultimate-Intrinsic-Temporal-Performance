function [fK] = retrieve_fK(LfK,MfK,NfK,dx,f)
% _________________________________________________________________________
% _________________________________________________________________________
%
%  Function to retrieve circulant if exists, or to generate it if not
%
% _________________________________________________________________________
% _________________________________________________________________________


fid = 1;

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

co = 299792458;
lambda  = co/f;
ko = 2*pi/lambda;

Tesla = f/(42.6e6);

% -------------------------------------------------------------------------
%                 get circulant
% -------------------------------------------------------------------------


name_fK = sprintf('..\\DATA_CIRCULANT\\fK_%ddot%dmm_%dx%dx%d_%dT.mat',floor(dx*1e3),round(10*(dx*1e3-floor(dx*1e3))),LfK,MfK,NfK,round(Tesla));

if  (~exist(name_fK,'file')) % the file with the circulant does not exists
    
    %     t1 = clock;
        
    % create grid for the circulanta
    xcir = 0:dx:(LfK-1)*dx;
    ycir = 0:dx:(MfK-1)*dx;
    zcir = 0:dx:(NfK-1)*dx;
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the fK Circulant:   %dx%dx%d\n\n',length(xcir),length(ycir),length(zcir));
    
    % define the 3D grid
    rcir = grid3d(xcir,ycir,zcir);
    
    clear xcir; clear ycir; clear zcir;
    
    % isOpen = (matlabpool('size') > 0);
    isOpen = ~isempty(gcp('nocreate'));
    if (~isOpen)
        mycluster = parcluster;
        parpool('local',mycluster.NumWorkers);
    end

    
    % Generate Green's Tensor Kop_mn
    [Kop_mn] = Generate_Kop(rcir,ko,dx);
   
    clear rcir;
    
    %     t2 = clock;
    %     fprintf(fid, '\n Green Tensor computed. elapsed time %f', etime(t2,t1));
    %
    %     t1 = clock;
    
    % Generate Circulant matrix
    [fK] = Generate_Circulant_Kop(Kop_mn);
       
    %save(name_fK, 'fK', '-V7.3');
    
    %     fprintf(1,'\n Circulant computed. Elapsed time %g [sec]\n', etime(clock,t1));
    
%     matlabpool close;
    
else % get the dimensions we need
  
    load(name_fK);
    fid=1;
    fprintf(fid, '\n     fK Circulant loaded from file: \n');
    fprintf(fid, '            %s \n\n\n', name_fK);
    
end


infocir = whos('fK');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n     FFT Space Dimensions:       %dx%dx%d',size(fK,1),size(fK,2), size(fK,3));
fprintf(fid, '\n     Circulant Memory:           %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');

