function [fG] = retrieve_fG(LfG,MfG,NfG,dx,f)
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

% Tesla = f/(42.6e6);

% -------------------------------------------------------------------------
%                 get circulant
% -------------------------------------------------------------------------


% name_fG = sprintf('..\\DATA_CIRCULANT\\fG_%ddot%dmm_%dx%dx%d_%dp%dT.mat',floor(dx*1e3),round(10*(dx*1e3-floor(dx*1e3))),LfG,MfG,NfG,floor(Tesla), round(10*(Tesla) - 10*floor(Tesla)));
% 
% if  (~exist(name_fG,'file')) % the file with the circulant does not exists
    
    t1 = tic;
        
    % create grid for the circulanta
    xcir = 0:dx:(LfG-1)*dx;
    ycir = 0:dx:(MfG-1)*dx;
    zcir = 0:dx:(NfG-1)*dx;
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n     Generating the fG Circulant:   %dx%dx%d\n\n',length(xcir),length(ycir),length(zcir));
    
    % define the 3D grid
    rcir = grid3d(xcir,ycir,zcir);
    
    clear xcir; clear ycir; clear zcir;
    
    % isOpen = (matlabpool('size') > 0);
    isOpen = ~isempty(gcp('nocreate'));
    if (~isOpen)
        mycluster = parcluster;
        parpool('local',mycluster.NumWorkers);
    end
    
    % Generate Green's Tensor G_mn
    [Gmn] = Generate_Nop(rcir,ko,dx);
    
    clear rcir;
    
    %     t2 = clock;
    %     fprintf(fid, '\n Green Tensor computed. elapsed time %f', etime(t2,t1));
    %
    %     t1 = clock;
    
    % Generate Circulant matrix
    [fG] = Generate_Circulant_Nop(Gmn);
    
   
%     save(name_fG, 'fG', '-V7.3');
%     
    fprintf(1,'\n Circulant computed. Elapsed time %g [sec]\n', toc(t1));
%     
%     matlabpool close;
%     
% else % get the dimensions we need
%   
%     load(name_fG);
%     fid=1;
%     fprintf(fid, '\n     fG Circulant loaded from file: \n');
%     fprintf(fid, '            %s \n\n\n', name_fG);
%     
% end


infocir = whos('fG');
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n     FFT Space Dimensions:       %dx%dx%d',size(fG,1),size(fG,2), size(fG,3));
fprintf(fid, '\n     Circulant Memory:           %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n\n ');

