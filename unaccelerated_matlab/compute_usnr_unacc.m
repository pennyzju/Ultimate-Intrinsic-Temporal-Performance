


function [usnr_map b1_map gsar_map]=compute_usnr_unacc(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)
% function [usnr_map b1_map gsar_map]=compute_usnr_unacc(nx,ny,nz,nchannels,datadir,b1target,mask,lossmatfilepref)

% read B1- map
b1maps=read_b1maps(datadir,1,nchannels,nx,ny,nz);

% read global SAR matrix
% eval( sprintf('load %s/%s.txt',datadir,lossmatfilepref) );
% eval( sprintf('lossmat=%s(:,1:2:end) + 1j*%s(:,2:2:end);',lossmatfilepref,lossmatfilepref) );
% lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
% q_inv=pinv(lossmat.');

lossmat=reshape( read_fdata( sprintf('%s/%s.dat',datadir,lossmatfilepref),-1 ),nchannels,2*nchannels );
lossmat=lossmat(:,1:2:end) + 1j*lossmat(:,2:2:end);
q=lossmat;
q_inv=pinv(q.');

% % make sure that it is PSD, i.e. remove small negative eigenvalues due to
% % machine accuracy
% [v d]=eig(lossmat);
% if min(diag(d))<0
%     lossmat=v*abs(d)/v; 
% end

% compute system matrix
usnr_map=zeros(nx,ny,nz);
if nargout>1
    b1_map=zeros(nx,ny,nz);
end
if nargout>2
    gsar_map=zeros(nx,ny,nz);
end
ind=find(mask>0);
[indx indy indz]=ind2sub( [nx ny nz],ind );
fprintf('\tComputing uSNR ');
sbuf=[];
for i=1:size(ind,1)
    
%     for dum=1:size(sbuf,2)
%         fprintf('\b');
%     end
%     sbuf=sprintf('[%d/%d]',i,size(ind,1));
%     fprintf('%s',sbuf);
    
    xpix=indx(i,1);
    ypix=indy(i,1);
    zpix=indz(i,1);
    A=compute_system_matrix(b1maps,xpix,ypix,zpix);
    x0=ones(2*nchannels,1);
   
    % test gradients
    % test_gradient(@(x)f0_v2(x,q),x0);
    % test_hessian(@(x)f0_v2(x,q),@(x)dualhess_v2(x,q),x0);
    
    % solve constrained optimization problem
    Aeq=[real(A) -imag(A);imag(A) real(A)];
    beq=[real(b1target);imag(b1target)];
    options = optimset('Largescale','off','Display','none','MaxFunEvals',10^5,'TolFun',1e-6,'TolCon',1e-6,'MaxIter',50000,...
            'Algorithm','interior-point','GradObj','on','GradConstr','on','SubproblemAlgorithm','ldl-factorization');
    options.Hessian='user-supplied';
    options.HessFcn=@(x,lambda)dualhess_v2(x,q,lambda);
    [x fval exitflag output lambda]=fmincon(@(x)f0_v2(x,q),x0,[],[],Aeq,beq,[],[],[],options);
    
    % results
    gsar=f0_v2(x,q);
    rf=x(1:nchannels,1) + 1j*x(nchannels+1:2*nchannels,1);
    b1=A*rf;
    usnr_map(ind(i))=abs(b1)/sqrt(gsar);
    if nargout>1
        b1_map(ind(i))=b1;
    end
    if nargout>2
        gsar_map(ind(i))=gsar;
    end
    
    % uSNR -- fast calc
    rf_f=q_inv*A'*pinv( A*q_inv*A' )*b1target;
    b1_f=A*rf_f;
    gsar_f=f0_v2([real(rf_f);imag(rf_f)],q);
    
    fprintf('\tuSNR -- pixel %d out of %d: slow=%e  fast=%e  slow-fast=%e\n',i,size(ind,1),usnr_map(ind(i)),abs(b1_f)/sqrt(gsar_f),usnr_map(ind(i))-abs(b1_f)/sqrt(gsar_f) );
    
end

fprintf('\n');




        
        
        




