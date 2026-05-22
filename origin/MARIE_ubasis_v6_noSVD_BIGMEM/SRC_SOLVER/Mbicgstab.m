function [x,flag,relres,iter,resvec] = Mbicgstab(A,b,tol,maxit,x0,idxSS,Ldim,Mdim,Ndim,ixcut,iycut,izcut)
%MODIFIED BICGSTAB   BiConjugate Gradients Stabilized Method.



% Check for an acceptable number of input arguments
if nargin < 2
    error(message('MATLAB:bicgstab:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:bicgstab:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:bicgstab:RSHsizeMatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:bicgstab:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:bicgstab:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 4 || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,NaN);
    end
    return
end


if ((nargin == 5) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:bicgstab:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end


% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x);
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(n/50),10,n-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    rho1 = rho;
    rho = rt' * r;
    if (rho == 0.0) || isinf(rho)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    if ii == 1
        p = r;
    else
        beta = (rho/rho1)*(alpha/omega);
        if (beta == 0) || ~isfinite(beta)
            flag = 4;
            break
        end
        p = r + beta * (p - omega * v);
    end

    ph = p;

    v = iterapp('mtimes',afun,atype,afcnstr,ph);
    rtv = rt' * v;
    if (rtv == 0) || isinf(rtv)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    alpha = rho / rtv;
    if isinf(alpha)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    if abs(alpha)*norm(ph) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    xhalf = x + alpha * ph;        % form the "half" iterate
    s = r - alpha * v;             % residual associated with xhalf
    normr = norm(s);
    normr_act = normr;
    resvec(2*ii) = normr;
    
    % ---------------------------------------------------------------------
    % plot image
    iterplot(xhalf,idxSS,Ldim,Mdim,Ndim,ixcut,iycut,izcut,resvec(1:2*ii)./n2b);
    % ---------------------------------------------------------------------
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        s = b - iterapp('mtimes',afun,atype,afcnstr,xhalf);
        normr_act = norm(s);
        resvec(2*ii) = normr_act;
        if normr_act <= tolb
            x = xhalf;
            flag = 0;
            iter = ii - 0.5;
            resvec = resvec(1:2*ii);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                x = xhalf;
                resvec = resvec(1:2*ii);
                break;
            end
        end
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii);
        break
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = xhalf;
        imin = ii - 0.5;
    end
    
    sh = s;

    t = iterapp('mtimes',afun,atype,afcnstr,sh);
    tt = t' * t;
    if (tt == 0) || isinf(tt)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    omega = (t' * s) / tt;
    if isinf(omega)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    if abs(omega)*norm(sh) < eps*norm(xhalf)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = xhalf + omega * sh;        % x = (x + alpha * ph) + omega * sh
    r = s - omega * t;
    normr = norm(r);
    normr_act = normr;
    resvec(2*ii+1) = normr;

    % ---------------------------------------------------------------------
    % plot image
    iterplot(x,idxSS,Ldim,Mdim,Ndim,ixcut,iycut,izcut,resvec(1:2*ii+1)./n2b);
    % ---------------------------------------------------------------------

    
    % check for convergence        
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - iterapp('mtimes',afun,atype,afcnstr,x);
        normr_act = norm(r);
        resvec(2*ii+1) = normr_act;
        if normr_act <= tolb
            flag = 0;
            iter = ii;
            resvec = resvec(1:2*ii+1);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                resvec = resvec(1:2*ii+1);
                break;
            end
        end        
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break
    end
    
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
    r = b - iterapp('mtimes',afun,atype,afcnstr,xmin);
    if norm(r) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicgstab',tol,maxit,ii,flag,iter,relres);
end
