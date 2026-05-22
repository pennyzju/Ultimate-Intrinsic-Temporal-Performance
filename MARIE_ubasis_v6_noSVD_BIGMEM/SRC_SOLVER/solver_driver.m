function [DATA_OUT] = solver_driver(DATA_IN)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that applies the solver on a set of different examples
%
%
%   INPUT:  DATA_IN     a vector of structures, each containing 
%
%                           DATA_IN.input           vector with excitation
%                           DATA_IN.circulant       circulant
%                           DATA_IN.epsilon         relative epsilon
%                           DATA_IN.domain          complete domain
%                           DATA_IN.frequency       frequency of operation
%                           DATA_IN.formulation 	formulation to use (1 or 2)
%                           DATA_IN.local           0 if all domain is used, 1 if only non-air is taken into account
%                           DATA_IN.method_code     method
%                                                   0 if Self coded 'GMRES_DR' (with deflated restart)
%                                                   1 if Self coded 'GPU Dedicated JVIE Solver based on GMRES Deflated Restarted'
%                                                   2 if Self coded 'GMRES_R' (with restart)
%                                                   3 if Self coded 'GPU Dedicated JVIE Solver based on GMRES Restarted'
%                                                   4 if Self coded 'GMRES'
%                                                   5 if Matlab built-in 'M_GMRES'
%                                                   6 if Matlab built-in 'M_GMRES_R' (with restart)
%                                                   7 if Matlab built-in 'M_QMR'
%                                                   8 if Matlab built-in 'M_TFQMR'
%                                                   9 if Matlab built-in 'M_BICG'
%                                                   10 if Matlab built-in 'M_BICGSTAB'
%                                                   11 if Self coded 'GPU Dedicated BICGSTAB'
%                           DATA_IN.tolerance       tolerance of the method
%                           DATA_IN.iterations      number of iterations
%                           DATA_IN.inner_it        number of iterations for restarted methods
%                           DATA_IN.outer_it        number of out cycles for restarted methods
%                           DATA_IN.ritz            size of subspace kept after each restart for deflated restarted methods
%                           DATA_IN.gpu_num         0 if no GPU is used
%                                                   if gpu_num > 0, it uses the designated gpu to accelerate the method (i.e. GPU 1 or GPU 2, or GPU 3...)
%                           DATA_IN.precond         0 if no preconditining
%                                                   1 if left preconditioning
%                                                   2 if right preconditioning
%                           DATA_IN.logfile         file to write info
%
%   OUTPUT: DATA_OUT    a vector of structures, each containing:
% 
%                           DATA_OUT.solution       vector with the solution
%                           DATA_OUT.flag           0 if converged, 1 if not
%                           DATA_OUT.relres         relative residue of the solution: norm(b - A*x)/norm(b)
%                           DATA_OUT.iterations     vector with [current internal iterations, current external iterations]
%                           DATA_OUT.resvec         vector containing norm of relative residual at each iteration
%                           DATA_OUT.time           time elapsed by the method
%
%
%
% _________________________________________________________________________
%
%   Computational Prototyping Group, RLE at MIT
% _________________________________________________________________________
% _________________________________________________________________________
%
%


% % % Exaple of structure, copy, paste and modified any field
% % currentdata.input = vExc;
% % currentdata.circulant = fG;
% % currentdata.epsilon = e_r;
% % currentdata.domain = r;
% % currentdata.frequency = freq;
% % currentdata.formulation = 2;
% % currentdata.local = 1;
% % currentdata.method = 0;
% % currentdata.tolerance = 1e-4;
% % currentdata.iterations = 300;
% % currentdata.inner_it = 30;
% % currentdata.outer_it = 20;
% % currentdata.ritz = 5;
% % currentdata.gpu_num = 0;
% % currentdata.precond = 0;
% % currentdata.logfile = './output.txt';
% % DATA_IN(1) = currentdata;


% get the number of experiments
num_experiments = length(DATA_IN);

% initialize the output structure
outdata.solution = [];
outdata.flag = [];
outdata.relres = [];
outdata.iterations = [];
outdata.resvec = [];
outdata.time = [];
% allocate the number of structures
DATA_OUT(num_experiments) = outdata;

for ii = 1:num_experiments
    
    % grab data and check
    EINC = DATA_IN(ii).input;
    FG = DATA_IN(ii).circulant;
    EPS = DATA_IN(ii).epsilon;
    R = DATA_IN(ii).domain;
    FREQ = DATA_IN(ii).frequency;
    FORMULATION = DATA_IN(ii).formulation;
    LOCAL = DATA_IN(ii).local;
    METHOD = DATA_IN(ii).method;
    TOL = DATA_IN(ii).tolerance;
    IT = DATA_IN(ii).iterations;
    INNER_IT = DATA_IN(ii).inner_it;
    OUTER_IT = DATA_IN(ii).outer_it;
    RITZ = DATA_IN(ii).ritz;
    GPU = DATA_IN(ii).gpu_num;
    PRECOND = DATA_IN(ii).precond;
    LOGFILE = DATA_IN(ii).logfile;
    
    % call the solver
    [JSOL,FLAG,RELRES,ITER,RESVEC,TIME] = SOLVE(EINC,FG,EPS,R,FREQ,FORMULATION,LOCAL,METHOD,TOL,IT,INNER_IT,OUTER_IT,RITZ,GPU,PRECOND,LOGFILE);
    
    % grab output data and encapsulate it
    DATA_OUT(ii).solution = JSOL;
    DATA_OUT(ii).flag = FLAG;
    DATA_OUT(ii).relres = RELRES;
    DATA_OUT(ii).iterations = ITER;
    DATA_OUT(ii).resvec = RESVEC;
    DATA_OUT(ii).time = TIME;
    
end
    
