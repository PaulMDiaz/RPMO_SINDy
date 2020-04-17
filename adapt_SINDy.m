function [Xi_SPGL1, Xi_SP, Xi_RR, Xi_IT, Sys_info] = adapt_SINDy(D,V,Lambda_cands,store)
% Sparse Identification of Non-linear Dynamics for a system
% Inputs:
% D*Xi = V 
% D is an MxP dictionary matrix 
% V is a Mxd matrix of derivative values
% solver_type 'OMP' 'SP' 'SPGL1'
% Output:
% Xi is a Pxd matrix of sparse coefficient vectors
%warning('off','MATLAB:rankDeficientMatrix');
%warning('off','MATLAB:nearlySingularMatrix');

if nargin < 4
   store = false; 
end
d = size(V,2);
[N,P] = size(D);
Xi = zeros(P,d);
weights = get_matrix_weights(D);
Dbar = D*weights; 

Xi_SPGL1 = Xi;
Xi_SP = Xi;
Xi_RR = Xi;
Xi_IT = Xi;

SPGL1_best_err = inf*ones(1,d);
SP_best_err = inf*ones(1,d);
RR_best_err = inf*ones(1,d);
IT_best_err = inf*ones(1,d);

RR_info =  struct('solver','RR','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
SP_info =  struct('solver','SP','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
SPGL1_info =  struct('solver','SPGL1','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
IT_info =  struct('solver','IT','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
RR_Opt_info = {}; 
SP_Opt_info = {}; 
SPGL1_Opt_info = {}; 
IT_Opt_info = {}; 


paramCG = struct('Method','pcg','display','off','LS',0,'TolX',1e-22,'TolFun',1e-22,'MaxIter',200,'optTol',1e-22);
lambda = 0;
for i = 1:length(Lambda_cands)+1
   if i == length(Lambda_cands)+1
      Phi_iter = eye(size(D,1));
      lambda = 0 + 1i;
   else
       lambda = Lambda_cands(i);
       param_Robust = struct('M',size(D,1),'Psi',Dbar,'L_iter',20,'lambda',lambda, 'isfigure', false);
       param_Robust.paramCG = paramCG;
       %fprintf('Optimizing Robust Projection Matrix for lambda = %0.3g ... \n', lambda);
       %tic; 
       [Phi_iter,~,~,~] = Projection_matrix_I(param_Robust);
       %time = toc;
       %fprintf('Robust Projection matrix optimized for lambda = %0.3g, completed in %0.3g seconds \n', lambda, time);
   end
   
   
    for l = 1:d
        
       %fprintf('lambda = %0.g , solving component %0.g\n', lambda, l);
       try
           %tic; 
           [xi_SPGL1, OPT_INFO] =  SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'SPGL1');
           xi_SPGL1 = weights*xi_SPGL1;
           SPGL1_Opt_info{l} = OPT_INFO;
           %time = toc;
           %fprintf('SPGL1 solved in = %0.3g seconds \n', time);
       catch ME
           fprintf( [ME.message , ' failed on lambda = %.g while solving component %.g \n'], lambda, l) 
           xi_SPGL1 = Xi(:,1);
       end
       %tic;
       [xi_SP, OPT_INFO] =  SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'SP');
       xi_SP = weights*xi_SP;
       SP_Opt_info{l} = OPT_INFO;
       %time = toc;
       %fprintf('SP solved in = %0.3g seconds \n', time);
       %tic;
       [xi_RR, OPT_INFO] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'RR');
       xi_RR = weights*xi_RR;
       RR_Opt_info{l} = OPT_INFO;
       %time = toc;
       %fprintf('RR solved in = %0.3g seconds \n', time);
       %tic;
       [xi_IT, OPT_INFO] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'IT');
       xi_IT = weights*xi_IT;
       IT_Opt_info{l} = OPT_INFO;
       %time = toc;
       %fprintf('IT solved in = %0.3g seconds \n', time);
        
       %tic
       SPGL1_err = norm(D*xi_SPGL1-V(:,l))/norm(V(:,l));
       SP_err = norm(D*xi_SP-V(:,l))/norm(V(:,l));
       RR_err = norm(D*xi_RR-V(:,l))/norm(V(:,l));
       IT_err = norm(D*xi_IT-V(:,l))/norm(V(:,l));
       
       if SPGL1_err < SPGL1_best_err(l)
            Xi_SPGL1(:,l) = xi_SPGL1;
            SPGL1_best_err(l) = SPGL1_err;
            W_SPGL1 = Phi_iter;
            SPGL1_info.lambda(l) = lambda;
       end
       if SP_err < SP_best_err(l)
            Xi_SP(:,l) = xi_SP;
            SP_best_err(l) = SP_err;
            W_SP = Phi_iter;
            SP_info.lambda(l) = lambda;
       end
       if RR_err < RR_best_err(l)
            Xi_RR(:,l) = xi_RR;
            RR_best_err(l) = RR_err;
            W_RR = Phi_iter;
            RR_info.lambda(l) = lambda;
       end
       if IT_err < IT_best_err(l)
            Xi_IT(:,l) = xi_IT;
            IT_best_err(l) = IT_err;
            W_IT = Phi_iter;
            IT_info.lambda(l) = lambda;
       end
       
       
       if i == length(Lambda_cands)+1
            RR_info = equiv_dict_stats(Dbar,W_RR,RR_info);
            SP_info = equiv_dict_stats(Dbar,W_SP,SP_info);
            IT_info = equiv_dict_stats(Dbar,W_IT,IT_info);
            SPGL1_info = equiv_dict_stats(Dbar,W_SPGL1,SPGL1_info);
            RR_info.reg_sol(:,l) = xi_RR;
            SP_info.reg_sol(:,l) = xi_SP;
            SPGL1_info.reg_sol(:,l) = xi_SPGL1;
            IT_info.reg_sol(:,l) = xi_IT;
       end
        
    end
    
    
    if store
        Opt_info = {SPGL1_Opt_info, SP_Opt_info, RR_Opt_info, IT_Opt_info};
        save(['SINDy_info_d',num2str(d),'_P',num2str(P),'_N',num2str(N)] , 'Opt_info');
    end
   %time = toc;
   %fprintf('Error processing time %0.3g seconds \n', time);
    
end



SPGL1_info.best_err = SPGL1_best_err;
SP_info.best_err = SP_best_err;
IT_info.best_err = IT_best_err;
RR_info.best_err = RR_best_err;
% Xi_SPGL1 = weights*Xi_SPGL1;
% Xi_SP = weights*Xi_SP;
% Xi_RR = weights*Xi_RR;
% Xi_IT = weights*Xi_IT;
Sys_info = {SPGL1_info,SP_info, RR_info, IT_info};    

end





