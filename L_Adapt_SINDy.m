function [Xi_SPGL1, Xi_SP, Xi_RR, Xi_IT, Sys_info] = L_Adapt_SINDy(D,V,Lambda_cands,plotting,store)
% Sparse Identification of Non-linear Dynamics for a system
% Inputs:
% D*Xi = V 
% D is an MxP dictionary matrix 
% V is a Mxd matrix of derivative values
% solver_type 'OMP' 'SP' 'SPGL1'
% Output:
% Xi is a Pxd matrix of sparse coefficient vectors

if nargin < 5
    store = false;
end
if nargin < 4
    plotting = false;
end
    
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');
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

M = size(D,2);  % use 2 for ETF Minimization and 1 for Identity matrix minimization (for over-determined systems)
% M = min(size(D)) % for underdetermined systems.



RR_info =  struct('solver','RR','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
SP_info =  struct('solver','SP','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
SPGL1_info =  struct('solver','SPGL1','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
IT_info =  struct('solver','IT','Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[],'best_err',[],'reg_sol',zeros(P,d));
%Sys_info = cell{d,4};


paramCG = struct('Method','pcg','display','off','LS',0,'TolX',1e-22,'TolFun',1e-22,'MaxIter',5e3,'optTol',1e-22);
%paramCG  = struct('Method','cg','display','on','LS',0,'TolX',1e-20,'TolFun',1e-20,'MaxIter',5e3,'optTol',1e-20);
Lambda_cands = [Lambda_cands , (0+1i)];

SP_err = zeros(length(Lambda_cands),d);
SPGL1_err = zeros(length(Lambda_cands),d);
RR_err = zeros(length(Lambda_cands),d);
IT_err = zeros(length(Lambda_cands),d);
SP_obj = zeros(length(Lambda_cands),d);
SPGL1_obj = zeros(length(Lambda_cands),d);
RR_obj = zeros(length(Lambda_cands),d);
IT_obj = zeros(length(Lambda_cands),d);
XI_SP = {zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands))};
XI_SPGL1 = {zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands))};
XI_RR =  {zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands))};
XI_IT = {zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands)),zeros(P,length(Lambda_cands))};




W_mats = zeros(M,size(D,1),length(Lambda_cands));

%W_mats = {};


Mu = zeros(length(Lambda_cands),1);
Mu_avg = zeros(length(Lambda_cands),1);
WD_fro = zeros(length(Lambda_cands),1);
W_fro = zeros(length(Lambda_cands),1);
Cond_W = zeros(length(Lambda_cands),1);
Cond_WD = zeros(length(Lambda_cands),1);

SP_solve_info = repmat({struct()}, d, length(Lambda_cands));
SPGL1_solve_info = repmat({struct()}, d, length(Lambda_cands));
RR_solve_info = repmat({struct()}, d, length(Lambda_cands));
IT_solve_info = repmat({struct()}, d, length(Lambda_cands));


for i = 1:length(Lambda_cands)
   
   if i == length(Lambda_cands)
      Phi_iter = eye(size(D,1));
      Dbar = D;
      Precond_info = [];
      W_mats(:,:,i) = eye;
   else
       close all
       Dbar = D*weights; 
       lambda = Lambda_cands(i);
       param_Robust = struct('M',M,'Psi',Dbar,'L_iter',20,'lambda',lambda, 'isfigure',plotting);
       param_Robust.paramCG = paramCG;
       %fprintf('Optimizing Robust Projection Matrix for lambda = %0.3g ... \n', lambda);
       %tic; 
       %[Phi_iter,~,~,~] = Projection_matrix_I(param_Robust);
       [Phi_iter,~,~,~] = Projection_matrix_ETF(param_Robust);
       %time = toc;
       %fprintf('Robust Projection matrix optimized for lambda = %0.3g, completed in %0.3g seconds \n', lambda, time);
       %pause
       %close all
       W_mats(:,:,i) = Phi_iter;
       Precond_info.weights = weights;
       Precond_info.D = D;
   end
       [Mu(i),Mu_avg(i)] = coherence(Phi_iter*Dbar);
       WD_fro(i) = norm(Phi_iter*Dbar,'fro');
       W_fro(i) = norm(Phi_iter,'fro');
       Cond_W(i) = cond(Phi_iter);
       Cond_WD(i) = cond(Phi_iter*Dbar);
   if plotting
       my_table = table(Lambda_cands(1:i)',Mu(1:i), Mu_avg(1:i), Cond_WD(1:i), WD_fro(1:i), Cond_W(1:i), W_fro(1:i),...
       'VariableNames',{'Lambda','Mu','Mu_avg','Cond_WD','Norm_WD','Cond_W','Norm_W'});
       disp(my_table);
   end
   
   
    for l = 1:d
       if i ~= length(Lambda_cands)
           Precond_info.V = V(:,l);
       end
        
       %fprintf('lambda = %0.g , solving component %0.g\n', lambda, l);
       try
           %tic; 
           [xi_spgl1, SPGL1_solve_info{l,i}] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'SPGL1',Precond_info);
           XI_SPGL1{l}(:,i) = xi_spgl1;
           %time = toc;
           %fprintf('SPGL1 solved in = %0.3g seconds \n', time);
       catch ME
           fprintf( [ME.message , ' failed on lambda = %.g while solving component %.g \n'], lambda, l);
           XI_SPGL1{l}(:,i) = ones(P,1);
       end
       %tic;
       try
           %tic; 
          [xi_sp, SP_solve_info{l,i}] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'SP',Precond_info);
           XI_SP{l}(:,i) = xi_sp;
           %time = toc;
           %fprintf('SPGL1 solved in = %0.3g seconds \n', time);
       catch ME
           fprintf( [ME.message , ' failed on lambda = %.g while solving component %.g \n'], lambda, l);
           XI_SP{l}(:,i) = ones(P,1);
       end 
        
        
       %time = toc;
       %fprintf('SP solved in = %0.3g seconds \n', time);
       %tic;
       try
           [xi_rr, RR_solve_info{l,i}] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'RR',Precond_info);
           XI_RR{l}(:,i) = xi_rr;
       catch ME
           fprintf( [ME.message , ' failed on lambda = %.g while solving component %.g \n'], lambda, l);
           XI_RR{l}(:,i) = ones(P,1);
       end 
       
       %time = toc;
       %fprintf('RR solved in = %0.3g seconds \n', time);
       %tic;
       try
           [xi_it, IT_solve_info{l,i}] = SINDy(Phi_iter*Dbar,Phi_iter*V(:,l),'IT',Precond_info);
           XI_IT{l}(:,i) = xi_it;
       catch ME
           fprintf( [ME.message , ' failed on lambda = %.g while solving component %.g \n'], lambda, l);
           XI_IT{l}(:,i) = ones(P,1);
       end 
       %time = toc;
       %fprintf('IT solved in = %0.3g seconds \n', time);
        
       %tic
       SPGL1_err(i,l) = norm(D*XI_SPGL1{l}(:,i)-V(:,l))/norm(V(:,l));
       SP_err(i,l) = norm(D*XI_SP{l}(:,i)-V(:,l))/norm(V(:,l));
       RR_err(i,l) = norm(D*XI_RR{l}(:,i) -V(:,l))/norm(V(:,l));
       IT_err(i,l) = norm(D*XI_IT{l}(:,i)-V(:,l))/norm(V(:,l));
       

       
       SPGL1_obj(i,l) = sum(abs(XI_SPGL1{l}(:,i)));
       %SPGL1_obj(i,l) = 0;
       %SPGL1_obj(i,l) = length(find(XI_SPGL1{l}(:,i) ~= 0));   
       SP_obj(i,l) = length(find(XI_SP{l}(:,i) ~= 0));   
       RR_obj(i,l) = length(find(XI_RR{l}(:,i) ~= 0));
       IT_obj(i,l) = length(find(XI_IT{l}(:,i) ~= 0));
       %IT_obj(i,l) = sum(abs(XI_IT{l}(:,i)));
       %SP_obj(i,l) = sum(abs(XI_SP{l}(:,i)));
       %RR_obj(i,l) = sum(abs(XI_RR{l}(:,i)));
       if i == length(Lambda_cands)
            RR_info.reg_sol(:,l) = XI_RR{l}(:,i);
            SP_info.reg_sol(:,l) = XI_SP{l}(:,i);
            SPGL1_info.reg_sol(:,l) = XI_SPGL1{l}(:,i); 
            IT_info.reg_sol(:,l) = XI_IT{l}(:,i);
       end
        
    end
   %time = toc;
   %fprintf('Error processing time %0.3g seconds \n', time);
       if plotting
          
          my_table = table([SPGL1_err(i,:)' SPGL1_obj(i,:)' ],[SP_err(i,:)' SP_obj(i,:)' ], [RR_err(i,:)' RR_obj(i,:)'] , [IT_err(i,:)' IT_obj(i,:)' ],...
               'VariableNames',{'SPGL1','SP','STRidge','STLS'});
          disp(my_table);
       end
    
end

if store
    save(['L_solver_info_d',num2str(d),'_P',num2str(P),'_N',num2str(N),'_M',num2str(M)]);
end

for l = 1:d
    [Xi_SPGL1(:,l), SPGL1_info.lambda(l), ind] =  L_curve_rpmo(SPGL1_err(:,l),SPGL1_obj(:,l),Lambda_cands,XI_SPGL1{l},plotting,'SPGL1',SPGL1_solve_info(l,:));
    %[Xi_SPGL1(:,l), SPGL1_info.lambda(l), ind] =  L_curve_rpmo(Lambda_cands,plotting,'SPGL1',SPGL1_solve_info(l,:));
    W_SPGL1 = W_mats(:,:,ind); SPGL1_info = equiv_dict_stats(Dbar,W_SPGL1,SPGL1_info); 
    
    [Xi_SP(:,l), SP_info.lambda(l), ind] =  L_curve_rpmo(SP_err(:,l),SP_obj(:,l),Lambda_cands,XI_SP{l},plotting,'SP',SP_solve_info(l,:));
    %[Xi_SP(:,l), SP_info.lambda(l), ind] =  L_curve_rpmo(Lambda_cands,plotting,'SP',SP_solve_info(l,:));
    W_SP = W_mats(:,:,ind); SP_info = equiv_dict_stats(Dbar,W_SP,SP_info); 
    
    [Xi_RR(:,l), RR_info.lambda(l), ind] =  L_curve_rpmo(RR_err(:,l),RR_obj(:,l),Lambda_cands,XI_RR{l},plotting,'RR',RR_solve_info(l,:));
    %[Xi_RR(:,l), RR_info.lambda(l), ind] =  L_curve_rpmo(Lambda_cands,plotting,'RR',RR_solve_info(l,:));
    W_RR = W_mats(:,:,ind); RR_info = equiv_dict_stats(Dbar,W_RR,RR_info); 
    
    [Xi_IT(:,l), IT_info.lambda(l), ind] =  L_curve_rpmo(IT_err(:,l),IT_obj(:,l),Lambda_cands,XI_IT{l},plotting,'IT',IT_solve_info(l,:));
    %[Xi_IT(:,l), IT_info.lambda(l), ind] =  L_curve_rpmo(Lambda_cands,plotting,'IT',IT_solve_info(l,:));
    W_IT = W_mats(:,:,ind); IT_info = equiv_dict_stats(Dbar,W_IT,IT_info); 
end



% RR_info = equiv_dict_stats(Dbar,W_RR,RR_info);
% SP_info = equiv_dict_stats(Dbar,W_SP,SP_info);
% IT_info = equiv_dict_stats(Dbar,W_IT,IT_info);
% SPGL1_info = equiv_dict_stats(Dbar,W_SPGL1,SPGL1_info);   
% 
% SPGL1_info.best_err = SPGL1_best_err;
% SP_info.best_err = SP_best_err;
% IT_info.best_err = IT_best_err;
% RR_info.best_err = RR_best_err;
% 
% if SPGL1_err < SPGL1_best_err(l)
%     Xi_SPGL1(:,l) = xi_SPGL1;
%     SPGL1_best_err(l) = SPGL1_err;
%     W_SPGL1 = Phi_iter;
%     SPGL1_info.lambda(l) = lambda;
% end
% if SP_err < SP_best_err(l)
%     Xi_SP(:,l) = xi_SP;
%     SP_best_err(l) = SP_err;
%     W_SP = Phi_iter;
%     SP_info.lambda(l) = lambda;
% end
% if RR_err < RR_best_err(l)
%     Xi_RR(:,l) = xi_RR;
%     RR_best_err(l) = RR_err;
%     W_RR = Phi_iter;
%     RR_info.lambda(l) = lambda;
% end
% if IT_err < IT_best_err(l)
%     Xi_IT(:,l) = xi_IT;
%     IT_best_err(l) = IT_err;
%     W_IT = Phi_iter;
%     IT_info.lambda(l) = lambda;
% end




% Xi_SPGL1 = weights*Xi_SPGL1;
% Xi_SP = weights*Xi_SP;
% Xi_RR = weights*Xi_RR;
% Xi_IT = weights*Xi_IT;
Sys_info = {SPGL1_info,SP_info, RR_info, IT_info};    

my_table = table(Lambda_cands',Mu, Mu_avg, Cond_WD, WD_fro, Cond_W, W_fro,...
       'VariableNames',{'Lambda','Mu','Mu_avg','Cond_WD','Norm_WD','Cond_W','Norm_W'});
       disp(my_table);


      
       
end





