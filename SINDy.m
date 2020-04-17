function [Xi, Opt_info] = SINDy(D,V,solver_type,Precond_info)
% Sparse Identification of Non-linear Dynamics for a system
% Inputs:
% D*Xi = V 
% D is an MxP dictionary matrix 
% V is a Mxd matrix of derivative values
% solver_type 'OMP' 'SP' 'SPGL1'
% Output:
% Xi is a Pxd matrix of sparse coefficient vectors
% warning('off','MATLAB:rankDeficientMatrix');
% warning('off','MATLAB:nearlySingularMatrix');
d = size(V,2);
[N,P] = size(D);
Xi = zeros(P,d);
%opts = spgSetParms('iterations',1e5,'verbosity',0,'optTol',1e-9,'bp');
if nargin < 4
    Precond_info = [];
end


for l = 1:d
    if strcmp(solver_type, 'OMP')
        Xi(:,l) = OMP(D,-inf,V(:,l),false);
    end
    if strcmp(solver_type,'SPGL1')
        Opt_info =  L_curve_bpdn(D,V(:,l),Precond_info);
        Xi(:,l) = Opt_info.solutions(:,Opt_info.min_ind);
    end
    if strcmp(solver_type,'LASSO')
        weights = get_matrix_weights(D);
        Dweights = D*weights;
        tau =  cross_val_tau(D,V(:,l));
        opts = spgSetParms('iterations',6000*size(D,2),'verbosity',0); %20*size(D,2)
        Xi(:,l) = weights*spg_lasso(Dweights,V(:,l),tau,opts);
    end
    if strcmp(solver_type,'SP')
        [Xi(:,l), Opt_info] =  SP(D,0,V(:,l),false,'L',Precond_info);
    end
    if strcmp(solver_type,'RR')
        Opt_info =  L_curve_RR(D,V(:,l),Precond_info);
        Xi(:,l) = Opt_info.solutions(:,Opt_info.min_ind);
    end
    if strcmp(solver_type,'IT')
        [Xi(:,l), Opt_info] = IT(D,V(:,l),[],Precond_info);
    end
    
end



    
end





