function [c_hat, Opt_info] = SP( Psi, K, u , weight,cross_val_type,Precond_info)
% Subspace Persuit
% Inputs: Psi  NxP measurement matrix
%       : K = approximate bound on signal sparsity such that K >= s.
%         set K = 0 if bound is unknown to cross validate K.
%       : u is an Nx1 vector of measurements (QOI samples) 
%       : weight == true for normalized support estimates
%       : cross_val_type, set this to 'random' for random resampling, or
%         set this to 'fold' for k-fold cross-validation
% Output: c Px1 vector of PCE coefficients
%INITIALIZATION:   
c_hat = zeros(size(Psi,2),1);
Opt_info = [];
stop = 0; num_iter=0; max_iter = 20*size(Psi,2);
if nargin < 5 
   cross_val_type = 'fold'; 
end
if nargin < 4 || weight == false
   W = eye(size(Psi,2)); 
   weight = false;
else 
   W = get_matrix_weights(Psi); 
   Psi =  Psi*W; 
end



if K == 0
    if strcmp('fold',cross_val_type)
        [K,~] = cross_val_SP_fold(Psi,u,weight);
    end
    if strcmp('random',cross_val_type)
        [K,~] = cross_val_SP(Psi,u,weight);
    end
    if strcmp('L',cross_val_type)
        [Opt_info] = L_curve_SP(Psi,u,weight,Precond_info);
        c_hat = Opt_info.solutions(:,Opt_info.min_ind);
        return
    end
end

%Initial support estimation
x = zeros(size(Psi,2),1);
corr = abs(Psi'*u); [vals,~] = sort(corr,'descend');
I = find(corr >= vals(K));
%Initial residual calculation
%x(I) =  Psi(:,I)\u;
x(I) =  pinv(Psi(:,I))*u;  
ur0 = u - Psi*W*x;

while stop == 0 
    corr = abs(Psi'*ur0);  [vals,~] = sort(corr,'descend'); %corr = abs(Psi'*ur0);
    II = find(corr >= vals(K));II = union(I,II);
    %LSP iteration. 
    x = zeros(size(Psi,2),1); x(II) = pinv(Psi(:,II))*u;  %x(II) = Psi(:,II)\u;  
    %Updated support estimation
    I0 = I;  [vals,~] = sort(abs(x),'descend');
    I =  find( abs(x) >= vals(K));  
    %Update Residual  
    x = zeros(size(Psi,2),1); x(I) = pinv(Psi(:,I))*u; %x(I) = Psi(:,I)\u;
    ur = u - Psi*W*x; 
    num_iter=num_iter+1;
    %Check stopping criteria
    if norm(ur,2) >= norm(ur0,2)                      
        I = I0;                                          
        stop = 1;
        %fprintf('SP stopped after %i iterations \n', num_iter)
    end
    if  num_iter == max_iter              
        stop = 1; 
        fprintf('SP hit maximum number of iterations \n')
    end
    ur0 = ur;
    %fprintf('SP residual = %f \n', norm(ur0))
   
end
%c_hat(I) = Psi(:,I)\u; %Alternatively Psi(:,I)\u;
c_hat(I) = pinv(Psi(:,I))*u; 
c_hat = W*c_hat;

end



function [Opt_info] = L_curve_SP(A,b, weight,Precond_info)
% L curve optimization function via Subspace Pursuit
% Inputs:
% Psi , NxP measurement matrix 
% u , Nx1 QoI vector
% Tol , tolernace for SP function



Opt_info = [];
[N,P] = size(A); %Number of basis functions
Opt_info.Ks = 1:P;%min(P,floor(N/2));
N_k = length(Opt_info.Ks);
Opt_info.fit_errs = zeros(N_k,1);
Opt_info.l1_norms = zeros(N_k,1);
Opt_info.l0_norms = zeros(N_k,1);
Opt_info.solutions = zeros(P,N_k);

if isempty(Precond_info)
    Dval = A;
    Vval = b;
   precond_weights = eye(P);
else
   Dval = Precond_info.D;
   Vval = Precond_info.V;
    precond_weights = Precond_info.weights;
end

if nargin < 3 
    weight = false;
end
for i = 1:N_k
    k = Opt_info.Ks(i);
    try 
        Xi = SP(A,k,b,weight);
        Xi = precond_weights*Xi;
    catch ME
        fprintf( [ME.message , 'SP failed on k = \n'], k);
        Xi = ones(P,1);
    end
    Opt_info.fit_errs(i) = norm(Dval*Xi - Vval)/norm(Vval);
    Opt_info.solutions(:,i) = Xi;
    Opt_info.l1_norms(i) = sum(abs(Xi));
    Opt_info.l0_norms(i) = length(find(abs(Xi)>0));
end
% figure
% plot(Opt_info.fit_errs(i),Opt_info.fit_errs(i),'-o')

Opt_info.l1_norm = norm(Opt_info.l1_norms);
Opt_info.l0_norm = norm(Opt_info.l0_norms);
Opt_info.fit_errs_norm = norm(Opt_info.fit_errs);

Opt_info.l1_norms = Opt_info.l1_norms/Opt_info.l1_norm;
Opt_info.l0_norms = Opt_info.l0_norms/Opt_info.l0_norm;
Opt_info.fit_errs = Opt_info.fit_errs/Opt_info.fit_errs_norm; 

[Opt_info.l0_norms , I] = sort(Opt_info.l0_norms);
Opt_info.fit_errs = Opt_info.fit_errs(I); %Opt_info.fit_errs = Opt_info.fit_errs/norm(Opt_info.fit_errs);
Opt_info.l1_norms = Opt_info.l1_norms(I); 
Opt_info.Ks = Opt_info.Ks(I);
Opt_info.solutions = Opt_info.solutions(:,I);

[~,min_ind] = min( sqrt((min(Opt_info.l0_norms)- Opt_info.l0_norms).^2 + (min(Opt_info.fit_errs) -Opt_info.fit_errs).^2));
K = Opt_info.Ks(min_ind);
Opt_info.K = K;
Opt_info.min_ind = min_ind;
% save(['LORENZ_SP_', num2str(N)])
%figure
%plot(Opt_info.l0_norms,Opt_info.fit_errs,'-o', %Opt_info.l0_norms(Opt_info.K_ind),Opt_info.fit_errs(Opt_info.K_ind),'rs');
%xlabel('$|| \xi ||_0$','interpreter','latex')
%ylabel('$|| V  -D\xi ||_2$','interpreter','latex')
%title('SP validation','interpreter','latex')
%xlim([0 1]); ylim([0 1]);
% Opt_info.fit_errs 
% Opt_info.l1_norms 
% Opt_info.l0_norms 
%K = Opt_info.K
%K_ind = Opt_info.K_ind
%pause 
%close all


end

