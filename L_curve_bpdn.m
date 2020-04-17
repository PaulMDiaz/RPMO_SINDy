function [Opt_info] = L_curve_bpdn(A,b,Precond_info)
% This function attempts to find an optimal value of sigma for the bpdn
% problem: 
%    minimize ||x||_1 subject to ||A*x-b||_2 <= sigma
% Inputs:
% A: is an M x P sensing matrix, where M is the number of
% data samples in time and P is the number of basis functions.
% y: is an M x 1 data bector corrupted with noise.
% Output:
% sigma is an approximate of the optimal tolerance for the spg_bpdn problem

Opt_info = [];
N_sigmas = 20; %Number of sigma tolerances to attempt 
[N,P] = size(A); %Number of basis functions
Opt_info.fit_errs = zeros(N_sigmas,1);
Opt_info.l1_norms = zeros(N_sigmas,1);
Opt_info.l0_norms = zeros(N_sigmas,1);
Opt_info.solutions = zeros(P,N_sigmas);
maxiter =  min(20*N,30*P);
opts = spgSetParms('iterations',maxiter ,'verbosity',0); 
%Opt_info.sigmas = linspace(0.0005,0.1,N_sigmas); %10.^(-linspace(2,0,N_sigmas));
Opt_info.sigmas = 10.^(-linspace(4,0,N_sigmas));

if isempty(Precond_info)
    Dval = A;
    Vval = b;
    precond_weights = eye(P);
else
    Dval = Precond_info.D;
    Vval = Precond_info.V;
    precond_weights = Precond_info.weights;
end

weights = get_matrix_weights(A);

for i = 1:N_sigmas
    sigma = Opt_info.sigmas(i);
    try
       Xi = weights*spg_bpdn(A*weights,b,sigma*norm(b),opts);
        if isnan(max(Xi))
            Xi = ones(P,1);
            fprintf('Spgl1 failed (NaN) failed on sigma = %.g max_iter =  %.g \n', sigma, maxiter)
        end
    catch ME
       fprintf( [ME.message , ' failed on sigma = %.g max_iter =  %.g \n'], sigma, maxiter) 
       Xi = ones(P,1);
    end
    Xi = precond_weights*Xi;
    
    Opt_info.fit_errs(i) = norm(Dval*Xi - Vval)/norm(Vval);
    Opt_info.solutions(:,i) = Xi;
    Opt_info.l1_norms(i) = sum(abs(Xi));
    Opt_info.l0_norms(i) = length(find(abs(Xi)>0));
end
% figure
% plot(Opt_info.fit_errs(i),Opt_info.fit_errs(i),'-o')

Opt_info.l1_norm = norm(Opt_info.l1_norms);
Opt_info.fit_errs_norm = norm(Opt_info.fit_errs);

Opt_info.l1_norms = Opt_info.l1_norms/Opt_info.l1_norm;
Opt_info.fit_errs = Opt_info.fit_errs/Opt_info.fit_errs_norm; 

[Opt_info.l1_norms , I] = sort(Opt_info.l1_norms);
Opt_info.l0_norms = Opt_info.l0_norms(I); 
% 
% [Opt_info.l0_norms, I] = sort(Opt_info.l0_norms);
% Opt_info.l1_norms = Opt_info.l1_norms(I); 

Opt_info.fit_errs = Opt_info.fit_errs(I); %Opt_info.fit_errs = Opt_info.fit_errs/norm(Opt_info.fit_errs);
Opt_info.sigmas = Opt_info.sigmas(I);
Opt_info.solutions = Opt_info.solutions(:,I);

[~,min_ind] = min( sqrt((min(Opt_info.l1_norms)- Opt_info.l1_norms).^2 + (min(Opt_info.fit_errs) -Opt_info.fit_errs).^2));
sigma = Opt_info.sigmas(min_ind);
Opt_info.sigma = sigma;
Opt_info.min_ind = min_ind;
%save(['LORENZ_BPDN_', num2str(N)])

%figure
% plot(Opt_info.l1_norms,Opt_info.fit_errs,'-o', %Opt_info.l1_norms(Opt_info.sigma_ind),Opt_info.fit_errs(Opt_info.sigma_ind),'rs');
%xlabel('$|| \xi ||_1$','interpreter','latex')
%ylabel('$|| V  -D\xi ||_2$','interpreter','latex')
%title('SPGl1 validation','interpreter','latex')
% Opt_info.fit_errs 
% Opt_info.l1_norms 
% Opt_info.l0_norms 
%sigma = Opt_info.sigma
%sigma_ind = Opt_info.sigma_ind
% pause 
% close all

end

