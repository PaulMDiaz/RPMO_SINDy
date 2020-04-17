function [Opt_info] = L_curve_RR(A,b,Precond_info)
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
N_sigmas = 60; %Number of sigma tolerances to attempt 
[N,P] = size(A); %Number of basis functions
Opt_info.fit_errs = zeros(N_sigmas,1);
Opt_info.l1_norms = zeros(N_sigmas,1);
Opt_info.l0_norms = zeros(N_sigmas,1);
Opt_info.l2_norms = zeros(N_sigmas,1);
Opt_info.solutions = zeros(P,N_sigmas);
Opt_info.sigmas = 10.^(linspace(-3,2,N_sigmas));

if isempty(Precond_info)
    Dval = A;
    Vval = b;
    precond_weights = eye(P);
else
    Dval = Precond_info.D;
    Vval = Precond_info.V;
    precond_weights = Precond_info.weights;
end


for i = 1:N_sigmas
    sigma = Opt_info.sigmas(i);
    try 
        Xi = STRidge(A,b,5e-5 , 10, sigma,2 );
        Xi = precond_weights*Xi;
    catch ME
        fprintf( [ME.message , 'RR failed on sigma = \n'], sigma);
        Xi = ones(P,1);
    end
    Opt_info.fit_errs(i) = norm(Dval*Xi - Vval)/norm(Vval);
    Opt_info.solutions(:,i) = Xi;
    Opt_info.l1_norms(i) = sum(abs(Xi));
    Opt_info.l0_norms(i) = length(find(abs(Xi)>0));
    Opt_info.l2_norms(i) = norm(Xi,2);
end
% figure
% plot(Opt_info.fit_errs(i),Opt_info.fit_errs(i),'-o')

Opt_info.l0_norm = norm(Opt_info.l0_norms);
Opt_info.fit_errs_norm = norm(Opt_info.fit_errs);

Opt_info.l0_norms = Opt_info.l0_norms/Opt_info.l0_norm;
Opt_info.fit_errs = Opt_info.fit_errs/Opt_info.fit_errs_norm; 

[Opt_info.l0_norms , I] = sort(Opt_info.l0_norms);
Opt_info.fit_errs = Opt_info.fit_errs(I); %Opt_info.fit_errs = Opt_info.fit_errs/norm(Opt_info.fit_errs);
%Opt_info.l0_norms = Opt_info.l0_norms(I); 
Opt_info.l1_norms = Opt_info.l1_norms(I); 
Opt_info.sigmas = Opt_info.sigmas(I);
Opt_info.solutions = Opt_info.solutions(:,I);

[~,min_ind] = min( sqrt((min(Opt_info.l0_norms)- Opt_info.l0_norms).^2 + (min(Opt_info.fit_errs) -Opt_info.fit_errs).^2));
sigma = Opt_info.sigmas(min_ind);
Opt_info.sigma = sigma;
Opt_info.min_ind = min_ind;
%save(['LORENZ_RR_', num2str(N)])
% figure
% plot(Opt_info.l0_norms,Opt_info.fit_errs,'-o', Opt_info.l0_norms(Opt_info.min_ind),Opt_info.fit_errs(Opt_info.min_ind),'rs');
% xlabel('$|| \xi ||_0$','interpreter','latex')
% ylabel('$|| V  -D\xi ||_2$','interpreter','latex')
% title('RR validation','interpreter','latex')
% Opt_info.fit_errs 
% Opt_info.l2_norms
% Opt_info.l1_norms 
% Opt_info.l0_norms 
% pause 
% close all

end

