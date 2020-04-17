function [Xi, Opt_info] = IT(D,V,lambda,Precond_info)
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');
Opt_info = [];
if isempty(lambda)
     %lambda = cross_val_lambda(D,V);
     Opt_info = L_curve_lambda(D,V,Precond_info);
     Xi = Opt_info.solutions(:,Opt_info.min_ind);
     return
end

Xi = D\V;
    for k = 1:10
        %biginds = setdiff(1:size(D,2),smallinds(:,ind));
        %Xi(biginds) = D(biginds)\V(biginds)
        for ind = 1:size(V,2)
            smallinds = find(abs(Xi(:,ind))<lambda);
            Xi(smallinds,ind) = 0;
            biginds = setdiff(1:size(D,2),smallinds); %~smallinds(:,ind); 
            Xi(biginds,ind) = D(:,biginds)\V(:,ind);
        end
    end
end


function lambda = cross_val_lambda(D,V)
R = 4; %Number of 'folds'
N = size(D,1);
N_r = floor(0.8*N); %Numer of reconstruction samples
N_lambdas = 30; %Number of tolerances to attempt
lambdas = linspace(0,-2,N_lambdas);
lambdas = 10.^(lambdas);
IT_errs = zeros(N_lambdas,1);
for i = 1:N_lambdas
    lambda = lambdas(i);
    val_errs = zeros(R,1);
    for r = 1:R
        %recon_inds = datasample(1:N,N_r,'Replace',false);
        %val_inds = setdiff(1:N,recon_inds);
        if r == R
            val_inds = (1+(r-1)*floor(R^(-1)*N)):N;
        else
            val_inds = (1+(r-1)*floor(R^(-1)*N)):r*floor(R^(-1)*N);
        end
        recon_inds = setdiff(1:N,val_inds);
        D_r = D(recon_inds,:);
        D_v = D(val_inds,:);
        V_r = V(recon_inds,:);
        V_v = V(val_inds,:);
        Xi = IT(D_r,V_r,lambda);
        val_errs(r) = norm(D_v*Xi- V_v)/norm(V_v);
    end
    IT_errs(i) = mean(val_errs);
end
[~,min_ind] = min(IT_errs);
lambda =  IT_errs(min_ind);

end 

function [Opt_info] = L_curve_lambda(A,b,Precond_info)
N_lambdas = 60; %Number of tolerances to attempt
%Opt_info.lambdas = 10.^(linspace(-3,2,N_lambdas));
Opt_info.lambdas = [0, 10.^(linspace(-3,2,N_lambdas))];
[N,P] = size(A); %Number of basis functionss
Opt_info.fit_errs = zeros(N_lambdas,1);
Opt_info.l1_norms = zeros(N_lambdas,1);
Opt_info.l0_norms = zeros(N_lambdas,1);
Opt_info.solutions = zeros(P,N_lambdas);
if isempty(Precond_info)
    Dval = A;
    Vval = b;
    precond_weights = eye(P);
else
    Dval = Precond_info.D;
    Vval = Precond_info.V;
    precond_weights = Precond_info.weights;
end

for i = 1:N_lambdas
    lambda = Opt_info.lambdas(i);
    try
        Xi = IT(A,b,lambda);
        Xi = precond_weights*Xi;
    catch ME
        fprintf( [ME.message , 'STLS failed on lambda = \n'], lambda);
        Xi = ones(P,1);
    end
    Opt_info.fit_errs(i) = norm(Dval*Xi - Vval)/norm(Vval);
    Opt_info.solutions(:,i) = Xi;
    Opt_info.l1_norms(i) = sum(abs(Xi));
    Opt_info.l0_norms(i) = length(find(abs(Xi)>0));
end

Opt_info.l1_norm = norm(Opt_info.l1_norms);
Opt_info.l0_norm = norm(Opt_info.l0_norms);
Opt_info.fit_errs_norm = norm(Opt_info.fit_errs);

Opt_info.l1_norms = Opt_info.l1_norms/Opt_info.l1_norm;
Opt_info.l0_norms = Opt_info.l0_norms/Opt_info.l0_norm;
Opt_info.fit_errs = Opt_info.fit_errs/Opt_info.fit_errs_norm; 

[Opt_info.l0_norms , I] = sort(Opt_info.l0_norms);
Opt_info.fit_errs = Opt_info.fit_errs(I); %Opt_info.fit_errs = Opt_info.fit_errs/norm(Opt_info.fit_errs);
Opt_info.l1_norms = Opt_info.l1_norms(I); 
Opt_info.lambdas = Opt_info.lambdas(I);
Opt_info.solutions = Opt_info.solutions(:,I);

[~,min_ind] = min( sqrt((min(Opt_info.l0_norms)- Opt_info.l0_norms).^2 + (min(Opt_info.fit_errs) -Opt_info.fit_errs).^2));
Opt_info.lambda = Opt_info.lambdas(min_ind);
Opt_info.min_ind = min_ind;

%save(['LORENZ_IT_', num2str(N)])
%figure
%plot(Opt_info.l0_norms,Opt_info.fit_errs,'-o', %Opt_info.l0_norms(Opt_info.lambda_ind),Opt_info.fit_errs(Opt_info.lambda_ind),'rs');
%xlabel('$|| \xi ||_0$','interpreter','latex')
%ylabel('$|| V  -D\xi ||_2$','interpreter','latex')
%title('STLS validation','interpreter','latex')
% Opt_info.fit_errs 
% Opt_info.l1_norms 
% Opt_info.l0_norms 
%lambda = Opt_info.lambda
%lambda_ind = Opt_info.lambda_ind
%pause 
%close all

end 