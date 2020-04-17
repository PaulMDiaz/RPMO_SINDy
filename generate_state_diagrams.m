
clear;%clc;
close all
addpath(genpath('./'))
rng('shuffle'); % Shuffles on local worker only
seed_offset = randi(floor(intmax/10),1);
timerval = tic 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% MODEL SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
func = @mean_field_model
ref_sol_fun = @MFM_ref_sol
X0 = [0.1; 0.1; 20]
TF = 5
M = 5000
Lambda = 1e-16
sample_inc = 1
Lambda_cands = 0.5*10.^linspace(-4,0,5)
sample_start = 1;
eta = 0.01 %Additive noise level for data and derivatives
ETA = 0.01
tik_lambda_cands_1 =  10.^(linspace(-16,-6,100)); 
tik_lambda_cands_2 =  10.^(linspace(-18,-6,100)); 
tik_lambda_cands_3 =  10.^(linspace(-18,-7,100)); 


% func = @mean_field_model
% ref_sol_fun = @MFM_ref_sol
% X0 = [0.1; 0.1; 20]
% TF = 4
% M =  200
% sample_start = 1;
% Lambda = 1e-16
% sample_inc = 20
% %Lambda_cands = 5*10.^(-3:1:-1) 
% Lambda_cands = 5*10.^(linspace(-7,-1,7)) %10.^(-8:0.25:-7)%linspace(10^(-8),10^(-1),8)%, 0.25, 0.5, 0.75]%10.^([-5:1:-3]) % 10.^(-2)
% %Lambda_cands = 0.00001



% 
% func = @MS
% ref_sol_fun = @MS_ref_sol
% TF = 4
% %X0 = [0,2,-1.2]';
% M = 120
% Lambda = 1e-16
% sample_inc = 20
% Lambda_cands =  0.5*10.^linspace(-4,0,50) %[0, 10.^linspace(-3,-2,10), linspace(0.02,1,30)]
% sample_start = 1;
% % X0 = [-0.01,0.01,0.6]';
% % tik_lambda_cands_1 =  10.^(linspace(-14,-4,100)); 
% % tik_lambda_cands_2 =  10.^(linspace(-14,-10,100)); 
% % %tik_lambda_cands_3 =  10.^(linspace(-16,-9,100)); 
% % %tik_lambda_cands_3 =  10.^(linspace(-12,0,200)); 
% % tik_lambda_cands_3 =  10.^(linspace(-14,-2,100)); 
% 
% eta = 0.01 %Additive noise level for data and derivatives
% tik_lambda_cands_1 =  10.^(linspace(-14,-6,100));
% tik_lambda_cands_2 =  10.^(linspace(-14,-6,100));
% tik_lambda_cands_3 =  10.^(linspace(-10,-3,200));
% 
% X0 = [-0.1,0.1,0.6]';

% eta = 0.005 %Additive noise level for data and derivatives
% tik_lambda_cands_1 =  10.^(linspace(-16,-6,100)); 
% tik_lambda_cands_2 =  10.^(linspace(-16,-6,100)); 
% tik_lambda_cands_3 =  10.^(linspace(-15,-4,100)); 




%400, 5, 10.^(linspace(-10,-8,30));
% func = @MS
% ref_sol_fun = @MS_ref_sol
% TF = 4
% X0 = [-1,0,0.3]';
% M = 19
% Lambda = 1e-17 
% sample_inc = 1


% 
% func = @lorenz63
% ref_sol_fun = @lorenz_ref_sol
% X0 = [-8;7;27];
% TF = 4
% M = 120
% Lambda = 1e-16
% eta = 0.01 %Additive noise level for data and derivatives
% sample_inc = 20
% Lambda_cands = 0.5*10.^linspace(-5,0,50)
% 
% tik_lambda_cands_1 =  10.^(linspace(-18,-8,100));
% tik_lambda_cands_2 =  10.^(linspace(-18,-8,100));
% tik_lambda_cands_3 =  10.^(linspace(-18,-8,100));
% sample_start = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotting = true;
store = false;

%%%%%
sample_sizes = M;
num_samps = 1;
n = 1;
R = 1;
p = 4 %/Total order of basis
d = 3 %Number of random variables
index_pc = nD_polynomial_array(d,p);
index_pc_M = nD_polynomial_array(d,p);

P = size(index_pc,1)
ETA = 0.01
%lambda =;%
%5*10.^(-5:0.5:2) %[5*10.^(-7:0.5:-1)]  %[0.005,0.05,0.1:0.025:0.5] %%10.^(-8:0.25:-7)%linspace(10^(-8),10^(-1),8)%, 0.25, 0.5, 0.75]%10.^([-5:1:-3]) % 10.^(-2)
%Lambda_cands =  5*10.^(-9:0.5:-5)
precond_type = 'none' %'robust', 'approx' (ksdensity), or 'none'% 

noisy_data = true
noisy_measurements = true
if eta ~= 0 && (noisy_data || noisy_measurements)
   noise = 'noisy'
else
   noise = 'noise-free'
end
gradient_type = 'exact' % TV, tik, or exact
data_sampling =  'continuous'  %'continuous' or  'random' 'dopt'  chr
noise_type = 'fixed' % 'added'

polytype = 'L' %Can set to 'L' but the coefficients will not match exact solution
T0 = 0


dt = 0.001


%w = warning('off','MATLAB:rankDeficientMatrix');
%close all
fprintf('Running...\n')
xi = 2*rand(10000,d) - 1;  
for isim = 1:10000
    psi_L(isim,:) = piset(xi(isim,:),index_pc);
    psi_M(isim,:) = piset_monomial(xi(isim,:),index_pc_M);
end
clear xi
T = [T0:dt:TF];



MuD_av = zeros(R,1);
MuD = zeros(R,1);
normD = zeros(R,1);
inf_normD = zeros(R,1);
Cond_D = zeros(R,1);

Mu_WSPGL1 = zeros(R,1);
Mu_WSPGL1_av = zeros(R,1);
cond_WSPGL1D = zeros(R,1);
norm_WSPGL1 = zeros(R,1);
fro_norm_W_SPGL1 = zeros(R,1);
inf_norm_WSPGL1 = zeros(R,1);
cond_WSPGL1 = zeros(R,1);
lambda_WSPGL1 = zeros(R,1);

Mu_WIT = zeros(R,1);
Mu_WIT_av = zeros(R,1);
cond_WITD = zeros(R,1);
norm_WIT = zeros(R,1);
fro_norm_W_IT = zeros(R,1);
inf_norm_WIT = zeros(R,1);
cond_WIT = zeros(R,1);
lambda_WIT = zeros(R,1);

Mu_WRR = zeros(R,1);
Mu_WRR_av = zeros(R,1);
cond_WRRD = zeros(R,1);
norm_WRR = zeros(R,1);
fro_norm_W_RR = zeros(R,1);
inf_norm_WRR = zeros(R,1);
cond_WRR = zeros(R,1);
lambda_WRR = zeros(R,1);

Mu_WSP = zeros(R,1);
Mu_WSP_av = zeros(R,1);
cond_WSPD = zeros(R,1);
norm_WSP = zeros(R,1);
fro_norm_W_SP = zeros(R,1);
inf_norm_WSP = zeros(R,1);
cond_WSP = zeros(R,1);
lambda_WSP = zeros(R,1);

SPGL1_err = zeros(R,1);
SPGL1_model_err = zeros(R,1);
SPGL1_stats = zeros(num_samps,12);
SPGL1_fit = zeros(R,1);

WSPGL1_err = zeros(R,1);
WSPGL1_model_err = zeros(R,1);
WSPGL1_stats = zeros(num_samps,15);
WSPGL1_fit = zeros(R,1);

IT_err = zeros(R,1);
IT_model_err = zeros(R,1);
IT_stats = zeros(num_samps,12);
IT_fit = zeros(R,1);

WIT_err = zeros(R,1);
WIT_model_err = zeros(R,1);
WIT_stats = zeros(num_samps,15);
WIT_fit = zeros(R,1);

RR_err = zeros(R,1);
RR_model_err = zeros(R,1);
RR_stats = zeros(num_samps,12);
RR_fit = zeros(R,1);

WRR_err = zeros(R,1);
WRR_model_err = zeros(R,1);
WRR_stats = zeros(num_samps,15);
WRR_fit = zeros(R,1);


SP_err = zeros(R,1);
SP_model_err = zeros(R,1);
SP_stats = zeros(num_samps,12);
SP_fit = zeros(R,1);

WSP_err = zeros(R,1);
WSP_model_err = zeros(R,1);
WSP_stats = zeros(num_samps,15);
WSP_fit = zeros(R,1);


timerval = tic;
r = R
       
        T = [T0:dt:TF];
        rng(seed_offset/r);
        
        meas_inds = [];
        m0 = 1;
        W = eye(M);    
        W_LSA = eye(M);
        V = zeros(M,P);
        if eta == 0 ||  noisy_measurements
            X0 = X0 + abs(X0).*randn(d,1)*0.0001;
        end
        [T,XX,V_Ref,D_M,D_L,scaling,A,B] = solve_ode(X0,T0:dt:TF,func,index_pc);
        dTA = [0; diff(T)];
        X_ref = XX;
        figure
        color_line3(XX(:,1),XX(:,2),XX(:,3),dTA,'LineWidth',4);
        view(27,16)
        axis equal;
        grid;
        xlabel('$x$','interpreter','latex','fontsize',24)
        ylabel('$y$','interpreter','latex','fontsize',24)

        zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
        title_string = 'Exact Solution'
        title(title_string,'interpreter','latex','fontsize',24);
        set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)
        set(gca,'FontSize',28)
        drawnow

        X0_ref_scaled = XX(1,:);
        X0_val = inv_shift_and_scale(XX(end,:),A,B);
        
        [T_val,X_val,~,D_val,~,~,~,~] = solve_ode(X0_val,T0:dt:TF,func,index_pc);

        %val_inds = (floor(size(D_M,1)/2)+2):size(D_M,1);
        %T_val = T(val_inds)';
        %D_val = D_M(val_inds,:);
        %V_val = V_Ref(val_inds,:);
        %X_val = XX(val_inds,:);
        
%         train_inds = 1:(floor(size(D_M,1)/2)+1);
%         T = T(train_inds);
%         D_M = D_M(train_inds,:); 
%         V_Ref = V_Ref(train_inds,:);
%         XX = XX(train_inds,:);
        

        Xi_ref = ref_sol_fun(A,B,P);
        D_val = [D_M(M+1:end,:);D_val];
        V_val = D_val*Xi_ref;
        % Add noise to the data
        
        dXX = zeros(size(XX));
        if noisy_data == true && eta ~= 0
            XX = inv_shift_and_scale(XX,A,B);
            if strcmp(noise_type,'added') 
                dXX = eta*abs(XX).*randn(size(XX));
                XX = XX+ dXX;
            end
            if strcmp(noise_type,'fixed')
                dXX = eta.*randn(size(XX));
                XX = XX+ dXX;
            end
            [XX,scaling] = shift_and_scale_data(XX);
        end
        if max(XX(:)) > 1 || min(XX(:)) < -1
            fprintf('Warning: Data lies outside [-1,1] support!')
        end
        if polytype == 'L'
           [~,D] = build_dictionary(XX,index_pc);
           [~,D_dXX] = build_dictionary(dXX,index_pc);
           
        else
           [D,~] = build_dictionary(XX,index_pc);
           [D_dXX,~] = build_dictionary(dXX,index_pc);
        end
        %D = D_L;
        sub_sample_inds = sample_start:sample_inc:size(XX,1);
        
        X_ref = X_ref(sub_sample_inds,:);
        dXX =dXX(sub_sample_inds,:);
        XX = XX(sub_sample_inds,:);
        D_dXX = D_dXX(sub_sample_inds,:);
        
        D_M = D_M(sub_sample_inds,:);
        D_L = D_L(sub_sample_inds,:);
        D = D(sub_sample_inds,:);
        D_true_norm = norm(D_L,2);
        D_noise_norm = norm(D_dXX,2);
        
        T = T(sub_sample_inds,:);
        V_Ref = V_Ref(sub_sample_inds,:);
        Tm = T(1:end-1)+ (T(2:end)- T(1:end-1))./2;
        

        dV = [];
        if strcmp(gradient_type,'TV')
            NN =  3:size(XX,1) - 2; 
            V_RRef = V_Ref(NN,:);
            figure
            plot(NN,V_RRef(:,1),NN,V_RRef(:,2),NN,V_RRef(:,3));
            hold on
            iter = 10;
            if eta == 0 
                 alph = 1e-4;
                 scale = 'small';
                 ep = 10000;
            else
                 alph = 1e-4;
                 scale = 'small';
                 ep = 10000;
            end
            dxdt_m =   TVRegDiff( XX(:,1), iter, alph, [], scale, ep, sample_inc*dt, 0, 0 );
            dxdt_m = dxdt_m(2:end-1);
            dydt_m =   TVRegDiff( XX(:,2), iter, alph, [], scale, ep, sample_inc*dt, 0, 0 );
            dydt_m = dydt_m(2:end-1);
            dzdt_m =   TVRegDiff( XX(:,3), iter, alph, [], scale, ep, sample_inc*dt, 0, 0 );
            dzdt_m = dzdt_m(2:end-1);
            
            dx_TV = pchip(Tm,dxdt_m,T); dx_TV = dx_TV(NN);
            dy_TV = pchip(Tm,dydt_m,T); dy_TV = dy_TV(NN);
            dz_TV = pchip(Tm,dzdt_m,T); dz_TV = dz_TV(NN);
            D = D(NN,:);
            
            V = [dx_TV, dy_TV,dz_TV ];
            plot(NN,V(:,1),NN,V(:,2),NN,V(:,3));
            dV = V - V_RRef;
            FD_approx_error = norm(dV,2)/norm(V_RRef,2)
            
            fprintf(' FD err =  %f \n', FD_approx_error)
%             toc(timerval)
        %     return
        end
        if strcmp(gradient_type,'tik')
            NN = size(XX,1) - 1; 
            V_RRef = V_Ref(1:NN,:);
            %figure
            dxdt_m = tik_diff(XX(:,1),T,[],plotting,tik_lambda_cands_1 );
            dxdt_m = dxdt_m(1:end-1);
            dydt_m = tik_diff(XX(:,2),T,[],plotting,tik_lambda_cands_2 );
            dydt_m = dydt_m(1:end-1);
            dzdt_m = tik_diff(XX(:,3),T,[],plotting,tik_lambda_cands_3 );
            dzdt_m = dzdt_m(1:end-1);
            Tm = Tm(1:end);
            T = T(1:end-1);
            D = D(1:end-1,:);
            dxdt = pchip(Tm,dxdt_m,T);
            dydt = pchip(Tm,dydt_m,T);
            dzdt = pchip(Tm,dzdt_m,T);
            V = [dxdt, dydt,dzdt];
            plot(1:NN,V_RRef(:,1),1:NN,V_RRef(:,2),1:NN,V_RRef(:,3));
            hold on
            plot(1:NN,V(:,1),1:NN,V(:,2),1:NN,V(:,3));
            drawnow
            dx1_err = norm(V(:,1)-V_RRef(:,1))/norm(V_RRef(:,1));
            dx2_err = norm(V(:,2)-V_RRef(:,2))/norm(V_RRef(:,2));
            dx3_err = norm(V(:,3)-V_RRef(:,3))/norm(V_RRef(:,3));
            FD_approx_error = mean([dx1_err, dx2_err, dx3_err])
            
            fprintf(' FD err =  %f, dxdt err = %f, dydt err = %f, dzdt err =  %f\n', FD_approx_error,  dx1_err ,  dx2_err ,  dx3_err )
            %pause
            %toc(timerval)
            %return
        end
      
        V_ref = V_Ref;
        if strcmp(gradient_type,'exact')
            V = V_ref;
        end
        
        tF = floor(size(V_Ref,1)/2); %number of timesteps in possible reconstruction data
        
        if strcmp(data_sampling,'continuous')
            meas_inds = 1:M; %chooses M contiguous in time samples
            %meas_inds = 1:sample_inc:size(D,2);
            %meas_inds = floor(linspace(1,tF,M)); %choose M samples uniformly
            %spaced in time
        end
        if  strcmp(data_sampling,'random')
            meas_inds = datasample(1:tF,M,'replace',false);
        end
        if strcmp(data_sampling,'dopt')
            meas_inds = rrqr_dopt_adapt(D,M,[]);
        end
        
        D = D(meas_inds,:);
        XX = XX(meas_inds,:);
        V = V(meas_inds,:);
        X_ref = X_ref(meas_inds,:);
        V_ref = V_ref(meas_inds,:);
        
        V_true_norm = norm(V,2);
        val_inds = setdiff(1:size(V_Ref,1),meas_inds); 
        
        
        if noisy_measurements == true && eta ~= 0
            %V = V + eta*abs(V).*randn(M,d);
            V = V.*(1+eta*randn(M,d));
            if strcmp(noise_type,'added') 
                
                dV = V*eta*randn(M,d);
                
                V = V + dV;
            end
            if strcmp(noise_type,'fixed')
                dV = eta*randn(M,d);
                V = V + dV;
            end
            dx1_err = norm(V(:,1)-V_Ref(meas_inds,1))/norm(V_Ref(meas_inds,1));
            dx2_err = norm(V(:,2)-V_Ref(meas_inds,2))/norm(V_Ref(meas_inds,2));
            dx3_err = norm(V(:,3)-V_Ref(meas_inds,3))/norm(V_Ref(meas_inds,3));
            FD_approx_error = mean([dx1_err, dx2_err, dx3_err])
        end
        V_noise_norm = norm(dV,2);
        
        
        L2_bound(r) = cond(D)*(D_noise_norm/D_true_norm + V_noise_norm/V_true_norm)
        
        if strcmp(precond_type,'approx')
            [W,bw]= opt_ksdensity(XX,D,'epanechnikov',[0.1 .3],M);
        end
        
        
        t0  = toc(timerval);      
        W_IT = eye(M); W_RR = eye(M);  W_SP = eye(M);  W_SPGL1 = eye(M); %W_LASSO = eye(M);
        weights = get_matrix_weights(D);
        Dbar = D*weights; 
        %*********** Solve for the precondtioned solutions
        Xi_WIT = zeros(size(Xi_ref)); Xi_WRR = zeros(size(Xi_ref)); Xi_WSPGL1 = zeros(size(Xi_ref));  Xi_WSP = zeros(size(Xi_ref));
        Xi_WIT_L = zeros(size(Xi_ref)); Xi_WRR_L = zeros(size(Xi_ref)); Xi_WSPGL1_L = zeros(size(Xi_ref));  Xi_WSP_L = zeros(size(Xi_ref));
        
        
        
        if strcmp(precond_type,'robust')

            %[Xi_WSPGL1_L, Xi_WSP_L, Xi_WRR_L, Xi_WIT_L, Sys_info] = adapt_SINDy(D,V,Lambda_cands)  
             [Xi_WSPGL1_L, Xi_WSP_L, Xi_WRR_L, Xi_WIT_L, Sys_info] = L_Adapt_SINDy(D,V,Lambda_cands,plotting,store)  
             Xi_WRR = psi_M\(psi_L*Xi_WRR_L);
             Xi_WSP = psi_M\(psi_L*Xi_WSP_L);
             Xi_WSPGL1 = psi_M\(psi_L*Xi_WSPGL1_L);
             Xi_WIT = psi_M\(psi_L*Xi_WIT_L);
%            Xi_WSP = pinv(psi_M)*(psi_L*Xi_WSP_L);
%            Xi_WSPGL1 = pinv(psi_M)*(psi_L*Xi_WSPGL1_L);
%            Xi_WIT = pinv(psi_M)*(psi_L*Xi_WIT_L);
%            Xi_WRR = pinv(psi_M)*(psi_L*Xi_WRR_L);

            Xi_RR_L = Sys_info{3}.reg_sol;
            Xi_SP_L = Sys_info{2}.reg_sol;
            Xi_SPGL1_L = Sys_info{1}.reg_sol;
            Xi_IT_L = Sys_info{4}.reg_sol;

            Xi_SP = psi_M\(psi_L*Xi_SP_L);
            Xi_SPGL1 = psi_M\(psi_L*Xi_SPGL1_L);
            Xi_IT = psi_M\(psi_L*Xi_IT_L);
            Xi_RR = psi_M\(psi_L*Xi_RR_L);
            %fprintf('r = %g, Adapt SINDy Solved in %g seconds \n',r, time);
        end
        
        
        


        
        %disp(Lambdas);
        G = Dbar'*Dbar;
        [MuD(r),MuD_av(r)] = coherence(D);
        normD(r) = norm(eye - G,'fro');
        inf_normD(r) = norm(Dbar(:),'inf');
        Cond_D(r) = cond(D);
        
        if strcmp(precond_type,'robust')    
            
            system_info =  structfun(@mean,Sys_info{1},'UniformOutput',false);
            Mu_WSPGL1(r) = system_info.Mu; 
            Mu_WSPGL1_av(r) = system_info.Mu_av; 
            inf_norm_WSPGL1(r) = system_info.inf_norm; 
            cond_WSPGL1D(r) = system_info.cond_WDbar; 
            norm_WSPGL1(r) = system_info.norm_WDbar;
            fro_norm_W_IT(r) = system_info.norm_W;
            cond_WSPGL1(r) = system_info.cond_W;
            inds =  find( imag(Sys_info{1}.lambda)  == 0)   
            lambda_WSPGL1(r) = get_avg_lambda(Sys_info{1}.lambda)
            
            system_info =  structfun(@mean,Sys_info{2},'UniformOutput',false);
            Mu_WSP(r) = system_info.Mu; 
            Mu_WSP_av(r) = system_info.Mu_av; 
            inf_norm_WSP(r) = system_info.inf_norm; 
            cond_WSPD(r) = system_info.cond_WDbar; 
            norm_WSP(r) = system_info.norm_WDbar;
            fro_norm_W_IT(r) = system_info.norm_W;
            cond_WSP(r) = system_info.cond_W;
            lambda_WSP(r) = get_avg_lambda(Sys_info{2}.lambda)
            
            system_info =  structfun(@mean,Sys_info{3},'UniformOutput',false);
            Mu_WRR(r) = system_info.Mu; 
            Mu_WRR_av(r) = system_info.Mu_av; 
            inf_norm_WRR(r) = system_info.inf_norm; 
            cond_WRRD(r) = system_info.cond_WDbar; 
            norm_WRR(r) = system_info.norm_WDbar;
            fro_norm_W_IT(r) = system_info.norm_W;
            cond_WRR(r) = system_info.cond_W;
            lambda_WRR(r) = get_avg_lambda(Sys_info{3}.lambda)
            
            system_info =  structfun(@mean,Sys_info{4},'UniformOutput',false);
            Mu_WIT(r) = system_info.Mu; 
            Mu_WIT_av(r) = system_info.Mu_av; 
            inf_norm_WIT(r) = system_info.inf_norm; 
            cond_WITD(r) = system_info.cond_WDbar; 
            norm_WIT(r) = system_info.norm_WDbar;
            fro_norm_W_IT(r) = system_info.norm_W;
            cond_WIT(r) = system_info.cond_W;
            lambda_WIT(r) = get_avg_lambda(Sys_info{4}.lambda)
            
            SPGL1_lambdas_iter = Sys_info{1}.lambda;
            fprintf('SPGL1 selected lambdas: \n')
            disp(Sys_info{1}.lambda)
            fprintf('SP selected lambdas: \n')
            disp(Sys_info{2}.lambda)
            fprintf('RR selected lambdas: \n')
            disp(Sys_info{3}.lambda)
            fprintf('IT selected lambdas: \n')
            disp(Sys_info{4}.lambda)
        end
        
        if strcmp(precond_type,'robust') == 0

            %*********** Solve for the SINDy solution
            [Xi_SPGL1_L, Xi_SP_L, Xi_RR_L, Xi_IT_L, Sys_info] = adapt_SINDy(D,V,[],store);
            Xi_SP = psi_M\(psi_L*Xi_SP_L);
            Xi_SPGL1 = psi_M\(psi_L*Xi_SPGL1_L);
            Xi_IT = psi_M\(psi_L*Xi_IT_L);
            Xi_RR = psi_M\(psi_L*Xi_RR_L);
            Xi_WSP = Xi_SP; Xi_WSPGL1 = Xi_SPGL1; Xi_WIT = Xi_IT; Xi_WRR = Xi_RR; % Xi_WLASSO = Xi_LASSO; 

            %fprintf('r = %g, SINDy Solved in %g seconds \n',r, time);
 
        end
        
        % Compute the validation errors

        WRR_fit(r)  = mean([norm(D*Xi_WRR_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_WRR_L(:,2) - V(:,2))/norm(V(:,2)), ...
            norm(D*Xi_WRR_L(:,3) - V(:,3))/norm(V(:,3))]);
        RR_fit(r)  = mean([norm(D*Xi_RR_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_RR_L(:,2) - V(:,2))/norm(V(:,2)), ...
            norm(D*Xi_RR_L(:,3) - V(:,3))/norm(V(:,3))]);
        IT_fit(r)  = mean([norm(D*Xi_IT_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_IT_L(:,2) - V(:,2))/norm(V(:,2)), ...
            norm(D*Xi_IT_L(:,3) - V(:,3))/norm(V(:,3))]);
        WIT_fit(r)  = mean([norm(D*Xi_WIT_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_WIT_L(:,2) - V(:,2))/norm(V(:,2)), ...
            norm(D*Xi_WIT_L(:,3) - V(:,3))/norm(V(:,3))]);
        SP_fit(r)  = mean([norm(D*Xi_SP_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_SP_L(:,2) - V(:,2))/norm(V(:,2))...
            norm(D*Xi_SP_L(:,3) - V(:,3))/norm(V(:,3))]);
        WSP_fit(r) = mean([norm(D*Xi_WSP_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_WSP_L(:,2) - V(:,2))/norm(V(:,2))...
            norm(D*Xi_WSP_L(:,3) - V(:,3))/norm(V(:,3))]);
        SPGL1_fit(r)  = mean([norm(D*Xi_SPGL1_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_SPGL1_L(:,2) - V(:,2))/norm(V(:,2))...
            norm(D*Xi_SPGL1_L(:,3) - V(:,3))/norm(V(:,3))]);
        WSPGL1_fit(r) = mean([norm(D*Xi_WSPGL1_L(:,1) - V(:,1))/norm(V(:,1)),...
            norm(D*Xi_WSPGL1_L(:,2) - V(:,2))/norm(V(:,2))...
            norm(D*Xi_WSPGL1_L(:,3) - V(:,3))/norm(V(:,3))]);

        WRR_err(r)  = mean([norm(D_val*Xi_WRR(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_WRR(:,2) - V_val(:,2))/norm(V_val(:,2)), ...
            norm(D_val*Xi_WRR(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        RR_err(r)  = mean([norm(D_val*Xi_RR(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_RR(:,2) - V_val(:,2))/norm(V_val(:,2)), ...
            norm(D_val*Xi_RR(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        IT_err(r)  = mean([norm(D_val*Xi_IT(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_IT(:,2) - V_val(:,2))/norm(V_val(:,2)), ...
            norm(D_val*Xi_IT(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        WIT_err(r)  = mean([norm(D_val*Xi_WIT(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_WIT(:,2) - V_val(:,2))/norm(V_val(:,2)), ...
            norm(D_val*Xi_WIT(:,3) - V_val(:,3))/norm(V_val(:,3))]);
%         LASSO_err(r)  = mean([norm(D_val*Xi_LASSO(:,1) - V_val(:,1))/norm(V_val(:,1)),...
%             norm(D_val*Xi_LASSO(:,2) - V_val(:,2))/norm(V_val(:,2))...
%             norm(D_val*Xi_LASSO(:,3) - V_val(:,3))/norm(V_val(:,3))]);
%         WLASSO_err(r) = mean([norm(D_val*Xi_WLASSO(:,1) - V_val(:,1))/norm(V_val(:,1)),...
%             norm(D_val*Xi_WLASSO(:,2) - V_val(:,2))/norm(V_val(:,2))...
%             norm(D_val*Xi_WLASSO(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        SP_err(r)  = mean([norm(D_val*Xi_SP(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_SP(:,2) - V_val(:,2))/norm(V_val(:,2))...
            norm(D_val*Xi_SP(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        WSP_err(r) = mean([norm(D_val*Xi_WSP(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_WSP(:,2) - V_val(:,2))/norm(V_val(:,2))...
            norm(D_val*Xi_WSP(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        SPGL1_err(r)  = mean([norm(D_val*Xi_SPGL1(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_SPGL1(:,2) - V_val(:,2))/norm(V_val(:,2))...
            norm(D_val*Xi_SPGL1(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        WSPGL1_err(r) = mean([norm(D_val*Xi_WSPGL1(:,1) - V_val(:,1))/norm(V_val(:,1)),...
            norm(D_val*Xi_WSPGL1(:,2) - V_val(:,2))/norm(V_val(:,2))...
            norm(D_val*Xi_WSPGL1(:,3) - V_val(:,3))/norm(V_val(:,3))]);
        
        WRR_model_err(r) = norm(Xi_WRR(:) - Xi_ref(:))/norm(Xi_ref(:));
        RR_model_err(r) = norm(Xi_RR(:) - Xi_ref(:))/norm(Xi_ref(:));
        IT_model_err(r) = norm(Xi_IT(:) - Xi_ref(:))/norm(Xi_ref(:));
        WIT_model_err(r) = norm(Xi_WIT(:) - Xi_ref(:))/norm(Xi_ref(:));
%         LASSO_model_err(r) = norm(Xi_LASSO(:) - Xi_ref(:))/norm(Xi_ref(:));
%         WLASSO_model_err(r) = norm(Xi_WLASSO(:) - Xi_ref(:))/norm(Xi_ref(:));
        SP_model_err(r) = norm(Xi_SP(:) - Xi_ref(:))/norm(Xi_ref(:));
        WSP_model_err(r) = norm(Xi_WSP(:) - Xi_ref(:))/norm(Xi_ref(:)); 
        SPGL1_model_err(r) = norm(Xi_SPGL1(:) - Xi_ref(:))/norm(Xi_ref(:));
        WSPGL1_model_err(r) = norm(Xi_WSPGL1(:) - Xi_ref(:))/norm(Xi_ref(:)); 
        if true 
%             fprintf('LASSO_err = %f \n',LASSO_model_err(r));
%             fprintf('WLASSO_err = %f \n',WLASSO_model_err(r));
            fprintf('IT_err = %f   , Val_err = %f \n',IT_model_err(r), IT_err(r));
            fprintf('WIT_err = %f   , Val_err = %f \n',WIT_model_err(r), WIT_err(r));
            fprintf('RR_err = %f   , Val_err = %f \n',RR_model_err(r), RR_err(r));
            fprintf('WRR_err = %f   , Val_err = %f \n',WRR_model_err(r), WRR_err(r));
            fprintf('SPGL1_err = %f   , Val_err = %f \n',SPGL1_model_err(r), SPGL1_err(r));
            fprintf('WSPGL1_err = %f   , Val_err = %f \n',WSPGL1_model_err(r), WSPGL1_err(r));
            fprintf('SP_err = %f   , Val_err = %f \n',SP_model_err(r), SP_err(r));
            fprintf('WSP_err = %f   , Val_err = %f \n',WSP_model_err(r), WSP_err(r));
        end
  
         %disp(Lambdas)
    mean_MuD = mean(MuD);mean_MuD_av = mean(MuD_av); mean_inf_normD = mean(inf_normD); mean_Cond_D = mean(Cond_D); mean_normD = mean(normD);
    
    RR_stats(n,:) = [M,mean(RR_model_err),mean(RR_err),mean(RR_fit),...
                    mean_MuD,mean_MuD_av,mean_inf_normD,mean_Cond_D,mean_normD,...
                    std(RR_model_err),std(RR_err),std(RR_fit)]; 
    IT_stats(n,:) = [M,mean(IT_model_err),mean(IT_err),mean(IT_fit),...
                    mean_MuD,mean_MuD_av,mean_inf_normD,mean_Cond_D,mean_normD...
                    std(IT_model_err),std(IT_err),std(IT_fit)]; 
    SPGL1_stats(n,:) = [M,mean(SPGL1_model_err),mean(SPGL1_err),mean(SPGL1_fit),...
                    mean_MuD,mean_MuD_av,mean_inf_normD,mean_Cond_D,mean_normD,...
                    std(SPGL1_model_err),std(SPGL1_err),std(SPGL1_fit)]; 
    SP_stats(n,:) = [M,mean(SP_model_err),mean(SP_err),mean(SP_fit),...
                    mean_MuD,mean_MuD_av,mean_inf_normD,mean_Cond_D,mean_normD,...
                    std(SP_model_err),std(SP_err),std(SP_fit)]; 
    
    mean_Mu_WIT = mean(Mu_WIT); mean_Mu_WIT_av = mean(Mu_WIT_av); mean_inf_norm_WIT = mean(inf_norm_WIT); mean_Cond_WITD = mean(cond_WITD); mean_norm_WIT = mean(norm_WIT); mean_fro_norm_W_IT = mean(fro_norm_W_IT); mean_Lambda_WIT = mean(lambda_WIT); mean_cond_WIT = mean(cond_WIT);
    WIT_stats(n,:) = [M,mean(WIT_model_err),mean(WIT_err),mean(WIT_fit)...
                    mean_Mu_WIT,mean_Mu_WIT_av,mean_inf_norm_WIT,mean_Cond_WITD,...
                    mean_norm_WIT,mean_cond_WIT,mean_fro_norm_W_IT,mean_Lambda_WIT,...
                    std(WIT_model_err),std(WIT_err),std(WIT_fit)]; 
    mean_Mu_WRR = mean(Mu_WRR); mean_Mu_WRR_av = mean(Mu_WRR_av); mean_inf_norm_WRR = mean(inf_norm_WRR); mean_Cond_WRRD = mean(cond_WRRD); mean_norm_WRR = mean(norm_WRR); mean_fro_norm_W_RR = mean(fro_norm_W_RR); mean_Lambda_WRR = mean(lambda_WRR); mean_cond_WRR = mean(cond_WRR);
    WRR_stats(n,:) = [M,mean(WRR_model_err),mean(WRR_err),mean(WRR_fit)...
                    mean_Mu_WRR,mean_Mu_WRR_av,mean_inf_norm_WRR,mean_Cond_WRRD,...
                    mean_norm_WRR,mean_cond_WRR,mean_fro_norm_W_RR,mean_Lambda_WRR,...
                    std(WRR_model_err),std(WRR_err),std(WRR_fit)];   
    mean_Mu_WSP = mean(Mu_WSP); mean_Mu_WSP_av = mean(Mu_WSP_av); mean_inf_norm_WSP = mean(inf_norm_WSP); mean_Cond_WSPD = mean(cond_WSPD); mean_norm_WSP = mean(norm_WSP); mean_fro_norm_W_SP = mean(fro_norm_W_SP); mean_Lambda_WSP = mean(lambda_WSP); mean_cond_WSP = mean(cond_WSP);
    WSP_stats(n,:) = [M,mean(WSP_model_err),mean(WSP_err),mean(WSP_fit)...
                    mean_Mu_WSP,mean_Mu_WSP_av,mean_inf_norm_WSP,mean_Cond_WSPD,...
                    mean_norm_WSP,mean_cond_WSP,mean_fro_norm_W_SP,mean_Lambda_WSP,...
                    std(WSP_model_err),std(WSP_err),std(WSP_fit)]; 
                    
    mean_Mu_WSPGL1 = mean(Mu_WSPGL1); mean_Mu_WSPGL1_av = mean(Mu_WSPGL1_av); mean_inf_norm_WSPGL1 = mean(inf_norm_WSPGL1); mean_Cond_WSPGL1D = mean(cond_WSPGL1D); mean_norm_WSPGL1 = mean(norm_WSPGL1); mean_fro_norm_W_SPGL1 = mean(fro_norm_W_SPGL1); mean_Lambda_WSPGL1 = mean(lambda_WSPGL1); mean_cond_WSPGL1 = mean(cond_WSPGL1);
    WSPGL1_stats(n,:) = [M,mean(WSPGL1_model_err),mean(WSPGL1_err),mean(WSPGL1_fit)...
                    mean_Mu_WSPGL1,mean_Mu_WSPGL1_av,mean_inf_norm_WSPGL1,mean_Cond_WSPGL1D,...
                    mean_norm_WSPGL1,mean_cond_WSPGL1,mean_fro_norm_W_SPGL1,mean_Lambda_WSPGL1,...
                    std(WSPGL1_model_err),std(WSPGL1_err),std(WSPGL1_fit)]; 
                
   
    RR = table(RR_stats(1:n,1),RR_stats(1:n,2),RR_stats(1:n,3),RR_stats(1:n,4),...
               RR_stats(1:n,5),RR_stats(1:n,6),RR_stats(1:n,7),RR_stats(1:n,8),RR_stats(1:n,9),...
               RR_stats(1:n,10),RR_stats(1:n,11),RR_stats(1:n,12),...
               'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu' ,'Mu_avg','max_entry','Cond_D','Dist','Ref_std','Val_std','Fit_std'})
    WRR = table(WRR_stats(1:n,1),WRR_stats(1:n,2),WRR_stats(1:n,3),WRR_stats(1:n,4),...
                WRR_stats(1:n,5),WRR_stats(1:n,6),WRR_stats(1:n,7),WRR_stats(1:n,8),...
                WRR_stats(1:n,9),WRR_stats(1:n,10),WRR_stats(1:n,11),WRR_stats(1:n,12),...
                WRR_stats(1:n,13),WRR_stats(1:n,14),WRR_stats(1:n,15),...
           'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu','Mu_avg','max_entry','Cond_WD','Dist','Cond_W','Dist_W','lambda','Ref_std','Val_std','Fit_std'})
    IT = table(IT_stats(1:n,1),IT_stats(1:n,2),IT_stats(1:n,3),IT_stats(1:n,4),...
               IT_stats(1:n,5),IT_stats(1:n,6),IT_stats(1:n,7),IT_stats(1:n,8),IT_stats(1:n,9),...
               IT_stats(1:n,10),IT_stats(1:n,11),IT_stats(1:n,12),...
               'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu' ,'Mu_avg','max_entry','Cond_D','Dist','Ref_std','Val_std','Fit_std'})
    WIT = table(WIT_stats(1:n,1),WIT_stats(1:n,2),WIT_stats(1:n,3),WIT_stats(1:n,4),...
                WIT_stats(1:n,5),WIT_stats(1:n,6),WIT_stats(1:n,7),WIT_stats(1:n,8),...
                WIT_stats(1:n,9),WIT_stats(1:n,10),WIT_stats(1:n,11),WIT_stats(1:n,12),...
                WIT_stats(1:n,13),WIT_stats(1:n,14),WIT_stats(1:n,15),...
                'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu','Mu_avg','max_entry','Cond_WD','Dist','Cond_W','Dist_W','lambda','Ref_std','Val_std','Fit_std'}) 
    SP = table(SP_stats(1:n,1),SP_stats(1:n,2),SP_stats(1:n,3),SP_stats(1:n,4),...
                SP_stats(1:n,5),SP_stats(1:n,6),SP_stats(1:n,7),SP_stats(1:n,8),SP_stats(1:n,9),...
                SP_stats(1:n,10),SP_stats(1:n,11),SP_stats(1:n,12),...
                'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu' ,'Mu_avg','max_entry','Cond_D','Dist','Ref_std','Val_std','Fit_std'})
    WSP = table(WSP_stats(1:n,1),WSP_stats(1:n,2),WSP_stats(1:n,3),WSP_stats(1:n,4),...
                WSP_stats(1:n,5),WSP_stats(1:n,6),WSP_stats(1:n,7),WSP_stats(1:n,8),...
                WSP_stats(1:n,9),WSP_stats(1:n,10),WSP_stats(1:n,11),WSP_stats(1:n,12),...
                WSP_stats(1:n,13),WSP_stats(1:n,14),WSP_stats(1:n,15),...
                'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu','Mu_avg','max_entry','Cond_WD','Dist','Cond_W','Dist_W','lambda','Ref_std','Val_std','Fit_std'})
    SPGL1 = table(SPGL1_stats(1:n,1),SPGL1_stats(1:n,2),SPGL1_stats(1:n,3),SPGL1_stats(1:n,4),...
                SPGL1_stats(1:n,5),SPGL1_stats(1:n,6),SPGL1_stats(1:n,7),SPGL1_stats(1:n,8),SPGL1_stats(1:n,9),...
                SPGL1_stats(1:n,10),SPGL1_stats(1:n,11),SPGL1_stats(1:n,12),...
                'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu' ,'Mu_avg','max_entry','Cond_D','Dist','Ref_std','Val_std','Fit_std'})
    WSPGL1 = table(WSPGL1_stats(1:n,1),WSPGL1_stats(1:n,2),WSPGL1_stats(1:n,3),WSPGL1_stats(1:n,4),...
                WSPGL1_stats(1:n,5),WSPGL1_stats(1:n,6),WSPGL1_stats(1:n,7),WSPGL1_stats(1:n,8),...
                WSPGL1_stats(1:n,9),WSPGL1_stats(1:n,10),WSPGL1_stats(1:n,11),WSPGL1_stats(1:n,12),...
                WSPGL1_stats(1:n,13),WSPGL1_stats(1:n,14),WSPGL1_stats(1:n,15),...
                'VariableNames',{'N','Ref_err','Val_err','Fit_rel_err','Mu','Mu_avg','max_entry','Cond_WD','Dist','Cond_W','Dist_W','lambda','Ref_std','Val_std','Fit_std'})

      if max(XX(:)) > 1 || min(XX(:)) < -1
            fprintf('Warning: Data lies outside [-1,1] support! \n')
      end
%     fprintf([save_name,'\n'])
%     save(save_name)
%     time = toc(timerval)
    
    % Getting Matlab memory usage from Linux/Unix/Mac 
    % get the parent process id
    [s,ppid] = unix(['ps -p $PPID -l | ' awkCol('PPID') ]);
    % get memory used by the parent process (resident set size)
    [s,memory_used] = unix(['ps -O rss -p ' strtrim(ppid) ' | awk ''NR>1 {print$2}'' ']); 
    % rss is in kB, convert to bytes 
    fprintf('Memory used on ppid %i for N-iteration:  %g (MB) \n',str2double(ppid),str2double(memory_used)*10^(-3));
    

    

    %save('lorenz63_benchmark_data.mat','D','V','D_M','D_L','X_ref','V_ref','Xi_ref','Xi_SPGL1','psi_M','psi_L','FD_approx_error','eta')
    
    
%%
   save_name = ['120_',func2str(func),'_state_data.mat',]
   save(save_name)
    
%%
close all
%clear all
%load('lorenz63_movie_data.mat')

T = [T0,2*TF]; 
dt = 0.00005;
T_val = [0:dt:TF];


%X0 = [-0.1,0.1,0.2]';

X0 = X0_ref_scaled;

[TT,X_data_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_ref);

    
[TT,X_data,V_ref,D_M,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_ref);
dTA = [0; diff(TT)];
%dta = zeros(size(TT))
%dTA = sqrt( sum((X_data_val - X_data_val).^2,2));



%F1 = figure
color_line3(X_data(:,1),X_data(:,2),X_data(:,3),dTA,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',28) 
ylabel('$y$','interpreter','latex','fontsize',28) 
zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',28)
title('Exact Solution','interpreter','latex','fontsize',28);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',28)

set(gca,'FontSize',28)
length(TT)
view([-45,45,45])
%%
save_name = ['Exact_T8_',func2str(func)]
save_name
%print(save_name,'-dpng','-r300')
saveas(F1,save_name,'epsc')


%%
close all
%F1 = figure

if strcmp(precond_type,'robust')
%T = [T0,TF]; 

[~,X_WSPGL1_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_WSPGL1);
WSPGL1_validation_error = norm(X_data_val-X_WSPGL1_val)/norm(X_data_val);
[TT_WSPGL1,X_WSPGL1,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_WSPGL1);


WSPGL1_Lyap = get_Lyap_exps(X_data_val,X_WSPGL1_val);

%dTA_WSPGL1 = sqrt( sum((X_WSPGL1_val - X_data_val).^2,2))';
dTA_WSPGL1 = [0; diff(TT_WSPGL1)];


subplot(2,4,1)
color_line3(X_WSPGL1(:,1),X_WSPGL1(:,2),X_WSPGL1(:,3),dTA_WSPGL1,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)

title_string = [{'BPDN$_{RPMO}$'}, {['$e_{ref}=$',num2str(WSPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(WSPGL1_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_WSPGL1)
view([-45,45,45])


[~,X_WRR_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_WRR);
WRR_validation_error = norm(X_data_val-X_WRR_val)/norm(X_data_val);

[TT_WRR,X_WRR,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_WRR);
%dTA_WRR = sqrt( sum((X_WRR_val - X_data_val).^2,2))';
dTA_WRR = [0; diff(TT_WRR)];

WRR_Lyap = get_Lyap_exps(X_data_val,X_WRR_val);

subplot(2,4,2)
color_line3(X_WRR(:,1),X_WRR(:,2),X_WRR(:,3),dTA_WRR,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STRidge$_{RPMO}$'}, {['$e_{ref}=$',num2str(WRR_stats(1:n,2))]},{['$e_{val}=$',num2str(WRR_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_WRR)
view([-45,45,45])




%T = [0,TF]; 
[~,X_WIT_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_WIT);
WIT_validation_error = norm(X_data_val-X_WIT_val)/norm(X_data_val);


[TT_WIT,X_WIT,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_WIT);
%dTA_WIT = sqrt( sum((X_WIT_val - X_data_val).^2,2))';
dTA_WIT = [0; diff(TT_WIT)];
WIT_Lyap = get_Lyap_exps(X_data_val,X_WIT_val);

subplot(2,4,3)
color_line3(X_WIT(:,1),X_WIT(:,2),X_WIT(:,3),dTA_WIT,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STLS$_{RPMO}$'}, {['$e_{ref}=$',num2str(WIT_stats(1:n,2))]}, {['$e_{val}=$',num2str(WIT_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_WIT)
view([-45,45,45])

%T = [0,TF]; 
[~,X_WSP_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_WSP);
WSP_validation_error = norm(X_data_val-X_WSP_val)/norm(X_data_val);
[TT_WSP,X_WSP,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_WSP);
dTA_WSP = [0; diff(TT_WSP)];
%dTA_WSP = sqrt( sum((X_WSP_val - X_data_val).^2,2))';

WSP_Lyap = get_Lyap_exps(X_data_val,X_WSP_val);

subplot(2,4,4)
color_line3(X_WSP(:,1),X_WSP(:,2),X_WSP(:,3),dTA_WSP,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'SP$_{RPMO}$'}, {['$e_{ref}=$',num2str(WSP_stats(1:n,2))]},{['$e_{val}=$',num2str(WSP_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_WSP)
view([-45,45,45])




%T = [0,TF]; 
[opt_T,X_SPGL1_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_SPGL1);
[TT_SPGL1,X_SPGL1,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_SPGL1);
%dTA_SPGL1 = sqrt( sum((X_SPGL1_val - X_data_val).^2,2))';
dTA_SPGL1 = [0; diff(TT_SPGL1)];
try
    SPGL1_validation_error = norm(X_data_val-X_SPGL1_val)/norm(X_data_val);
    SPGL1_Lyap = get_Lyap_exps(X_data_val,X_WRR_val);
catch ME
    fprintf( [ME.message ' Setting validation error to Inf \n']) 
    SPGL1_validation_error = Inf;
    SPGL1_Lyap = get_Lyap_exps(X_data_val,X_data_val);
  
end



subplot(2,4,5)
color_line3(X_SPGL1(:,1),X_SPGL1(:,2),X_SPGL1(:,3),dTA_SPGL1,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'BPDN'}, {['$e_{ref}=$',num2str(SPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(SPGL1_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_SPGL1)
view([-45,45,45])

%T = [0,TF]; 
[~,X_RR_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_RR);
RR_validation_error = norm(X_data_val-X_RR_val)/norm(X_data_val);
[TT_RR,X_RR,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_RR);

RR_Lyap = get_Lyap_exps(X_data_val,X_RR_val);
dTA_RR = [0; diff(TT_RR)];

subplot(2,4,6)
color_line3(X_RR(:,1),X_RR(:,2),X_RR(:,3),dTA_RR,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STRidge'}, {['$e_{ref}=$',num2str(RR_stats(1:n,2))]}, {['$e_{val}=$',num2str(RR_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_RR)
view([-45,45,45])

[~,X_IT_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_IT);
IT_validation_error = norm(X_data_val-X_IT_val)/norm(X_data_val);


[TT_IT,X_IT,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_IT);

IT_Lyap = get_Lyap_exps(X_data_val,X_IT_val);
dTA_IT = [0; diff(TT_IT)];

subplot(2,4,7)
color_line3(X_IT(:,1),X_IT(:,2),X_IT(:,3),dTA_IT,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STLS'}, {['$e_{ref}=$',num2str(IT_stats(1:n,2))]}, {['$e_{val}=$',num2str(IT_validation_error)]} ];
title(title_string,'interpreter','latex','fontsize',24);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_IT)
view([-45,45,45])


%T = [0,TF]; 
[~,X_SP_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_SP);
SP_validation_error = norm(X_data_val-X_SP_val)/norm(X_data_val);
[TT_SP,X_SP,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_SP);
SP_Lyap = get_Lyap_exps(X_data_val,X_SP_val);
dTA_SP = [0; diff(TT_SP)];
%dTA_SP = sqrt( sum((X_SP_val - X_data_val).^2,2))';

subplot(2,4,8)
color_line3(X_SP(:,1),X_SP(:,2),X_SP(:,3),dTA_SP,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'SP'}, {['$e_{ref}=$',num2str(SP_stats(1:n,2))]}, {['$e_{val}=$',num2str(SP_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_SP)
view([-45,45,45])

else 
     

[~,X_SPGL1_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_SPGL1);
SPGL1_validation_error = norm(X_data_val-X_SPGL1_val)/norm(X_data_val);    
[TT_SPGL1,X_SPGL1,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_SPGL1);


SPGL1_Lyap = get_Lyap_exps(X_data_val,X_SPGL1_val);
dTA_SPGL1 = [0; diff(TT_SPGL1)];

%dTA_SPGL1 = sqrt( sum((X_SPGL1_val - X_data_val).^2,2))';


subplot(1,4,1)
color_line3(X_SPGL1(:,1),X_SPGL1(:,2),X_SPGL1(:,3),dTA_SPGL1,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'BPDN'}, {['$e_{ref}=$',num2str(SPGL1_stats(1:n,2))]}, {['$e_{val}=$',num2str(SPGL1_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_SPGL1)
view([-45,45,45])

%T = [0,TF]; 
[~,X_RR_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_RR);
RR_validation_error = norm(X_data_val-X_RR_val)/norm(X_data_val);
  
[TT_RR,X_RR,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_RR);



RR_Lyap = get_Lyap_exps(X_data_val,X_RR_val);
dTA_RR = [0; diff(RR)];
%dTA_RR = sqrt( sum((X_RR_val - X_data_val).^2,2))';


subplot(1,4,2)
color_line3(X_RR(:,1),X_RR(:,2),X_RR(:,3),dTA_RR,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STRidge'}, {['$e_{ref}=$',num2str(RR_stats(1:n,2))]}, {['$e_{val}=$',num2str(RR_validation_error)]} ];
title(title_string,'interpreter','latex','fontsize',24);

set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_RR)
view([-45,45,45])


%T = [0,TF]; 
[~,X_IT_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_IT);
IT_validation_error = norm(X_data_val-X_IT_val)/norm(X_data_val);
  
[TT_IT,X_IT,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_IT); 


IT_Lyap = get_Lyap_exps(X_data_val,X_IT_val);
dTA_IT = [0; diff(TT_IT)];
%dTA_IT = sqrt( sum((X_IT_val - X_data_val).^2,2))';


subplot(1,4,3)
color_line3(X_IT(:,1),X_IT(:,2),X_IT(:,3),dTA_IT,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'STLS'}, {['$e_{ref}=$',num2str(IT_stats(1:n,2))]}, {['$e_{val}=$',num2str(IT_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)

set(gca,'FontSize',28)
length(TT_IT)
view([-45,45,45])


[~,X_SP_val,~,~,~,~,~,~] = solve_ode(X0,T_val,@SINDy_func,index_pc,Xi_SP);
SP_validation_error = norm(X_data_val-X_SP_val)/norm(X_data_val);
[TT_SP,X_SP,~,~,~,~,~,~] = solve_ode(X0,T,@SINDy_func,index_pc,Xi_SP);


SP_Lyap = get_Lyap_exps(X_data_val,X_SP_val);
%dTA_SP = sqrt( sum((X_SP_val - X_data_val).^2,2))';
dTA_SP = [0; diff(TT_SP)];


subplot(1,4,4)
color_line3(X_SP(:,1),X_SP(:,2),X_SP(:,3),dTA_SP,'LineWidth',4);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)

ylabel('$y$','interpreter','latex','fontsize',24)

zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
title_string = [{'SP'}, {['$e_{ref}=$',num2str(SP_stats(1:n,2))]},  {['$e_{val}=$',num2str(SP_validation_error)]}];
title(title_string,'interpreter','latex','fontsize',24);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)


set(gca,'FontSize',28)
length(TT_SP)
view([-45,45,45])


end

%%
save_name = ['recon_500_',func2str(func)]
save_name
%print(save_name,'-dpng','-r300')
saveas(F1,save_name,'epsc')



%%

save_name = ['120_',func2str(func),'_movie_data.mat']
save(save_name)
%%
save_name = ['200_',func2str(func)]
print(save_name,'-dpng','-r300')

%%
save_name = ['800_prediction_',func2str(func),'_M',num2str(M)]

save_name
print(save_name,'-dpng','-r300') 



%%
close all
F1 = figure
plot3(XX(:,1),XX(:,2),XX(:,3), 'bx','markersize',16);
hold on
view(27,16)
axis equal;
grid;
xlabel('$x$','interpreter','latex','fontsize',24)
ylabel('$y$','interpreter','latex','fontsize',24)
zlabel('$z$ \hspace{0.5cm}','interpreter','latex','fontsize',24)
%title_string = [{'SP'}, {['$e_{ref}=$',num2str(SP_stats(1:n,2))]} ,{['$e_{val}=$',num2str(SP_validation_error) ]}];
%title('Data $\mathbf{X}$','interpreter','latex','fontsize',24);
set(get(gca,'ZLabel'),'Rotation',pi/2,'FontSize',20)
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
set(gca,'FontSize',28)
%length(TT)
view([-45,45,45])

%%

save_name = ['Data_large_',func2str(func)]
save_name
%print(save_name,'-dpng','-r300')
saveas(F1,save_name,'epsc')

%%
close all
figure
subplot(3,2,1)
plot(1:M,XX(:,1),'bx');
subplot(3,2,3)
plot(1:M,XX(:,2),'bx')
subplot(3,2,5)
plot(1:M,XX(:,3),'bx')
subplot(3,2,2)
plot(1:M,V(:,1),'LineWidth',4);
subplot(3,2,4)
plot(1:M,V(:,2),'LineWidth',4);
subplot(3,2,6)
plot(1:M,V(:,3),'LineWidth',4);
set(gca,'FontSize',28)
