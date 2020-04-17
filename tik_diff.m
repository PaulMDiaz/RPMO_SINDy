function [dx ] = tik_diff(x,T,lambda , plotting,lambda_cands)
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');
[G,n,dt] = diff_mat(T);
d = x-x(1);


D3 = -1*diag(ones(n-1,1),1) + 0.5*diag(ones(n-2,1),2) + diag(ones(n-1,1),-1) -  0.5*diag(ones(n-2,1),-2);

D2 = full(gallery('tridiag',n,-1,2,-1))/(dt^2);
D1 = full(gallery('tridiag',n,0,-1,1))/dt;


%D1(end,end-1) = 1;
L = D1'*D2';
%L = D1'*D2'*D3';
%L = D1;

%G

%L = [D1', D2'];
%L = [eye(n); D1; D2;];
%L = D2;
%pinv_mat = geninv(A'*A + lambda*(Dt*Dt'));
%dx = pinv_mat*(A'*y);
%dx = pinv(A'*A + lambda*(Dt*Dt'))*(A'*y);

if nargin < 4
    plotting = false;
end

% size(G)
% size(L)

if length(lambda) == 1
    [U,sm, Y,v] =  cgsvd(G,L);
    p = floor(length(sm/2));
    Lam = diag(sm(1:p));
    M = diag(sm(p+1:end));
    x = (Lam'*Lam + lambda^2*M'*M)\(Lam'*U'*d);
    dx = Y*x;
end

if isempty(lambda) || nargin < 3
    dx = L_curve_lambda(G,L,d,plotting,lambda_cands);
end



end


function [A,n,dt] = diff_mat(T)
n = length(T);
dt = T(2)-T(1);
a = T(1);
b = T(end);
A = zeros(n);
t = zeros(n,1);

for j = 1:n+1
    t(j) = a +(j-1)*dt;
end

for j = 1:n
    for i = 1:n
        if T(i) <= t(j)
            A(i,j) = 0;
        end
        if (t(j) < T(i)) && ( T(i)  < t(j+1))
           A(i,j) = T(i)-t(j); 
        end
        if t(j+1) <= T(i)
           A(i,j) = dt; 
        end
        
    end
    
end
end

function dx = L_curve_lambda(G,L,d,plotting,lambda_cands)
 %inc 20
%lambda_cands =  10.^(linspace(-6,-5,100)); %[ 10.^(linspace())];
%lambda_cands =  10.^(linspace(-6,-5,1000)); %[ 10.^(linspace())];
%lambda_cands = linspace( 6e-6,2e-6,100); %MS for eta 0.01 sample inc 20
%lambda_cands = linspace( 4.9e-6,2e-6,100); 


% ambda_cands = linspace( 6e-6,1e-8,300);


[U,sm, Y,~,~] =  cgsvd(G,L);
p = floor(length(sm/2));
Lam = diag(sm(1:p));
M = diag(sm(p+1:end));
dxdt = []; semi_norm = []; residual = [];
for i = 1 : length(lambda_cands)
    lambda = lambda_cands(i);
    x = (Lam'*Lam + lambda^2*M'*M)\(Lam'*U'*d);
    dxdt(:,i) = Y*x;
    semi_norm(i) = norm(L*dxdt(:,i) ,2);
    residual(i) = norm(G*dxdt(:,i) -d,2);
end

residual = residual/norm(residual);
semi_norm = semi_norm/norm(semi_norm);
[semi_norm , I] = sort(semi_norm);
residual = residual(I);
dxdt = dxdt(:,I);
lambda_cands = lambda_cands(I);
    


[~,min_ind] = min( sqrt((min(residual)- residual).^2 + (min(semi_norm) -semi_norm).^2));
%My_mat = [semi_norm', residual', (sqrt((min(residual)- residual).^2 + (min(semi_norm) -semi_norm).^2))' ]
%My_mat(min_ind,:)
%lambda = lambda_cands(min_ind)


if plotting
    figure
    plot(semi_norm,residual,'-o',[min(semi_norm),semi_norm(min_ind)],[min(residual),residual(min_ind)],'-rs');
    ylabel('$||Gm -d ||_2$','interpreter','latex')
    xlabel('$|| Lm ||_2$','interpreter','latex')
    title('Tikhonov Regularization','interpreter','latex')
    ylim([0.9*min(min(semi_norm),min(residual)), 1.1*max(max(semi_norm),max(residual))]);
    ylim([0.9*min(min(semi_norm),min(residual)), 1.1*max(max(semi_norm),max(residual))]);
    %semi_norm
    %residual
    fprintf('tik diff lambda = %g  \n', lambda_cands(min_ind))
    pause 
    close all
end

dx = dxdt(:,min_ind);

end