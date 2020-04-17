function [F1, F2, F3 ,Xi, lambda, ind] =  L_curve_rpmo_plotter(err,obj,lambdas,sols,plotting,solver,solver_info)

n_lambdas = length(lambdas);
n_params = length(solver_info{1}.fit_errs);
n_non_precond_params = length(solver_info{end}.fit_errs);
P = size(sols,1);

Lambdas = repmat(lambdas,n_params,1);
Lambdas = Lambdas(:);
Err = zeros(n_params*(n_lambdas-1)+n_non_precond_params,1);
Obj = zeros(n_params*(n_lambdas-1)+n_non_precond_params,1);
Sols = zeros(P,n_params*(n_lambdas-1)+n_non_precond_params);

if plotting
    close all
end

xl = [];
yl = [];

p1 = [];
p2 = [];
p3 = [];

F1 = figure;
for l = n_lambdas:-1:1
    if l == n_lambdas
        spec_inds = ((l-1)*n_params+1):(((l-1)*n_params+1)+n_non_precond_params-1);
    else    
        spec_inds = ((l-1)*n_params+1):l*n_params;
    end
    Opt_info = solver_info{l};   
    OBJ = [];
   if strcmp('SPGL1',solver)
        OBJ = Opt_info.l1_norms.*Opt_info.l1_norm;
        xlab = '$|| \mathbf{\xi}_i ||_{1}$';
   else
        OBJ = Opt_info.l0_norms.*Opt_info.l0_norm;
        xlab = '$|| \mathbf{\xi}_i ||_{0}$';
   end
   if max(isnan(OBJ)) 
       OBJ = NaN*ones(size(Opt_info.l0_norms));
       Opt_info.fit_errs = NaN*ones(size(Opt_info.fit_errs));
       Opt_info.fit_errs_norm = 1;
       fprintf('%s failed (NaN) failed on parameter = %.g \n',solver, lambdas(l))
   end
   %spec_inds
   %OBJ
   %size(spec_inds)
   %size(OBJ)
   Obj(spec_inds) = OBJ;
   
   Err(spec_inds) =  Opt_info.fit_errs.*Opt_info.fit_errs_norm;
   Sols(:,spec_inds) = Opt_info.solutions;
   
    if plotting
        hold on
        p1 = semilogy(OBJ,Opt_info.fit_errs*Opt_info.fit_errs_norm,'--o','markersize',13,'LineWidth',3);
        p1.Color(4) = 0.2;
        set(0,'DefaultLegendAutoUpdate','off')
        p2 = semilogy([min(OBJ),OBJ(Opt_info.min_ind)],[min(Opt_info.fit_errs)*Opt_info.fit_errs_norm,Opt_info.fit_errs(Opt_info.min_ind)*Opt_info.fit_errs_norm],':r+','markersize',12,'LineWidth',2);
        set(0,'DefaultLegendAutoUpdate','off')
        ylabel('$|| \dot{\mathbf{x}}_i  -\mathbf{\Theta}\mathbf{\xi}_i ||_2$','interpreter','latex')
        xlabel(xlab,'interpreter','latex')
        xl = xlim;
        yl = ylim;
        %set(gca,'FontSize',22)
        grid on
        hold on
        box on
        if l == n_lambdas
            set(p1,'Visible','off')
            set(p2,'Visible','off')
            set(gca, 'YScale', 'log','fontsize',22)
            p3 = semilogy(OBJ,Opt_info.fit_errs*Opt_info.fit_errs_norm,'-ko','LineWidth',3,'markersize',13,'MarkerFaceColor',[1,1,1]);
            p3.Color(4) = 0.2;
            p2 = semilogy([min(OBJ),OBJ(Opt_info.min_ind)],[min(Opt_info.fit_errs)*Opt_info.fit_errs_norm,Opt_info.fit_errs(Opt_info.min_ind)*Opt_info.fit_errs_norm],':r+','markersize',12,'LineWidth',2);
            %L = legend('Candidate Solutions (colors differ)','L-curve Criterion','Standard SINDy');
        end
    end
end
L = legend([p1 p2 p3],{'Cand. Sols. (colors differ)','L-curve Criterion','Standard SINDy'}) 
set(L,'interpreter','latex','fontsize',22);
set(0,'DefaultLegendAutoUpdate','off')


good_inds = find(~isnan(Obj));
Obj = Obj(good_inds); Err = Err(good_inds); Sols = Sols(:,good_inds); Lambdas = Lambdas(good_inds);
Obj_norm = norm(Obj);
Obj = Obj./Obj_norm;
Err_norm = norm(Err);
Err = Err./Err_norm;
[Obj , I ] = sort(Obj);
Err = Err(I);
Sols = Sols(:,I);
Lambdas = Lambdas(I);
[~,Min_ind] = min( sqrt((min(Err)- Err).^2 + (min(Obj) -Obj).^2));
%Xi = Sols(:,Min_ind);
Lambda = Lambdas(Min_ind);

obj_norm = norm(obj);
obj = obj./obj_norm;
err_norm = norm(err);
err = err./err_norm;

[obj , I ] = sort(obj);
err = err(I);
sols = sols(:,I);
lambdas = lambdas(I);

[~,min_ind] = min( sqrt((min(err)- err).^2 + (min(obj) -obj).^2));


lambda = lambdas(min_ind);
reg_ind = find(lambdas == 0+1i);
reg_Ind = find(Lambdas == 0+1i);

precoherenced = find(Lambdas ~= 0+1i);
Xi = sols(:,min_ind);

if plotting
    F2 = figure;
    disp(Lambda)
    %semilogy(Obj, Err,'bo',[min(Obj),Obj(Min_ind)],[min(Err),Err(Min_ind)],'-r*',Obj(reg_Ind),Err(reg_Ind),'-k+','markersize',8,'LineWidth',1);
    p1 = semilogy(Obj(precoherenced)*Obj_norm, Err(precoherenced)*Err_norm,'bo','markersize',13,'LineWidth',3);
    hold on
    p2 = semilogy(Obj(reg_Ind)*Obj_norm,Err(reg_Ind)*Err_norm,'-ko','LineWidth',3,'markersize',13,'MarkerFaceColor',[1,1,1]);
    ylabel('$|| \dot{\mathbf{x}}_i  -\mathbf{\Theta}\mathbf{\xi}_i ||_2$','interpreter','latex')
    xlabel(xlab,'interpreter','latex')
    set(gca,'FontSize',24)
    xlim(xl)
    ylim(yl) 
    L = legend('Precoherenced','Standard SINDy');
    set(L,'interpreter','latex','fontsize',24);
    %set(0,'DefaultLegendAutoUpdate','off')
    grid on
    box on
    %drawnow
    F3 = figure;
    disp(lambda)
    
    p1 = semilogy(obj*obj_norm,err*err_norm,'-b+','markersize',13,'LineWidth',3,'MarkerEdgeColor',[1, 0 ,0]);
    hold on
    p2 = semilogy(obj(reg_ind)*obj_norm,err(reg_ind)*err_norm,'ko','LineWidth',3,'markersize',13,'MarkerFaceColor',[1,1,1]);
    p3 = semilogy([min(obj),obj(min_ind)].*obj_norm,[min(err),err(min_ind)].*err_norm,':r+','markersize',13,'LineWidth',3);
    ylabel('$|| \dot{\mathbf{x}}_i  -\mathbf{\Theta}\mathbf{\xi}_i ||_2$','interpreter','latex')
    xlabel(xlab,'interpreter','latex')
    set(gca,'FontSize',24)
    xlim(xl)
    ylim(yl) 
    L = legend([p1 p3 p2],{'Candidate Solutions','L-curve Criterion','Standard SINDy'});
    set(L,'interpreter','latex','fontsize',24);
    %set(0,'DefaultLegendAutoUpdate','off')
    grid on
    box on
    %pause 
    disp(solver)
end
%ind = find(lambdas == lambda);
ind = find(lambdas == Lambda);

end