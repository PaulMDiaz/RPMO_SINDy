function [T,X,V_ref,D_M,D_L,scaling,A,B] = solve_ode(X0,T,func,index_pc,Xi)
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    if nargin == 5
        [T, X] = ode45(@(T,X) func(T,X,Xi,index_pc), T, X0,options);
        V_exact = func(mean(T), X, Xi, index_pc);
    else
        [T, X] = ode45(@(T,X) func(T,X), T, X0,options);
        V_exact = func(mean(T), X);
    end
    

    A = min(X);
    B = max(X);
    [X,scaling] = shift_and_scale_data(X); 
 
    [D_M,D_L] = build_dictionary(X,index_pc);
    %Xi_ref_M  = MS_ref_sol(A,B,size(D_L,2));
    %V_ref = [V_exact(:,1)*scaling(1),V_exact(:,2)*scaling(2),V_exact(:,3)*scaling(3)];
   
    V_ref = zeros(size(V_exact));
    for i = 1:size(X,2);
        V_ref(:,i) = V_exact(:,i)*scaling(i);
    end
        
%     approx_truncation_error = norm(D_M*Xi_ref_M-V_ref)/norm(V_ref);
%     if approx_truncation_error > 1e-10
%         fprintf('Warning, Ode45 relative validation error is %f of ref soln  results may be inaccurate. \n',approx_truncation_error)
%     end
end