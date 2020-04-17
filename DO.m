function dxdt = DO(t, x)

xi =  0.0519;
omega = 6.63;
epsilon = -.5;


    
if size(x,2) == 1
%dydt = [sigma.*(y(2)-y(1)); rho.*y(1)-y(2)-y(1).*y(3); y(1).*y(2)-beta.*y(3)];
dxdt = [x(2); -2*omega*xi*x(2) - omega^2*(x(1) + epsilon*x(1)^3)];
else 
%dydt = [sigma.*(y(:,2)-y(:,1)), rho.*y(:,1)-y(:,2)-y(:,1).*y(:,3), y(:,1).*y(:,2)-beta.*y(:,3)];
dxdt = [x(:,2), -2.*omega.*xi.*x(:,2) - omega.^2.*(x(:,1) + epsilon*x(:,1).^3)];
end


end