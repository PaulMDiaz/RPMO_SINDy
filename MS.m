function [dydt] = MS(t,y)
%Function for the lorenz 63 system
R = 70;
G = 40;

if size(y,2) == 1
%dydt = [sigma.*(y(2)-y(1)); rho.*y(1)-y(2)-y(1).*y(3); y(1).*y(2)-beta.*y(3)];
dydt = [y(2); -y(2)+R.*y(1)-G.*(y(1)+y(3))-R.*y(1).*y(3)^2; y(1)];
else 
%dydt = [sigma.*(y(:,2)-y(:,1)), rho.*y(:,1)-y(:,2)-y(:,1).*y(:,3), y(:,1).*y(:,2)-beta.*y(:,3)];
dydt = [y(:,2), -y(:,2)+R.*y(:,1)-G.*(y(:,1)+y(:,3))-R.*y(:,1).*y(:,3).^2, y(:,1)];
end


end

