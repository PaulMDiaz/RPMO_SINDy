function [dydt] = mean_field_model(t,y)
%Function for the mean field model
% mu = 2;
% omega = 15;
% A = -0.1;
% lambda = 10;

mu = 2;
omega = 20;
A = -0.01;
lambda = 10;

if size(y,2) == 1
dydt = [mu.*y(1) - omega.*y(2) + A.*y(1).*y(3); omega.*y(1)+ mu.*y(2) + A.*y(2).*y(3);-lambda.*(y(3)-y(1)^2-y(2)^2)];
else 
dydt = [mu.*y(:,1) - omega.*y(:,2) + A.*y(:,1).*y(:,3),omega.*y(:,1)+ mu.*y(:,2) + A.*y(:,2).*y(:,3), -lambda.*(y(:,3)-y(:,1).^2-y(:,2).^2)];
end


end

