function Xi_ref = lorenz_ref_sol(A,B,P)
beta = 8/3;
sigma = 10;
rho = 28;
d = 3;

Xi_ref = zeros(P,d);
Bx = (B(1)-A(1))/2;
Ax = (B(1)+A(1))/2;
By = (B(2)-A(2))/2;
Ay = (B(2)+A(2))/2;
Bz = (B(3)-A(3))/2;
Az = (B(3)+A(3))/2;
%dwx/dt 
Xi_ref(1,1) = Ay-Ax;              % c0 
Xi_ref(2,1) = -Bx;                 % cwx 
Xi_ref(3,1) = By;                 % cwy
Xi_ref(:,1) = (sigma/Bx).*Xi_ref(:,1);

%dwydt
Xi_ref(1,2) = rho*Ax-Ax*Az - Ay;  % c0
Xi_ref(2,2) = rho*Bx-Bx*Az;              % cwx 
Xi_ref(3,2) = -By;                % cwy 
Xi_ref(4,2) = -Ax*Bz;             % cwz 
Xi_ref(8,2) = -Bx*Bz;             % cwxwz
Xi_ref(:,2) = Xi_ref(:,2)/By;
%dwzdt 
Xi_ref(1,3) = Ax*Ay-beta*Az;      % c0
Xi_ref(2,3) = Bx*Ay;              % cwx
Xi_ref(3,3) = Ax*By;              % cwy
Xi_ref(4,3) = -beta*Bz;           % cwz
Xi_ref(6,3) = Bx*By;              % cwxwy
Xi_ref(:,3) = Xi_ref(:,3)./Bz;
end