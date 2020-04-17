function Xi_ref = MFM_ref_sol(A,B,P)

% mu = 2;
% omega = 15;
% a = -0.1;
% lambda = 10;
mu = 2;
omega = 20;
a = -0.01;
lambda = 10;
d = 3;

Xi_ref = zeros(P,d);
Bx = (B(1)-A(1))/2;
Ax = (B(1)+A(1))/2;
By = (B(2)-A(2))/2;
Ay = (B(2)+A(2))/2;
Bz = (B(3)-A(3))/2;
Az = (B(3)+A(3))/2;

%dwx/dt 
Xi_ref(1,1) = mu*Ax-omega*Ay+a*Ax*Az;           % c0 
Xi_ref(2,1) = mu*Bx+a*Bx*Az;                    % cwx 
Xi_ref(3,1) = -omega*By;                        % cwy
Xi_ref(4,1) = a*Ax*Bz;                          % cwz
Xi_ref(8,1) = a*Bx*Bz;                          % cwxwz
Xi_ref(:,1) = Xi_ref(:,1)./Bx;

%dwy/dt
Xi_ref(1,2) = omega*Ax+mu*Ay+a*Ay*Az;           % c0
Xi_ref(2,2) = omega*Bx;                         % cwx 
Xi_ref(3,2) = mu*By+a*Az*By;                    % cwy 
Xi_ref(4,2) = a*Ay*Bz;                          % cwz 
Xi_ref(9,2) = a*By*Bz;                          % cwywz
Xi_ref(:,2) = Xi_ref(:,2)./By;

%dwz/dt
Xi_ref(1,3) = Az-Ax^2-Ay^2;                    % c0
Xi_ref(2,3) = -2*Bx*Ax;                           % cwx
Xi_ref(3,3) = -2*By*Ay;                           % cwy
Xi_ref(4,3) = Bz;                              % cwz
Xi_ref(5,3) = -Bx^2;                           % cwx^2
Xi_ref(7,3) = -By^2;                           % cwy^2
Xi_ref(:,3) = Xi_ref(:,3).*(-lambda/Bz);





end