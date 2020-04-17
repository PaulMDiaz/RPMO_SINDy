function Xi_ref = MS_ref_sol(A,B,P)
R = 70;
G = 40;

d = 3;
Xi_ref = zeros(P,d);
Bx = (B(1)-A(1))/2;
Ax = (B(1)+A(1))/2;
By = (B(2)-A(2))/2;
Ay = (B(2)+A(2))/2;
Bz = (B(3)-A(3))/2;
Az = (B(3)+A(3))/2;
%dwx/dt 
Xi_ref(1,1) = Ay/Bx;  % c0 
Xi_ref(3,1) = By/Bx;  % cwx 

%dwydt
Xi_ref(1,2) = -Ay/By + R*Ax/By -G*Ax/By - G*Az/By - R*Ax*Az^2/By;  % c0
Xi_ref(2,2) = R*Bx/By - G*Bx/By -R*Bx*Az^2/By;                     % cwx 
Xi_ref(3,2) = -1;                                                   % cwy 
Xi_ref(10,2) = -R*Ax*Bz^2/By;                                      % cwz^2
Xi_ref(8,2) = -2*R*Bx*Bz*Az/By;                                    % cwxwz
Xi_ref(4,2) = -G*Bz/By - 2*R*Ax*Bz*Az/By;                          % cwz 
Xi_ref(18,2) = -R*Bx*Bz^2/By;                                      % cwxwz^2                       

%dwzdt 
Xi_ref(1,3) = Ax/Bz; % c0
Xi_ref(2,3) = Bx/Bz; % cwx


end