function Xi_ref = DO_ref_sol(A,B,P)
xi =  0.0519;
omega = 6.63;
epsilon = -.5;

d = 2;


Xi_ref = zeros(P,d);
Bx = (B(1)-A(1))/2;
Ax = (B(1)+A(1))/2;
By = (B(2)-A(2))/2;
Ay = (B(2)+A(2))/2;

%dwx/dt 
Xi_ref(1,1) = Ay/Bx;  % c0 
Xi_ref(3,1) = By/Bx;  % cwx 

%dwydt
Xi_ref(1,2) = -2*omega*xi*Ay - omega^2*(Ax+epsilon*Ax^3);                    % c0
Xi_ref(2,2) = -omega^2*(Bx+epsilon*3*Ax^2*Bx);                            % cwx 
Xi_ref(3,2) = -2*omega*xi*By;                                             % cwy
Xi_ref(4,2) =  -omega^2*epsilon*3*Ax*Bx^2;                                % cx^2
Xi_ref(7,2) = -omega^2*Bx^3*epsilon;                                     % cx^3
Xi_ref(:,2) = Xi_ref(:,2)./By;





end