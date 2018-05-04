function u = rk6(f,u,t,dt)

        % Runge-Kutta 6 from Alshina07 
        s = 7;
        b = [1/12 0 0 0 5/12 5/12 1/12];
        a = zeros(7,6);
        a(2,1) = 4/7; 
        a(3,1) = 115/112; a(3,2) = -5/16;
        a(4,1) = 589/630; a(4,2) = 5/18; a(4,3) = -16/45;
        a(5,1) = 229/1200 - 29/6000*sqrt(5); a(5,2) = 119/240 - 187/1200*sqrt(5); a(5,3) = -14/75 + 34/375*sqrt(5); a(5,4) = -3/100*sqrt(5);
        a(6,1) = 71/2400 - 587/12000*sqrt(5); a(6,2) = 187/480 - 391/2400*sqrt(5); a(6,3) = -38/75 + 26/375*sqrt(5); a(6,4) = 27/80 - 3/400*sqrt(5); a(6,5) = (1+sqrt(5))/4;
        a(7,1) = -49/480 + 43/160*sqrt(5); a(7,2) = -425/96 + 51/32*sqrt(5); a(7,3) = 52/15 - 4/5*sqrt(5); a(7,4) = -27/16 + 3/16*sqrt(5); a(7,5) = 5/4 - 3/4*sqrt(5); a(7,6) = 5/2 - 1/2*sqrt(5);

        c = [0, 4/7, 5/7, 6/7, (5-sqrt(5))/10, (5+sqrt(5))/10, 1];

N = size(u,1);
K = zeros(N,s); 
% Matris för att lagra stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%% Time-stepping %%%%%%%%%%%%%
for mm = 1:s
    u_temp = u;
    for l = 1:mm-1
        u_temp = u_temp + dt*a(mm,l)*K(:,l);
    end
    K(:,mm) = f(u_temp,t+dt*c(mm));
end

% Update
for mm=1:s
    u = u + dt*b(mm)*K(:,mm);
end
        