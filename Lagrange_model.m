%Euler Lagrange Model
clc
clear

syms M11 M12 M13 M22 M23 M33 K1 K2 K3 R1 R2 R3 Rg C1 C2 C3 s











%(K2*K3 - M23^2*s^4 + K3*M22*s^2 + K2*M33*s^2 + M22*M33*s^4 + C2*M22*M33*s^3 + C3*M22*M33*s^3 + C2*K3*M22*s + C3*K2*M33*s + C2*C3*M22*M33*s^2)



% M = [ M22/K2 M23/K2;
%       M23/K3 M33/K3];
% 
% C = [ C2*M22/K2 0;
%       0 C3*M33/K3];
% K = [ 1 0;
%       0 1];
% G = inv(M*s^2 +C*s + K)*[1;0]

% mass = ur5e_inertia_matrix_5kg(0,0,-(k-1)*pi/9,0,0,0);
% M11 = mass(1,1);
% M12 = mass(1,2);
% M22 = mass(2,2);
% M23 = mass(2,3);
% M33 = mass(3,3);
% M13 = mass(1,3);
% 
% K1 = 67*7.65;
% K2 = 150*M22;
% K3 = 300*M33;
% C1 = 6;
% C2 = 6;
% C3 = 6;
% 
% 
% p = [M22*M33-M23^2
%     C2*M22*M33+C3*M22*M33
%     K3*M22+K2*M33+C2*C3*M22*M33
%     C2*K3*M22+C3*K2*M33 
%     K2*K3] ;
% roots(p)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% base prediction
for i = 0:8
M = [ M11/K1 M12/K1 M13/K1;
      M12/K2 M22/K2 M23/K2;
      M13/K3 M23/K3 M33/K3];

C = [ C1*M11/K1 0 0;
      0 C2*M22/K2 0;
      0 0 C3*M33/K3];
K = [ 1 0 0;
      0 1 0;
      0 0 1];
G = inv(M*s^2 +C*s + K)*[0;1;0]

mass = ur5e_inertia_matrix_5kg(0,0*pi/9,-i*pi/9,0,0,0)
M11 = mass(1,1);
M12 = mass(1,2);
M22 = mass(2,2);
M23 = mass(2,3);
M33 = mass(3,3);
M13 = mass(1,3);

K1 = 67*7.65;
K2 = 550;
K3 = 1200;
C1 = 6;
C2 = 8;
C3 = 10;

p=[-M11*M23^2-M13^2*M22-M12^2*M33+2*M12*M13*M23+M11*M22*M33  
    -C1*M11*M23^2-C2*M13^2*M22-C3*M12^2*M33+C1*M11*M22*M33+C2*M11*M22*M33+C3*M11*M22*M33     
    -K3*M12^2-K1*M23^2-K2*M13^2+K3*M11*M22+K2*M11*M33+K1*M22*M33+C1*C2*M11*M22*M33+C1*C3*M11*M22*M33+C2*C3*M11*M22*M33    
    +C1*K3*M11*M22+C2*K3*M11*M22+C1*K2*M11*M33+C3*K2*M11*M33+C2*K1*M22*M33+C3*K1*M22*M33+C1*C2*C3*M11*M22*M33 
    +K2*K3*M11+K1*K3*M22+K1*K2*M33+C1*C2*K3*M11*M22+C1*C3*K2*M11*M33+C2*C3*K1*M22*M33  
    +C1*K2*K3*M11+C2*K1*K3*M22+C3*K1*K2*M33       
    K1*K2*K3]

complex_root = roots(p)
base_root = complex_root(5)
c = -2*real(base_root)
d = real(base_root)^2 + imag(base_root)^2


% q = [M22*M33-M23^2  C2*M22*M33+C3*M22*M33  K3*M22+K2*M33+C2*C3*M22*M33     C2*K3*M22+C3*K2*M33  K2*K3]
% roots(q)


sqrt(K1/M11)
end










