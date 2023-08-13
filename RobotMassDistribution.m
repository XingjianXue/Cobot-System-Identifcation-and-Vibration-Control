% Numerical computation of the kinematics and dynamics of the UR5e
% Nosa Edoimioya
% 09/12/2022

close all
clear all
clc

% Declare symbols
syms g q_1 q_2 q_3 q_4 q_5 q_6...
    dq_1 dq_2 dq_3 dq_4 dq_5 dq_6
%% DH and Inertia Parameters of UR5e - make sure all joint angles start from 0 rad
alpha = [0,pi/2,0,0,pi/2,-pi/2,0]; %[rad]
a = [0,0,-0.425,-0.3922,0,0,0]; %[m]
d = [0.1625,0,0,0.1333,0.0997,0.0996,0]; %[m]
m = [3.761,8.058,2.846,1.37,1.3,0.365,0.720]; %[kg]

p_c1 = [0,-25.61,1.93]*10^(-3); %[m]
p_c2 = [-212.5,0,113.36]*10^(-3);
p_c3 = [-150,0,26.5]*10^(-3);
p_c4 = [0,-1.8,16.34]*10^(-3);
p_c5 = [0,1.8,16.34]*10^(-3);
p_c6 = [0,0,-1.159]*10^(-3);

p_c1 = [0,1.93,-25.61]*10^(-3); %[m]
p_c2 = [0,-24.201,212.5]*10^(-3);
p_c3 = [0,26.5,119.93]*10^(-3);
p_c4 = [0,110.9,16.34]*10^(-3);
p_c5 = [0,1.8,110.9]*10^(-3);
p_c6 = [0,1.159,0]*10^(-3);

p_c1 = [0,-25.61,1.93]*10^(-3); %[m]
p_c2 = [-212.5,0,113.36]*10^(-3);
p_c3 = [-119,0,26.5]*10^(-3);
p_c4 = [0,16.34,-1.8]*10^(-3);
p_c5 = [0,1.8,16.34]*10^(-3);
p_c6 = [0,0,-1.159]*10^(-3);
p_c7 = [10,-32,27]*10^(-3); % nozzle
In_1 = [84,0,0;0,64,0;0,0,84]*10^(-4); % these are not measured
In_2 = [78,0,0;0,21,0;0,0,21]*10^(-4);
In_3 = [16,0,0;0,462,0;0,0,462]*10^(-4);
In_4 = [16,0,0;0,16,0;0,0,9]*10^(-4);
In_5 = [16,0,0;0,16,0;0,0,9]*10^(-4);
In_6 = eye(3)*10^(-4);
In_7 = [8.90,0,0;0,8.9,0;0,0,8.9]*10^(-4);
%% OTHER PARAMETERs AND SYMBOLs
alpha_0=alpha(1);alpha_1=alpha(2);alpha_2=alpha(3);alpha_3=alpha(4);alpha_4=alpha(5);alpha_5=alpha(6);alpha_6=alpha(7);
a_0=a(1);a_1=a(2);a_2=a(3);a_3=a(4);a_4=a(5);a_5=a(6);a_6=a(7);
d_1=d(1);d_2=d(2);d_3=d(3);d_4=d(4);d_5=d(5);d_6=d(6);d_7=d(7);
p_cx1=p_c1(1);p_cy1=p_c1(2);p_cz1=p_c1(3);
p_cx2=p_c2(1);p_cy2=p_c2(2);p_cz2=p_c2(3);
p_cx3=p_c3(1);p_cy3=p_c3(2);p_cz3=p_c3(3);
p_cx4=p_c4(1);p_cy4=p_c4(2);p_cz4=p_c4(3);
p_cx5=p_c5(1);p_cy5=p_c5(2);p_cz5=p_c5(3);
p_cx6=p_c6(1);p_cy6=p_c6(2);p_cz6=p_c6(3);
p_cx7=p_c7(1);p_cy7=p_c7(2);p_cz7=p_c7(3);
m_1=m(1);m_2=m(2);m_3=m(3);m_4=m(4);m_5=m(5);m_6=m(6);m_7=m(7);

%% ROTATION MATRICEs
R_1=[cos(q_1) -sin(q_1) 0;
     sin(q_1)*cos(alpha_0) cos(q_1)*cos(alpha_0) -sin(alpha_0);
     sin(q_1)*sin(alpha_0) cos(q_1)*sin(alpha_0)  cos(alpha_0)];
R_2=[cos(q_2) -sin(q_2) 0;
     sin(q_2)*cos(alpha_1) cos(q_2)*cos(alpha_1) -sin(alpha_1);
     sin(q_2)*sin(alpha_1) cos(q_2)*sin(alpha_1)  cos(alpha_1)];
R_3=[cos(q_3) -sin(q_3) 0;
     sin(q_3)*cos(alpha_2) cos(q_3)*cos(alpha_2) -sin(alpha_2);
     sin(q_3)*sin(alpha_2) cos(q_3)*sin(alpha_2)  cos(alpha_2)];
R_4=[cos(q_4) -sin(q_4) 0;
     sin(q_4)*cos(alpha_3) cos(q_4)*cos(alpha_3) -sin(alpha_3);
     sin(q_4)*sin(alpha_3) cos(q_4)*sin(alpha_3)  cos(alpha_3)];
R_5=[cos(q_5) -sin(q_5) 0;
     sin(q_5)*cos(alpha_4) cos(q_5)*cos(alpha_4) -sin(alpha_4);
     sin(q_5)*sin(alpha_4) cos(q_5)*sin(alpha_4)  cos(alpha_4)];
R_6=[cos(q_6) -sin(q_6) 0;
     sin(q_6)*cos(alpha_5) cos(q_6)*cos(alpha_5) -sin(alpha_5);
     sin(q_6)*sin(alpha_5) cos(q_6)*sin(alpha_5)  cos(alpha_5)];
R_7=[cos(0) -sin(0) 0;
     sin(0)*cos(alpha_6) cos(0)*cos(alpha_6) -sin(alpha_6);
     sin(0)*sin(alpha_6) cos(0)*sin(alpha_6)  cos(alpha_6)];

%% POSITION VECTORs
p_1=[a_0;-sin(alpha_0)*d_1;cos(alpha_0)*d_1];
p_2=[a_1;-sin(alpha_1)*d_2;cos(alpha_1)*d_2];
p_3=[a_2;-sin(alpha_2)*d_3;cos(alpha_2)*d_3];
p_4=[a_3;-sin(alpha_3)*d_4;cos(alpha_3)*d_4];
p_5=[a_4;-sin(alpha_4)*d_5;cos(alpha_4)*d_5];
p_6=[a_5;-sin(alpha_5)*d_6;cos(alpha_5)*d_6];
p_7=[a_6;-sin(alpha_6)*d_7;cos(alpha_6)*d_7];
%% TRANSLATION MATRICES AND FORWARD KINEMATICS
T_1 = [R_1,p_1;zeros(1,3),1];
T_2 = [R_2,p_2;zeros(1,3),1];
T_3 = [R_3,p_3;zeros(1,3),1];
T_4 = [R_4,p_4;zeros(1,3),1];
T_5 = [R_5,p_5;zeros(1,3),1];
T_6 = [R_6,p_6;zeros(1,3),1];
T_7 = [R_7,p_7;zeros(1,3),1];
T = T_1*T_2*T_3*T_4*T_5*T_6*T_7;
T_61 = T_2*T_3*T_4*T_5*T_6;
%% COMs' POSITION VECTORs
p_c1=p_1+R_1*[p_cx1;p_cy1;p_cz1];
p_c2=p_1+R_1*(p_2+R_2*[p_cx2;p_cy2;p_cz2]);
p_c3=p_1+R_1*(p_2+R_2*(p_3+R_3*[p_cx3;p_cy3;p_cz3]));
p_c4=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*[p_cx4;p_cy4;p_cz4])));
p_c5=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*(p_5+R_5*([p_cx5;p_cy5;p_cz5])))));
p_c6=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*(p_5+R_5*(p_6+R_6*[p_cx6;p_cy6;p_cz6])))));
p_c7=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*(p_5+R_5*(p_6+R_6*(p_7+R_7*[p_cx7;p_cy7;p_cz7]))))));

%% SYSTEM's STATEs
q=[q_1;q_2;q_3;q_4;q_5;q_6];
% dq=[dq_1;dq_2;dq_3;dq_4;dq_5;dq_6];

%% LINEAR PART of JACOBIANs
J_v1=jacobian(p_c1,q);
J_v2=jacobian(p_c2,q);
J_v3=jacobian(p_c3,q);
J_v4=jacobian(p_c4,q);
J_v5=jacobian(p_c5,q);
J_v6=jacobian(p_c6,q);
J_v7=jacobian(p_c7,q);
%% ROTATION MATRICEs FROM BASE
R_20=R_1*R_2;
R_30=R_20*R_3;
R_40=R_30*R_4;
R_50=R_40*R_5;
R_60=R_50*R_6;
R_70=R_60*R_7;
%% ANGULAR PART of JACOBIANs
%o=zeros(3,7);
J_o1=[R_1(:,3),zeros(3,5)];
J_o2=[R_1(:,3),R_20(:,3),zeros(3,4)];
J_o3=[R_1(:,3),R_20(:,3),R_30(:,3),zeros(3,3)];
J_o4=[R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),zeros(3,2)];
J_o5=[R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),R_50(:,3),zeros(3,1)];
J_o6=[R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),R_50(:,3),R_60(:,3)];
%% JACOBIAN MATRIX OF THE END-EFFECTOR
Jacobi = [J_v7;J_o6];
%% ROBOT's INERTIA (MASS) MATRIX
M=J_v1.'*m_1*eye(3)*J_v1+J_o1.'*R_1*In_1*R_1.'*J_o1...
 +J_v2.'*m_2*eye(3)*J_v2+J_o2.'*R_20*In_2*R_20.'*J_o2...
 +J_v3.'*m_3*eye(3)*J_v3+J_o3.'*R_30*In_3*R_30.'*J_o3...
 +J_v4.'*m_4*eye(3)*J_v4+J_o4.'*R_40*In_4*R_40.'*J_o4...
 +J_v5.'*m_5*eye(3)*J_v5+J_o5.'*R_50*In_5*R_50.'*J_o5...
 +J_v6.'*m_6*eye(3)*J_v6+J_o6.'*R_60*In_6*R_60.'*J_o6...
 +J_v7.'*m_7*eye(3)*J_v7+J_o6.'*R_70*In_7*R_70.'*J_o6;
%%
q_1 = 0;
q_2 = -pi;
q_3 = 0;
q_4 = -pi/2;
q_5 = 0;
q_6 = 0;

double(subs(M))

%% SAVE AS MATLAB FUNCTION

% matlabFunction(M,'File','ur5e_inertia_matrix_nozzle');
% matlabFunction(Jacobi,'File','ur5e_jacobian_nozzle');
% matlabFunction(T_61,'File','ur5e_T61_transformation_matrix');
%% CORIOLIS and CENTRIFUGAL MATRIX
for k=1:6
   for s=1:6
      C(k,s)=.5*((diff(M(k,s),q_1)+diff(M(k,1),q(s,1))-diff(M(1,s),q(k,1)))*dq_1...
                +(diff(M(k,s),q_2)+diff(M(k,2),q(s,1))-diff(M(2,s),q(k,1)))*dq_2...
                +(diff(M(k,s),q_3)+diff(M(k,3),q(s,1))-diff(M(3,s),q(k,1)))*dq_3...
                +(diff(M(k,s),q_4)+diff(M(k,4),q(s,1))-diff(M(4,s),q(k,1)))*dq_4...
                +(diff(M(k,s),q_5)+diff(M(k,5),q(s,1))-diff(M(5,s),q(k,1)))*dq_5...
                +(diff(M(k,s),q_6)+diff(M(k,6),q(s,1))-diff(M(6,s),q(k,1)))*dq_6);
   end
end
%% POTENTIAL ENERGIES and GRAVITY VECTOR
P1=m_1*[0,0,g]*p_c1;
P2=m_2*[0,0,g]*p_c2;
P3=m_3*[0,0,g]*p_c3;
P4=m_4*[0,0,g]*p_c4;
P5=m_5*[0,0,g]*p_c5;
P6=m_6*[0,0,g]*p_c6;
P7=m_7*[0,0,g]*p_c7;
P=P1+P2+P3+P4+P5+P6+P7;
g_1=diff(P,q_1);
g_2=diff(P,q_2);
g_3=diff(P,q_3);
g_4=diff(P,q_4);
g_5=diff(P,q_5);
g_6=diff(P,q_6);
G=[g_1;g_2;g_3;g_4;g_5;g_6];

%% DYNAMICAL EQUATIONs of MOTION
% % % % %           M(q)*ddq + C(q,dq)dq + G(q) = u
save ('UR5e_num.mat','T','Jacobi','M','C','G');

fid = fopen('UR5eT_num.txt', 'w');
fwrite(fid, char(T), 'char');
fclose(fid);

fid = fopen('UR5eM_num.txt', 'w');
fwrite(fid, char(M), 'char');
fclose(fid);

fid = fopen('UR5eC_num.txt', 'w');
fwrite(fid, char(C), 'char');
fclose(fid);

fid = fopen('UR5eG_num.txt', 'w');
fwrite(fid, char(G), 'char');
fclose(fid);

fid = fopen('UR5eJ_num.txt', 'w');
fwrite(fid, char(Jacobi), 'char');
fclose(fid);

clear 