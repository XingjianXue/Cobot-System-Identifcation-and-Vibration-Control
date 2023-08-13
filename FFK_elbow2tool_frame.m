function r = FFK_elbow2tool_frame(q_ref)
%% Develop by Iago - 06/12/2023
% This script will be useful to calculate foward kinematics of the robot.
% r is the forward kinematics from the elbow to the tool-frame.
% The script was writen/validated using two references
% Equation (n-1)T_n - First transformation matrix from:
% https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters 
% Parameters and case to validate:
% https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/

%% D-H parameters - Reference: https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
%       a[m]    d[m]    alpha[rad] -        a = Ur5e(joint_number,1)
Ur5e = [0    0.1625     pi/2;    % Joint 1  d = Ur5e(joint_number,2)
    -0.425     0         0;      % Joint 2  alpha = Ur5e(joint_number,3)   
    -0.3922    0         0;      % Joint 3
        0    0.1333     pi/2;    % Joint 4
        0    0.0997    -pi/2;    % Joint 5
        0    0.0996      0];     % Joint 6

%q_ref  = [0,0,0,0,0,0]; % reference configuration
%q_ref  = [302,190,59,199,120,90]*pi/180; % validation configuration
T_Ur5e = zeros(4,4,6);

for i = 1:6
    
    a             = Ur5e(i,1);
    d             = Ur5e(i,2);
    alpha         = Ur5e(i,3);
    T_Ur5e(:,:,i) = [ cos(q_ref(i)) -sin(q_ref(i))*cos(alpha) sin(q_ref(i))*sin(alpha)  a*cos(q_ref(i));
                      sin(q_ref(i)) cos(q_ref(i))*cos(alpha)  -cos(q_ref(i))*sin(alpha) a*sin(q_ref(i));
                            0                sin(alpha)                cos(alpha)              d       ;
                            0                    0                         0                   1      ];  

end

r = T_Ur5e(:,:,4)*T_Ur5e(:,:,5)*T_Ur5e(:,:,6); % Rigid body transformation matrix from the elbow joint to the tool frame

