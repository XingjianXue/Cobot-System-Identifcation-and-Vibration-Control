% Read system identification data from automatic system ID
% Nosa Edoimioya
% 10/14/2022

% Edited by Iago - 12-Jun-2023
% Nosa has done the code to extract the FRFs from the shoulder and base
% joints. Now, I will add the elbow joint.

% Edited by Iago - 28-Jun-2023
% In this new version, I am calculating and saving 3 types of FRFs
% End-effector position/Rotory encoder - [rad]/[rad] - H
% End-effector position/desired joint trajectory - [rad]/[rad] - H_full
% Rotory encoder/desired joint trajectory - [rad]/[rad] - H_joint
 
% Edited by Xingjian - July-August 2023
% Modelling, fitting and prediction

close all
clear all
clc

global Ts

%% Parameters
fs = 500; % sampling rate [Hz]
Ts = 1/fs; % computation step size [s]
p  = 10; % number of sine cycles

freq_range          = 2:0.5:60; % [Hz]
max_position_stroke = 0.045; % [rad]
var_acc_max         = 10; % [rad/s] - Data from automatic system ID scrip with joint controller - python 

%% Reading files
dir_date_str = 'base_final'; % Directory where the experimental data is
date_str     = 'save'; % Directory where the FRFs will be saved 
data_dir_string = ['C:\Users\xjys\Desktop\Code\',dir_date_str] % Folder where the data is placed
frf_dir_string  = ['C:\Users\xjys\Desktop\Code\',date_str] % Folder to sabe the FRFs
% lines = readlines(fullfile(data_dir_string,'*.csv'));
files           = dir(fullfile(data_dir_string,'*.csv'));
filenames       = {files.name};

%% Filter data

fc_lp = 80; % Low pass filter cutoff frequency [Hz]
[B,A] = butter(4,fc_lp*2/fs);

%% Loop through files and process data

joint_positions = zeros(length(filenames),6); % for storing joint positions

for k = 1:length(filenames)
    filename = char(fullfile(data_dir_string,filenames{k}));
    % remove .csv
    end_filename_idx = strfind(filename,'.csv');
    filename = filename(1:end_filename_idx-1);

    % find the joint positions by parsing the string
    pos_idx = strfind(filename,'pos'); % returns a vector of the starting index
    end_idx = strfind(filename,'joint_');
    sub_str = filename(pos_idx+4:end_idx(2)-1);
    sub_str = replace(sub_str,'pt','.');

    V_idx     = strfind(filename,'_V'); % returns a vector of the starting index 
    V_idx_end = strfind(filename,'rad');
    R_idx     = strfind(filename,'_R');
    R_idx_end = strfind(filename,'m_');
    
    V_rad     = filename(V_idx+2:V_idx_end-1);
    V_rad     = str2double(replace(V_rad,'pt','.'));
    
    R_m     = filename(R_idx+2:R_idx_end-1);
    R_m     = str2double(replace(R_m,'pt','.'));

    spaces = strfind(sub_str,'_');
    q_ref = zeros(1,7); % 6 joints on UR5e but there is 1 repeted joint in the file name
    for j = 1:7
        if j == 1
            q_ref(j) = str2double(sub_str(1:spaces(j)-1));
        elseif j==7
            q_ref(j) = str2double(sub_str(spaces(6)+1:end));
        else
            q_ref(j) = str2double(sub_str(spaces(j-1)+1:spaces(j)-1));
        end
    end
    
    q_ref(5) = [];
    joint_positions(k,:) = q_ref;
    % get the commanded joint (zero-indexed so add 1)
    commaded_idx = str2double(filename(end))+1;
    
    % rotation, transformation, and jacobian matrices
    R = ur5e_R60_rotation_matrix(q_ref(1),q_ref(2),q_ref(3),q_ref(4),q_ref(5),q_ref(6)); % w.r.t the universal frame
    if commaded_idx == 1
        T         = ur5e_T60_transformation_matrix(q_ref(1),q_ref(2),q_ref(3),q_ref(4),q_ref(5),q_ref(6));
        joint_str = 'base';
        r             = norm(T(1:3,4)); % distance of tool position from commanded joint axis [m]
    elseif commaded_idx == 2
        T         = ur5e_T61_transformation_matrix(q_ref(2),q_ref(3),q_ref(4),q_ref(5),q_ref(6));
        joint_str = 'shoulder';
        r             = norm(T(1:3,4)); % distance of tool position from commanded joint axis [m]
    elseif commaded_idx == 3
        T         = FFK_elbow2tool_frame(q_ref);
        joint_str = 'elbow';
        r          = norm(T(1:3,4)); % distance of tool position from commanded joint axis [m]
    end

    J = ur5e_jacobian(q_ref(1),q_ref(2),q_ref(3),q_ref(4),q_ref(5));
    
    % read data and store in arrays
    M = readmatrix(filename);
    % 1 = Time [s]
    % 2 = Commanded Base Joint Position [rad] - SineSweep after controlbox interpolation
    % 3 = Motor Output Base Position [rad]
    % 4 = Tool acceleration (x-axis) [m/s^2]
    % 5 = Tool acceleration (y-axis) [m/s^2]
    % 6 = Tool acceleration (z-axis) [m/s^2]
    % 7 = Desired Base Joint Position [rad] - SineSweep from the Python code

    t_des         = M(:,1); % [s]
    joint_des     = M(:,7); %[rad]
    joint_out     = M(:,3); %[rad]
    joint_command = M(:,2); %[rad]
    tool_acc_x    = filtfilt(B,A,M(:,4)); % [m/s^2]
    tool_acc_y    = filtfilt(B,A,M(:,5)); % [m/s^2]
    tool_acc_z    = filtfilt(B,A,M(:,6)); % [m/s^2]

    %add
    % joint1_out    = filtfilt(B,A,M(:,9));
    % joint2_out    = filtfilt(B,A,M(:,11));

    joint1_out    = M(:,9);
    joint2_out    = M(:,11);

    joint1_des    = M(:,8);
    joint2_des    = M(:,10);


    clear M
    
    % Pre-process data
    g = -9.81; % gravitational acceleration in universal frame
    G = [0; 0; g];
    
    partial_joint = J(:,commaded_idx); % respective column of Jacobian matrix, i.e. partial derivatives of tool position wrt to commanded joint
    joint_acc     = zeros(size(t_des));
    
    for i = 1:length(t_des)
        tool_acc = [tool_acc_x(i); tool_acc_y(i); tool_acc_z(i)]; % acceleration in tool frame [m/s^2]
    
        uni_acc  = R*tool_acc - G; % acceleration in universal frame [m/s^2]

        % tangential acceleration of tool position relative to the commanded joint [m/s^2]
        % assumes 0 acceleration in the orientation directions, i.e., in
        % ddot{[x,y,z,phi,theta,psi]}, ddot{[phi,theta,psi]} = 0
        tang_acc     = ([uni_acc; 0; 0; 0].'*partial_joint)/norm(partial_joint); 
        joint_acc(i) = tang_acc/r; % [rad/s^2]
    end
    
    % process data with FFT
    INPUT      = joint_des - joint_des(1); % [rad]
    INPUT_temp = INPUT;
    INPUT2     = joint_out; 
    
    OUTPUT     = joint_acc; % [rad/s^2]

    %OUTPUT     = joint1_out; %base


    

    f_start_init = 0; % f_start initial guess
    for freq_i = 1:length(freq_range)
        sin_frequency = freq_range(freq_i);
        % find pure sine index
        temp = [];
        for j = 1:ceil((1/sin_frequency*p)/Ts)
            temp(j) = (j-1)*Ts;
        end

        max_accel = ceil(2*pi*sin_frequency*sqrt(max_position_stroke));
        if max_accel < var_acc_max
            sin_amplitude = max_accel;
        else
            sin_amplitude = var_acc_max;
        end

        q_temp = -(sin_amplitude./(2*pi*sin_frequency)^2).*sin(2*pi*sin_frequency*temp); % displacement [rad]

        [c1,lags1] = xcorr(INPUT_temp, q_temp); % correlation to find the points in the desired trajectory
        sin_index1 = find(c1 == max(c1));
        sin_index1 = lags1(sin_index1)+1; % start of sine input
        f_start    = sin_index1(1); % start index number for this frequency
        f_end      = sin_index1(1) + length(q_temp) - 1; % end index number for this frequency
        f_start_init = f_start;
        INPUT_temp(1:f_start_init) = 0; % zero out the used values
        index_sin_fft = [f_start:f_end];

        % DFT
        if length(INPUT) < f_end
            input  = INPUT(f_start:end);
            input2 = INPUT2(f_start:end);
            output = OUTPUT(f_start:end);
            index_sin_fft = [f_start:length(INPUT)];
        else
            input  = INPUT(f_start:f_end);
            input2 = INPUT2(f_start:f_end);
            output = OUTPUT(f_start:f_end);
        end

        L      = length(input);
        NFFT   = L;

        in_fft = fft(input,NFFT)/length(input); % input
        in_fft = in_fft(1:NFFT/2+1);

        in_fft2 = fft(input2,NFFT)/length(input2); % input2
        in_fft2 = in_fft2(1:NFFT/2+1);

        out_fft = fft(output,NFFT)/length(output); % output
        out_fft = out_fft(1:NFFT/2+1);

        f_fft = fs/2*linspace(0,1,NFFT/2+1);

        % length(in_fft2)
        % length(out_fft)
        % length(f_fft)


        f_fft = f_fft';
        fPt   = find(abs(freq_range(freq_i) - f_fft) == min(abs(freq_range(freq_i) - f_fft)));

        w   = sin_frequency*2*pi; % [rad/s]
        s   = sqrt(-1)*w;
        Gdi = 1/s;
        
        % Note that we need to integrate the FRF two time to ensure we will
        % have and FRF that gives the output as position rather than
        % acceleration.

        % End-effector position/Rotory encoder - [rad]/[rad] - H
        % End-effector position/desired joint trajectory - [rad]/[rad] - H_full
        % Rotory encoder/desired joint trajectory - [rad]/[rad] - H_joint

        if commaded_idx == 1 % Exciting the system using the base joint
            H_Y1(freq_i)       = out_fft(fPt)/in_fft2(fPt); % output/input [rad/s^2]/[rad]
           
            H_Y1(freq_i)       = H_Y1(freq_i)*Gdi*Gdi;      % integrate twice -> output/input [rad]/[rad]
            H_Y1_full(freq_i)  = out_fft(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad]
            H_Y1_full(freq_i)  = H_Y1_full(freq_i)*Gdi*Gdi; % integrate twice -> output/input [rad]/[rad]
            H_Y1_joint(freq_i) = in_fft2(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad] end-effector
           

        elseif commaded_idx == 2 % Exciting the system using the shoulder joint
            H_Y2(freq_i)       = out_fft(fPt)/in_fft2(fPt); % output/input [rad/s^2]/[rad]
            H_Y2(freq_i)       = H_Y2(freq_i)*Gdi*Gdi;      % integrate twice -> output/input [rad]/[rad]
            H_Y2_full(freq_i)  = out_fft(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad]
            H_Y2_full(freq_i)  = H_Y2_full(freq_i)*Gdi*Gdi; % integrate twice -> output/input [rad]/[rad]
            H_Y2_joint(freq_i) = in_fft2(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad]
            

        elseif commaded_idx == 3 % Exciting the system using the elbow joint
            H_Y3(freq_i)       = out_fft(fPt)/in_fft2(fPt); % output/input [rad/s^2]/[rad]
            H_Y3(freq_i)       = H_Y3(freq_i)*Gdi*Gdi;      % integrate twice -> output/input [rad]/[rad]
            H_Y3_full(freq_i)  = out_fft(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad]
            H_Y3_full(freq_i)  = H_Y3_full(freq_i)*Gdi*Gdi; % integrate twice -> output/input [rad]/[rad]
            H_Y3_joint(freq_i) = in_fft2(fPt)/in_fft(fPt);  % output/input(desired position) [rad/s^2]/[rad]
         
        end

        % plot_process = 1;
        % if plot_process == 1
        %     figure(100)
        %     subplot(2,2,1)
        %     plot(t_des(index_sin_fft),INPUT(index_sin_fft),'b','LineWidth',2)
        %     hold on
        %     plot(t_des(index_sin_fft),q_temp,'r')
        %     hold off
        %     grid on
        %     xlabel('time [s]')
        %     ylabel('joint position [rad]')
        %     title(['Input Frequency = ' num2str(freq_range(freq_i)) ' Hz'])
        % 
        %     subplot(2,2,3)
        %     plot(t_des(index_sin_fft),OUTPUT(index_sin_fft),'b','LineWidth',2)
        %     grid on
        %     xlabel('time [s]')
        %     ylabel('base acceleration [rad/s^2]')
        % 
        %     subplot(2,2,2)
        %     plot(f_fft,2*abs(in_fft))
        %     xlabel('Frequency [Hz]')
        %     ylabel('FFT Mag. Input')
        %     hold on
        %     plot(f_fft(fPt),2*abs(in_fft(fPt)),'ko')
        %     xlim([1 60])
        %     hold off
        %     grid on
        % 
        %     subplot(2,2,4)
        %     plot(f_fft,2*abs(out_fft))
        %     xlabel('Frequency [Hz]')
        %     ylabel('FFT Mag. Onput')
        %     hold on
        %     plot(f_fft(fPt),2*abs(out_fft(fPt)),'ko')
        %     xlim([1 60])
        %     hold off
        %     grid on
        % 
        %     pause(0.0001);
        % end

    end

       
figure(k)
subplot(2,1,1)
semilogx(freq_range,mag2db(abs(H_Y1_full)),'LineWidth',2)

    yline(-3,'r-.','-3 dB','LineWidth',2);
    ylabel('Mag [dB]')
    grid on
    title('End-effector Displacement to Base Desired')
    %legend('Pose1','Pose2','Pose3','Pose4','Pose5','Pose6','Pose7','Pose8','Pose9','-3db')

subplot(2,1,2)
semilogx(freq_range,unwrap(angle(H_Y1_full))*180/pi,'LineWidth',2)
    yline(-3,'r-.','-3 dB','LineWidth',2);
    ylabel('Mag [dB]')
    grid on
    title('End-effector Displacement to Base Desired')
    %legend('Pose1','Pose2','Pose3','Pose4','Pose5','Pose6','Pose7','Pose8','Pose9','-3db')
end




   

%% Save joint position vector
% joint_positions_filename = ['joint_positions_matrix_',date];
% joint_positions_filename = replace(joint_positions_filename,'-','_');
% save(fullfile(frf_dir_string,joint_positions_filename),'joint_positions')