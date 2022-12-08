% calulating and smoothing speeds
clear;

Rat=2; % which rat?
if Rat==1 % Bubble
    path='/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data';
elseif Rat==11 % Bubble 2nd dataset
    path='/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data 2';
elseif Rat==2 % Roger
    path='/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Roger_data';
elseif Rat==3 % Coconut
    path='/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Coconut_data';
elseif Rat==4 % Rolf
    path='/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Rolf_data';
end

load(fullfile(path,'data.mat')); % loading the dataset

for L_D=0:1 % running on dark and light conditions
    
    
    Time = data.task(L_D+1).tracking.timestamp - double(data.task(L_D+1).start); % corresponding time [sec]
    dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]

%     dt = 1/double(data.task(L_D+1).tracking_fs); % sampling frequency [Hz]

    x=data.task(L_D+1).tracking.x; % x position
    y=data.task(L_D+1).tracking.y; % y position
    
    % Using matlab gradient (1 dt difference)
    vx = gradient(x,dt); % x speed vs time [cm/sec]
    vy = gradient(y,dt); % y speed vs time [cm/sec]

    V = sqrt( (vx.^2) + (vy.^2)); % velocity absolute value [cm/sec] vs time
    %% Smoothing
    % Unclear how exactly to perform smooting, smoothing with std=200ms
    [V,~] = smoothdata(V,'gaussian',120); % correspond to 200ms standart deviation  (since smoothdata time_win=1/5 of total time window)

    %% Allocating
    if L_D==0
        V_dark=V;
    elseif L_D==1
        V_light=V;
    end
end