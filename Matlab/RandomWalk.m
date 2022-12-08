function [ X_R_traj, Y_R_traj] = RandomWalk(Rat,L_D,Lim,Tt,cluster)

% Generating random walk based on time and speed from actual trajectory in
% variying artificial enivoemnet size
% Inputs:
% Rat = 1 for Bubble, 2 for Roger
% L_D = 0 for dark and 1 for light
% Lim = arena radi size, [cm]

ori_dir=pwd;

a=60; % baseline arena radi in which good coverage is acheived, in larger arena use larger speed and vice-versa
%% Load real trajctory speeds and time
if Rat==1 % Bubble
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Rat_Bubble_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data'); % when run locally
    end
elseif Rat==11 % Bubble 2nd dataset
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2/Rat_Bubble_data 2'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data 2'); % when run locally
    end
elseif Rat==2 % Roger
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger/Rat_Roger_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Roger_data'); % when run locally
    end
elseif Rat==3 % Coconut
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut/Rat_Coconut_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Coconut_data'); % when run locally
    end
elseif Rat==4 % Rolf
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf/Rat_Rolf_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Rolf_data'); % when run locally
    end
end

load('data.mat');
cd(ori_dir);
%% Recorded speed
speed=data.task(L_D+1).tracking.MY_speed; % recorded speed

p_stay_direction=.9; % probability for speed to stay in the same direction
p_stay_magni=.5; % probability for speed to stay in the same magnitude, starts with uniform .5
%% Recorded time
Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
start_time =double( data.task(L_D+1).start);
Time = Time-start_time;
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
% Tt = round(numel(Time)); % # of time points
%% Generating random walk
X_R_traj=zeros(1,Tt); % will contain the x coordinates of the random walk
Y_R_traj=zeros(1,Tt); % will contain the y coordinates of the random walk

X_R_traj(1)=0; % start from the origin
Y_R_traj(1)=0; % start from the origin

Direction = zeros(2); % will contain the movemoent direction for the previous and current steps
Direction(:,1) = sign(rand(1,2)-0.5); % sample initial direction of movements
    
V = datasample(speed,2); % initial speed magnitudes, sample absoulute value of speed from experiment (overall speed)
V=V./sqrt(2); % speed in x and y directions (need to divide by sqrt(2) )

for i=2:Tt % running on time stpes
  
    %% Change speed magnitude?
    change_magni = rand(2,1); % do we change speed magnitude? 
    idx_change = find( change_magni>p_stay_magni);
    
    for cc=idx_change % running on coordinates which need to change speed magnitude
        V(cc) = datasample(speed,1)/sqrt(2); % speed vector
    end
    
    if norm(V)<3*(Lim/a)^2 % if the speed absoulte value is below 3 cm/sec than we change it with Pr=1
        p_stay_magni=0;
    end
    if norm(V)>=3*(Lim/a)^2 % if the speed absoulte value is above 3 cm/sec than we keep it with Pr=0.8
        p_stay_magni=0.8;
    end
    
    %% Setting directions of movement
    change_direc = rand(2,1); % changing direction probablities
    idx_change = find( change_direc>p_stay_direction); % coordinates to change direction in
    idx_stay = find( change_direc<=p_stay_direction); % coordinates to leave direction as is
    
    change_direc(idx_change) = -1;  % prbablity to change direction
    change_direc(idx_stay) = 1; % prbablity to change direction
    
    Direction(:,2) = Direction(:,1) .*change_direc; % updating direction based on previous direction
    %% Update position
    X_R_traj(i) = X_R_traj(i-1) + Direction(1,2)*V(1)*dt; % position in x [cm]
    Y_R_traj(i) = Y_R_traj(i-1) + Direction(2,2)*V(2)*dt; % position in y [cm]
    
    radi = sqrt( X_R_traj(i)^2 + Y_R_traj(i)^2 ); % radi of current location
    
    while radi>=Lim-1 % if we exceed the bounadary of the arena we need to do update positions again
        
        change_magni = rand(2,1); % do we change speed magnitude? 
        idx_change = find( change_magni>p_stay_magni);
    
        for cc=idx_change % running on coordinates which need to change speed magnitude
            V(cc) = datasample(speed,1)/sqrt(2);
        end
    
        
        change_direc = rand(2,1); % changing direction probablities
        idx_change = find( change_direc>p_stay_direction); % coordinates to change direction in
        idx_stay = find( change_direc<=p_stay_direction); % coordinates to leave direction as is

        change_direc(idx_change) = -1;  % prbablity to change direction
        change_direc(idx_stay) = 1; % prbablity to change direction

        Direction(:,2) = Direction(:,1) .*change_direc; % updating direction based on previous direction
    
        
        X_R_traj(i) = X_R_traj(i-1) + Direction(1,2)*V(1)*dt; % position in x
        Y_R_traj(i) = Y_R_traj(i-1) + Direction(2,2)*V(2)*dt; % position in y
    
        radi = sqrt( X_R_traj(i)^2 + Y_R_traj(i)^2 ); % radi of current location
    end
    
    Direction(:,1) = Direction(:,2); % current direction is now the previous direction
end

end