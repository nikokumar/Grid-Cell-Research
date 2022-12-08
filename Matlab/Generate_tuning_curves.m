close all; clear;
%% Setting the variance of the gaussian kernel
% sig=3^2; % Gaussian variance [cm^2] - ORI
sig=5^2; % Gaussian variance [cm^2]
Sigma=[sig, 0; 0, sig]; % Covariance matrix [cm^2]
ld=1; % set =0 for using dark data, and =1 for using light data
speed_threshold=3; % threshold for speed [cm/sec],  using >= this speed
%% New datasets
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

load(fullfile(path,'data.mat')); % loading the data
%% Thresholding speed
speed=data.task(ld+1).tracking.MY_speed; % my light/dark trajectory speed
speed_idc = find(speed>=speed_threshold);

x_Pos=data.task(ld+1).tracking.x; % recorded x position
y_Pos=data.task(ld+1).tracking.y; % recorded y position

Time = data.task(ld+1).tracking.timestamp; % corresponding time [sec]
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]

x_Pos = x_Pos(speed_idc); % thresholding speed, x tracking
y_Pos = y_Pos(speed_idc); % thresholding speed, y tracking
Time = Time(speed_idc); % thresholding speed, time
%% The arena
rounding=0; % [cm], we will work in 1cm unit resulotion
dr=10^rounding; % spatial resulotion [cm]

Lim=75; % [cm] radius of the circled arena
L = numel(-Lim:dr:Lim);  % size of arena [pixels] = [cm] resulotion

[X,Y] = meshgrid(-Lim:dr:Lim, -Lim:dr:Lim);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
%% Generating time map (same for all neurons)
Time_map = zeros(L, L); % will contain the time map of the arena

for k=1:numel(Time) % running on relevant time points

    gaussian = mvnpdf([X(:) Y(:)],[x_Pos(k),y_Pos(k)],Sigma);
    gaussian = reshape(gaussian, size(X,1), size(X,2));
    
    Time_map = Time_map+gaussian; % adding a gaussian around current location to time map
end

idc = sub2ind(size(Time_map), out_y,out_x); % linear indices of points outside of circled arena
%% Moving on to the spike maps
if ld==1 % if we generate tuning curves from light
    load(fullfile(path,'ST_L.mat')); % loading light spike trains
elseif ld==0 % if we generate tuning curves from dark
    load(fullfile(path,'ST_D.mat')); % loading dark spike trains
end
ST = ST(:,speed_idc); % thresholding speed, spike trains
N=size(ST,1); % # of grid cells (units)

Tuning_Curves = zeros(L, L, N); % will contain the tuning curves for each neuron

for n=1:N % running on all grid cells to genertae their tuning curves
    
    rate_map = zeros(L, L); % will contain spike map of the n'th neuron
       
    st=ST(n,:); % spike train of the n'th neuron
    spiking_indices=find(st>0); % # of indices (times) where the neuron fired
    
    for ss=1:numel(spiking_indices) % running on all the times the neuron fired
        sp=st(spiking_indices(ss)); % # of emitted spikes
        
        x_pos=x_Pos(spiking_indices(ss)); % x location of the spike
        y_pos=y_Pos(spiking_indices(ss)); % y location of the spike
        
        for ll=1:sp % running and setting gaussian the same times as the # of spikes emitted (usually 1)
            gaussian = mvnpdf([X(:) Y(:)],[x_pos,y_pos],Sigma);
            gaussian = reshape(gaussian, size(X,1), size(X,2));
        
            rate_map = rate_map+gaussian; % updating the spike map of the n'th neuron
        end        
    end
    
    tc=(1/dt).*(rate_map./Time_map); % tuning curve of the n'th neuron [Hz]
    tc(idc)=0; % points outside of circle have never been visited thus 0 firing rate
    
    Tuning_Curves(:,:,n)=tc; % registering it
end  