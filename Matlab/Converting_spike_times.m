% Converting the given timing of spikes to a time vector of the whole
% experimnet with 0 for no spikes and 1 for spikes
clear;
Rat=4; % which rat?
L_D=1; % =1 for Light conditions and 0 for Dark conditions


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

load(fullfile(path,'data.mat'));

%% Load recorded spike times
Spike_times = data.task(L_D+1).spike_timestamp; % times of all the spikes from dark
Spike_units = data.task(L_D+1).spike_cluster_id; % corresponding relevant unit that fired
%% Load recorded tracking times
Time = data.task(L_D+1).tracking.timestamp; % Time vector of trajectory tracking [sec]


Units = unique(Spike_units); % the distinct units (MUST be exactly equal to ID vector, ID vector should be sorted)
N=numel(Units); % # of neurons

dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
Tp=numel(Time); % # of time points

ST=zeros(N,Tp); % will contain the number of spikes emitted for each neuron and for all time bin

%% Setting ones for spikes in time vector

for n=1:N % running on all neurons
    
    u=find(Spike_units==Units(n)); % relevant unit spikes indices of spikes
    st=Spike_times(u); % corresponding times of spikes for the relevant neurons
    
%     st=Spike_times(n).spike_times; % times of spikes of the n'th neuron
    D=floor(st./dt); % time bins coordinates in the time vector where spikes occured
    D(D==0)=1; % for spikes happening right before even 1 dt have passed

    for k=1:numel(D) % running over spikes
        ST(n,D(k))=ST(n,D(k))+1; % add 1 for the relevant time bin
    end
end