% Generating a guessed trajectory for each trial and measure its error from
% the corresponding mle for the 0 shifts case

Rat=11; % what rat data do we analyze? 1 for Bubble, 11 for Bubble 2nd dataset, 2 for Roger, 3 for Coconut, 4 for Rolf

speed_thresh=3; % speed threshold used in simulations
L_D=0:1; % dark and light (0 for dark, 1 for light)
tr = 1:30; % trials that were used, there is variation of spike dilution for each trial

ori_dir=pwd; % original directory of upper bounds
%% Paramters of arena as in decoding results for guessing a trajectory
rounding=0; % [cm], we will work in 1cm unit resulotion
Lim=75;
margin=5; % [cm] margin addition to the radius for not losing (and reflecting) probability leaks from outside the real arena back in
%% Getting speed
if Rat==1 % Bubble
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Rat_Bubble_data'); % going into rat 1 (Bubble) data directory
elseif Rat==11 % Bubble 2nd dataset
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2/Rat_Bubble_data 2'); % going into rat 1 (Bubble) 2nd data directory
elseif Rat==2 % Roger
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger/Rat_Roger_data'); % going into rat 2 (Roger) data directory
elseif Rat==3 % Coconut
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut/Rat_Coconut_data'); % going into rat 3 (Cocnut) data directory
elseif Rat==4 % Rolf
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf/Rat_Rolf_data'); % going into rat 4 (Cocnut) data directory
end
load('data.mat'); % loading the data file
cd(ori_dir);

%% Directory where multi mle trajectories are
if Rat==1 % Bubble
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble'); % going into rat 1 (Bubble) data directory
elseif Rat==11 % Bubble 2nd dataset
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2'); % going into rat 1 (Bubble) 2nd data directory
elseif Rat==2 % Roger
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger'); % going into rat 2 (Roger) data directory
elseif Rat==3 % Coconut
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut'); % going into rat 3 (Cocnut) data directory
elseif Rat==4 % Rolf
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf'); % going into rat 4 (Cocnut) data directory
end

cd('Diltued dark AND light spikes (+by speed threshold)');
cd('Spatial shifts'); % we use 0 shift anyways so it doesn't matter if we use identical/independent shifts
shifts_dir = pwd; % parent shifts directory
%%
Guessed_error_Multi = zeros(numel(L_D),numel(tr)); % will contain the mean errors for each illumination an trials

for ld = L_D % running on dark and light
    speed=data.task(ld+1).tracking.MY_speed; % recorded speed in this experiment
    speed_idc=find(speed<speed_thresh); % indices below speed threshold to be removed

    Tt = numel(speed); % # of total time points before filtering
    %% Go into dark/light
    cd(shifts_dir);
    if ld==0
        Na='Dark';
    elseif ld==1
        Na='Light';
    end
    cd(Na); % going into results trajectory
    r_dir = pwd; % results dir where estimated trajectories are
    
    %% Get all times to remove from all trials
    rm_times = []; % will contain all time to remove from all trials
    for i=tr
        % Loading the 0 shift corresponding multi-module mle
        load([Na,sprintf('_spatial_shift=0_tr=%d.mat',i)]);

        no = Simulation_Results.Norms; % log(z)
        C=find(no==0); % due to continue counter process (spikes totally disagreed and this times step is digarded
        rm_times=[rm_times,C]; % times to remove for all trials eventually
    end

    rm_times=unique([rm_times,speed_idc']); % times below speed threshold and without decoding results from all trials
    %%
    t_counter=1;
    for i=tr % running on all trials in this illumination conditions

        %% First generate a guessed trajectory
        cd(ori_dir);
        [X_guessed, Y_guessed] = Guessed_trajectory_Multi(Tt,Lim,rounding,i); % the guessed trajectory for this trial
        X_guessed = X_guessed + Lim+1+margin; % transform to arena coordinates [margin,Lim
        Y_guessed = Y_guessed + Lim+1+margin; % transform to arena coordinates

        %% Loading the 0 shift corresponding multi-module mle
        cd(r_dir); % going into relevant results directory
        load([Na,sprintf('_spatial_shift=0_tr=%d.mat',i)]);

        mle = Simulation_Results.MLE; % the multi module mle

        g_e = sqrt( (mle(1,:)-X_guessed').^2 + (mle(2,:)-Y_guessed').^2 ); % the chance level error vs time
        g_e(rm_times)=[]; % removing times
        Guessed_error_Multi(ld+1,t_counter) = mean(g_e); % mean chance level error for this illumination and trial

        t_counter = t_counter+1;
    end
end

%% Saving
cd(ori_dir);
cd('Multi'); % multi decoding directory
save Guessed_error_Multi Guessed_error_Multi;