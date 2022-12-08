clear;
Rat=2; % what rat data do we analyze? 1 for, 11 for Bubble 2nd dataset, 2 for Roger, 3 for Coconut, 4 for Rolf
identical=1; % plot control (identical shifts to all) or shifts by modules (0 for modules-wise results , 1 for identical shift results)
speed_thresh=3; % speed threshold used in simulations
L_D=0:1; % dark and light (0 for dark, 1 for light)
tau=100; % time of exp kernel we analyze [ms]

Modules=1:3;
if Rat==3 % Coconut has only 2 modules
    Modules=1:2;
end
n_uni_pairs = factorial(numel(Modules))/2; % # of uni module pairs

ori_dir=pwd; % original directory
%% Shifts and trials used
spatial_shift=0:1:25; % shifts that were used
tr=1:30; % used trials

per_autocorr_decrease = .15; % threshold value to which autocorrelation should drop to signal stop accumulating variance
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
%% Running on dark and light
for ld=L_D % run on dark and light data
    %% Getting speed
    speed=data.task(ld+1).tracking.MY_speed; % recorded speed in this experiment
    speed_idc=find(speed<speed_thresh); % indices below speed threshold to be removed
    tt = numel(speed); % # of total time steps for this illimination condition
    
    
    U_to_M_mean = zeros(numel(Modules),numel(spatial_shift)); % will contain Uni to Multi mean distances for each spatial shift
    U_to_M_SEM = U_to_M_mean; % will contain corresponding SEM    
    
    U_to_U_mean = zeros(n_uni_pairs,numel(spatial_shift)); % will contain Uni to Uni mean distances for each spatial shift
    U_to_U_SEM = U_to_U_mean;  % will contain corresponding SEM
    
    M_to_T_mean_and_SEM = zeros(2,numel(spatial_shift)); % will contain Multi to True mean distances (MAE) for each spatial shift
    
    
    shift_counter=1;
    for j=spatial_shift % running on shifts
        %% Files for the j'th spatial shift
        rm_times=[]; % will contian all times to remove for the j'th shift (from all trials)
        
        Unis_to_Multi = zeros(numel(Modules),tt); % will contain the average distances over time of uni to multi 
        Unis_to_Unis = zeros(n_uni_pairs,tt); % will contain the average distances over time of uni to uni 
        Multi_Err = zeros(1,tt); % will contain the average distance over time of multi to true 
               
        t_count=0; % trials counter
        for i=tr % running on all trials of this spatial shift
            
            % Get the relevant distances vs time for this simulation
            [Distances_uni_to_multi, Distances_between_unis, multi_err, Remove] = Uni_distances_measure_cluster(tau,Rat,ld,j,i,identical);
            %% Summing Values vs time
            Unis_to_Multi = Unis_to_Multi + Distances_uni_to_multi; % adding the next trial, vs time
            Unis_to_Unis = Unis_to_Unis + Distances_between_unis; % adding the next trial, vs time
            Multi_Err = Multi_Err + multi_err; % adding the next trial, vs time
            
            rm_times = unique([rm_times,Remove]); % times to be eventually removed
            t_count = t_count+1;
        end
        rm_times=unique([rm_times,speed_idc']); % times below speed threshold and without decoding results from all trials
        Unis_to_Multi(:,rm_times) = []; % remove times that should be removed
        Unis_to_Unis(:,rm_times) = []; % remove times that should be removed
        Multi_Err(rm_times) = []; % remove times that should be removed        
        
        Unis_to_Multi = Unis_to_Multi./t_count; % Average vs time over trials stochaticity of shifts and dilution of spikes (have distances vs time)
        Unis_to_Unis = Unis_to_Unis./t_count; % Average vs time over trials stochaticity of shifts and dilution of spikes (have distances vs time)
        Multi_Err = Multi_Err./t_count; % Average  vs time over trials stochaticity of shifts and dilution of spikes (have distances vs time)
        
        
        U_to_M_mean(:,shift_counter) = mean(Unis_to_Multi,2); % mean distances between uni-mles to multi-mle for the j'th spatial shift
        for k=1:size(U_to_M_SEM,1) % running on modules to generate corresponding SEMs
            U_to_M_SEM(k,shift_counter) = SEM_time_series(Unis_to_Multi(k,:),per_autocorr_decrease); % corresponding SEMs
        end
        
        U_to_U_mean(:,shift_counter) = mean(Unis_to_Unis,2); % mean distances between uni-mles pairs for the j'th spatial shift
        for k=1:size(U_to_U_SEM,1) % running on pairs of modules to generate corresponding SEMs
            U_to_U_SEM(k,shift_counter) = SEM_time_series(Unis_to_Unis(k,:),per_autocorr_decrease); % corresponding SEMs
        end
        
        M_to_T_mean_and_SEM(1,shift_counter) = mean(Multi_Err); % MAE for the j'th spatial shift
        M_to_T_mean_and_SEM(2,shift_counter) = SEM_time_series(Multi_Err,per_autocorr_decrease); % corresponding SEMs
        
         
        shift_counter=shift_counter+1;
    end
    %% Placing results in a structure
    if ld==0 % dark
        Dark.U_to_M_mean = U_to_M_mean;
        Dark.U_to_M_SEM = U_to_M_SEM;
        Dark.U_to_U_mean = U_to_U_mean;
        Dark.U_to_U_SEM = U_to_U_SEM;
        Dark.M_to_T_mean_and_SEM = M_to_T_mean_and_SEM;
    
    elseif ld==1 % light
        Light.U_to_M_mean = U_to_M_mean;
        Light.U_to_M_SEM = U_to_M_SEM;
        Light.U_to_U_mean = U_to_U_mean;
        Light.U_to_U_SEM = U_to_U_SEM;
        Light.M_to_T_mean_and_SEM = M_to_T_mean_and_SEM;
    end
end

cd(ori_dir);
%% Saving
name_d=sprintf('Dark,Tau=%d.mat',tau);
name_l=sprintf('Light,Tau=%d.mat',tau);
save(name_d,'Dark');
save(name_l,'Light');