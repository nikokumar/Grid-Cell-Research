Rat=2; % what rat data do we analyze? 1 for, 11 for Bubble 2nd dataset, 2 for Roger, 3 for Coconut, 4 for Rolf

identical=0; % plot control (identical rotations to all) or rotations by modules (0 for modules-wise results , 1 for identical shift results)
speed_thresh=3; % speed threshold used in simulations
L_D=0:1; % dark and light (0 for dark, 1 for light)

ori_dir=pwd; % original directory
%% Shifts and trials used
spatial_rotation = 0:5:90; % spatial rotation that were used [degrees]
tr = 1:30; % trials that were used

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
for ld=L_D
    %% Getting speed
    speed=data.task(ld+1).tracking.MY_speed; % recorded speed in this experiment
    speed_idc=find(speed<speed_thresh); % indices below speed threshold to be removed
    %% Iterate on results    
    Norms_mean_and_SEM = zeros(2,numel(spatial_rotation)); % will contain <log(Z)> values and corresponding SEM for each spatial rotation
    Error_mean_and_SEM = zeros(2,numel(spatial_rotation)); % will contain MAE values and corresponding SEM for each spatial rotation
    Radi_err_mean_and_SEM = zeros(2,numel(spatial_rotation)); % will contain mean radi error values and corresponding SEM for each spatial rotation
    
    rotation_counter=1; % counter of spatial rotations
    for j=spatial_rotation
        %% Locating into results directory
        if ld==1
            Na='Light';
        elseif ld==0
            Na='Dark';
        end
        cd(Na); % going into light/dark results directory  
        
        %% Files for the j'th spatial shift
        rm_times=[]; % will contain times to remove for the j'th shift, since there was no result from decoding
        
        No = zeros(1,numel(speed)); % will contain the average over trials of norms
        Eo = zeros(1,numel(speed)); % will contain the average over trials of mae
        Ro = zeros(1,numel(speed)); % will contain the average over trials of radi error

        t_count=0; % trials counter
        for i=tr % running on all trials of this spatial rotation     
            %% Loading file
            if identical==0
                name=[Na,sprintf('_spatial_rotation=%d_tr=%d.mat',j,i)]; % name of file
            elseif identical==1
                name=[Na,sprintf('_identical_spatial_rotation=%d_tr=%d.mat',j,i)]; % name of file for identical shifts case
            end
            load(name); % loading it
            %% We will remove speed at the end for all trials together vs time
            % Mark all points that we skipped during decoding because spikes
            % completly disagree, can only be if norm(t) stayed=0
            
            no = Simulation_Results.Norms; % log(z)
            C=find(no==0); % due to continue counter process (spikes totally disagreed and this times step is digarded
            rm_times=[rm_times,C]; % times to remove for all trials eventually             
            %% Summing Values vs time
            No = No + no; % sum norms vs time
            Eo = Eo + Simulation_Results.ERR; % sum error vs time
            Ro = Ro + Simulation_Results.R_distance; % sum radi error vs time
            
            t_count = t_count+1;
        end       
        %% Going back to original directory after finish iterating of trials for this spatial shift
        cd(ori_dir);
        %% Continue analysis
        rm_times=unique([rm_times,speed_idc']); % times below speed threshold and without decoding results from all trials
        No(rm_times)=[]; % remove rm times norms
        Eo(rm_times)=[]; % remove rm times absolute error
        Ro(rm_times)=[]; % remove rm times radi error
        
        % Averaging across time, end with 1 signal vs time for the j'th spatial shift
        No = No./t_count; % Average vs time norms over trials of stochatic of rotations and dilution of spikes (have norms vs time)
        Eo = Eo./t_count; % Average vs time mae over trials of stochatic of rotations and dilution of spikes (have mae vs time)
        Ro = Ro./t_count; % Average vs time radi error over trials of stochatic of rotations and dilution of spikes (have radi error vs time)
        
        
        Norms_mean_and_SEM(1,rotation_counter) = mean(No); % <log(Z)> for the j'th spatial rotation
        Norms_mean_and_SEM(2,rotation_counter) = SEM_time_series(No,per_autocorr_decrease); % corresponding SEM
        
        Error_mean_and_SEM(1,rotation_counter) = mean(Eo); % MAE for the j'th spatial rotation
        Error_mean_and_SEM(2,rotation_counter) = SEM_time_series(Eo,per_autocorr_decrease); % corresponding SEM
        
        Radi_err_mean_and_SEM(1,rotation_counter) = mean(Ro); % mean radi error for the j'th spatial rotation
        Radi_err_mean_and_SEM(2,rotation_counter) = SEM_time_series(Ro,per_autocorr_decrease); % corresponding SEM

        rotation_counter = rotation_counter+1;
    end
    %% Placing results in a structure
    if ld==1
        Light.NO = Norms_mean_and_SEM;
        Light.EO = Error_mean_and_SEM;
        Light.RO = Radi_err_mean_and_SEM;
    elseif ld==0
        Dark.NO = Norms_mean_and_SEM;
        Dark.EO = Error_mean_and_SEM;
        Dark.RO = Radi_err_mean_and_SEM;
    end
end
%% Saving
save Dark Dark;
save Light Light;