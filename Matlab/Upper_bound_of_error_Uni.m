% Load the multi module mle and pick a location at random within the
% circular unit cell for an upper bound of error, for the 0 shifts case

% There is no variation of spike dilution,
Rat=4; % 1 for Bubble, 11 for Bubble 2nd dataset, 2 for Roger, 3 for Coconut, 4 for Rolf

ori_dir = pwd; % original directory
tr = 1:30; % trials

L_D=0:1; % 0 for dark and 1 for light
Tau=100; % exp time we worked with [ms]
%% Arena
Lim=75; % arena radius [cm]
dr=1; % spatial resuotion [cm]
speed_thresh=3; % speed threshold used in simulations

L = numel(-Lim:dr:Lim);  % size of arena+margin [pixels] = [cm] resulotion
[X,Y] = meshgrid(-Lim:dr:Lim, -Lim:dr:Lim);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena
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
%% Uni - module circle unit cell radius
% Make sure we use the correct radius for each rat exactly as used in the decoding
if Rat==1 || Rat==11 % Bubble and Bubble 2nd dataset
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Exp Kernel/Multi-module with identical spatial shifts');
    if Rat==11
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2/Exp Kernel/Multi-module with identical spatial shifts');
    end
    r_real_loc = [21,28,40];
%     if module==1
%         used_neurons=m1_id'; % neurons that belong to the 1st module
%         r_real_loc = 21; % radi of a circle unit cell [cm]
%     elseif module==2
%         used_neurons=m2_id'; % neurons that belong to the 2nd module
%         r_real_loc = 28; % radi of a circle unit cell [cm]
%     elseif module==3
%         used_neurons=m3_id'; % neurons that belong to the 3rd module
%         r_real_loc = 40; % radi of a circle unit cell [cm]
%     end
end

if Rat==2 % Roger
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger/Exp Kernel/Multi-module with identical spatial shifts');
    r_real_loc = [25,40,52];
%     if module==1
%         used_neurons=m1_id'; % neurons that belong to the 1st module
%         r_real_loc = 25; % radi of a circle unit cell [cm]
%     elseif module==2
%         used_neurons=m2_id'; % neurons that belong to the 2nd module
%         r_real_loc = 40; % radi of a circle unit cell [cm]
%     elseif module==3
%         used_neurons=m3_id'; % neurons that belong to the 3rd module
%         r_real_loc = 52; % radi of a circle unit cell [cm]
%     end
end

if Rat==3 % Coconut (only 2 modules)
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut/Exp Kernel/Multi-module with identical spatial shifts');
    r_real_loc = [26,45];
%     if module==1
%         used_neurons=m1_id'; % neurons that belong to the 1st module
%         r_real_loc = 26; % radi of a circle unit cell [cm]
%     elseif module==2
%         used_neurons=m2_id'; % neurons that belong to the 2nd module
%         r_real_loc = 45; % radi of a circle unit cell [cm]
%     end
end

if Rat==4 % Rolf
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf/Exp Kernel/Multi-module with identical spatial shifts');
    r_real_loc = [22,28,43];
%     if module==1
%         used_neurons=m1_id'; % neurons that belong to the 1st module
%         r_real_loc = 22; % radi of a circle unit cell [cm]
%     elseif module==2
%         used_neurons=m2_id'; % neurons that belong to the 2nd module
%         r_real_loc = 28; % radi of a circle unit cell [cm]
%     elseif module==3
%         used_neurons=m3_id'; % neurons that belong to the 3rd module
%         r_real_loc = 43; % radi of a circle unit cell [cm]
%     end
end
%% Get the multi module mle trajectory
ld_dir = pwd; % parent directory of light/dark results we use
Guessed_Error_Uni = cell(1,numel(L_D));
for ld=L_D   
    speed=data.task(ld+1).tracking.MY_speed; % recorded speed in this experiment
    speed_idc=find(speed<speed_thresh); % indices below speed threshold to be removed

    cd(ld_dir);
    if ld==0
        Na='Dark';
    elseif ld==1
        Na='Light';
    end
    cd(Na);
    
    load([Na,sprintf(':Tau=%dms,identical_spatial_shift=0,tr=1.mat',Tau)]); % we look at 0 spatial shift, and all trials are the same so we pick the first trial
    multi_mle = Simulation_Results.MLE; % multi module mle
    % need to set once multi_mle to be in the range [-Lim,Lim] like X and Y trajectories
    multi_mle = multi_mle-Lim; % there was not margin in the exp decoder
        
    Rel_times = Simulation_Results.Rel_times; % time we started to decode from
    rm_times = Simulation_Results.rm_times; % times to remove where no multi mle was generated
    rm_times=unique([rm_times,speed_idc']);
    rm_times(rm_times==0)=[]; % remove 0 time

    guessed_error_uni = zeros(numel(r_real_loc),size(multi_mle,2),numel(tr)); % will contain the guessed error for every module vs time and trial for this illumination

    %% Running on time steps and drawing random position within the inner circle each time step
    for j = Rel_times+1:size(multi_mle,2) % running on times where we have decoded position
        
        if sum(j==rm_times)==0 % as long as it's a time point that sholdn't been removed
        
            XX = X - multi_mle(1,j);
            YY = Y - multi_mle(2,j);
    
            r_sq = XX.^2 + YY.^2;

            for t=tr % running on independent trials
    
                for m=1:numel(r_real_loc) % running on modules
        
                    Log_P = ones(L); % will contain possible locations in the unit cell for the m'th module
                    
                    [local_out_y,local_out_x]=find(r_sq>=r_real_loc(m)^2); % coordinates to set to 0 firing rate since they are outside of circle
                    local_idc = sub2ind(size(X), local_out_y,local_out_x);
        
                    Log_P(local_idc)=0; % cropping around the unit cell
                    Log_P(idc)=0; % cropping what's outside the arena
                    
                    [draw_y,draw_x] = find(Log_P==1); % possible positions to sample randomly from the inner circle
                    
                    num_of_points = numel(draw_x); % # of points we can sample from
                    ind = round((num_of_points-1).*rand) + 1; % index from 1 to # of points we can sample from
                    
                    Uni_guessed_x = draw_x(ind) -Lim-1; % x guessed for the j'th time step and the m'th module
                    Uni_guessed_y = draw_y(ind) -Lim-1; % y guessed for the j'th time step and the m'th module
                    
                    guessed_error_uni(m,j,t) = sqrt( (Uni_guessed_x-multi_mle(1,j))^2 +(Uni_guessed_y-multi_mle(2,j))^2 );
                end
            end
        end
    end

    guessed_error_uni(:,rm_times,:) = []; % removing rm times

    guessed_error_uni = mean(guessed_error_uni,2); % average over time
    guessed_error_uni = squeeze(guessed_error_uni); % squeeze to mean for each module and trial

    Guessed_Error_Uni{ld+1} = guessed_error_uni; % registering the results for this illuminatation
end

%% Saving
cd(ori_dir);
cd('Uni');
save Guessed_Error_Uni Guessed_Error_Uni;