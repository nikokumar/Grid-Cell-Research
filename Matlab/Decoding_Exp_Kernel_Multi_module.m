function [ ] = Decoding_Exp_Kernel_Multi_module( Rat,L_D,Tau, max_shift,tr)
% Decoding Tergeir data set (Light/Dark) using estimated tuning curves and
% recorded spikes
% Use Biological estimated tuning curves and Biological spikes to decode


ori_dir=pwd;
speed_threshold=3; % diluting spikes trains  - but according to speed thresohld
%%
cluster=1; % run on cluster? set to 1

Rel_Tau=10; % how many Tau's back is it relevant to look at on the spike trains?
% Beyond that spikes are weighted by a 0 instead of the kernal
%% Check if result already exist
if L_D==0
    Na='Dark';
elseif L_D==1
    Na='Light';
end

if cluster==1
    cd(Na);
end

name = [Na,sprintf(':Tau=%dms,identical_spatial_shift=%d,tr=%d.mat',Tau,max_shift,tr)]; % name of file
if exist(name,'file')~=0 % if this file exists we don't run the script
    return;
end
%% Set conditions
identical_time_shift=0; % do we present an indentical time shift to all neurons? if yes set to 1
% this should check if the likelihood doesn't decrease
time_shift=0; % do we present time shift to distinct modules? if yes set to 1

identical_spatial_shift=1; % do we present an indentical spatial shift to all neurons? if yes set to 1
% this should check if the likelihood doesn't decrease
spatial_shift=0; % do we present a spaital shift to distinct modules? if yes set to 1

Poisson_spikes=0; % use generated artificial poisson spikes? if yes set to 1

dilute_spikes=0; % do we filter spikes, if yes set to 1

permute_spikes=0; % do we permute the spike trains (for each neuron)? if yes set to 1

Dilute_dark_AND_light=0; % do we dilute both dark and light spikes ? if yes set to 1


percentage_of_neurons=100; % the percentage of used neurons [0-100] - check that MSE increases when using less neurons

groups_of_diff_neuron=0; % if we decode with identical sized group of neurons, but differnet neurons

upper_bound_of_error=0; % check what is upper bound of decoding error? if yes set to 1
%% Loading data strcut and extracting trajectory for reference
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
rat_dir=pwd;

load('data.mat'); % loading the data file

rounding=0; % [cm], we will work in 1cm unit resulotion

X_traj =roundn( data.task(L_D+1).tracking.x, rounding); % x location of the animal's trajectory [cm]
% rounded to 1cm resulotion
Y_traj = roundn( data.task(L_D+1).tracking.y, rounding); % y location of the animal's trajectory [cm]
% rounded to 1cm resulotion
    
Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
start_time =double( data.task(L_D+1).start);
    
Time = Time-start_time;

dr=10^rounding; % spatial resulotion [cm]

Lim=75; % [cm] radius of the circled arena
margin=0; % [cm] margin addition to the radius for not losing (and reflecting) probability leaks from outside the real arena back in
% Don't need margin in this case since there is no leak of probability

L = numel(-Lim-margin:dr:Lim+margin);  % size of arena+margin [pixels] = [cm] resulotion

[X,Y] = meshgrid(-Lim-margin:dr:Lim+margin, -Lim-margin:dr:Lim+margin);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
Tt = numel(Time); % # of time points

%% Exp kernel time constant
Tau=Tau/1000; % translate from [ms] to [sec]
Tau_indices=round(Tau/dt); % translating time constant [sec] to indices
Rel_times=round(Rel_Tau*Tau_indices); % # of time points that we look at backward
%% Load tuning curves and spikes
load('Tuning_Curves.mat'); % loading estimated tuning curves
TC=zeros(L,L,size(Tuning_Curves,3)); % embedding inside arena with margins
TC(margin+1:margin+1+2*Lim,margin+1:margin+1+2*Lim,:)=Tuning_Curves; % embedding
Tuning_Curves=TC; clear TC; % renaming

Tuning_Curves=Tuning_Curves+10^-5; % adding a bit firing rate to avoid log(0) = -inf

load('ID.mat'); % corrseponding ID's of neurons
load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
N=size(Tuning_Curves,3); % # of neurons (grid cells)

if L_D==1 % loading Light spike trains
    load('ST_L.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory
elseif L_D==0 % loading Dark spike trains
    load('ST_D.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory   
end
%% If we generate random trajectory from uniform posterior
if upper_bound_of_error==1
    [X_guessed, Y_guessed ] = Upper_bound_of_Error(Tt,Lim,rounding,tr);
end
%% Allocation to modules
m1_id=Neurons_sorted_according_to_modules.M1; % id of neuorons from first module
m2_id=Neurons_sorted_according_to_modules.M2; % id of neuorons from second module
m3_id=Neurons_sorted_according_to_modules.M3; % id of neuorons from third module

for m=1:numel(m1_id) % running on 1st module neurons
    m1_id(m)=find(m1_id(m)==ID); % translating unit number to serial neuron # (from units # to 1 to N)
end
for m=1:numel(m2_id) % running on 2nd module neurons
    m2_id(m)=find(m2_id(m)==ID);
end
for m=1:numel(m3_id) % running on 3rd module neurons
    m3_id(m)=find(m3_id(m)==ID);
end
cd(ori_dir); % go back to original directory where we came from
%% Conditions
if Dilute_dark_AND_light==1 % we use new spike trains diluted based on light and dark, and on speed
    
    if time_shift==1 || identical_time_shift==1 % can't use >0 speed threshold for temporal shift, only for spatial shifts
        speed_threshold = 0; % can't use >0 speed threshold for temporal shift, only for spatial shift
    end
    
    [ NEW_STD, NEW_STL ] = Dilute_Spikes_neuron_wise( Rat, speed_threshold,tr ); % Dliuted Dark and light spike trains
    if L_D==0 % if we decode dark
        ST=NEW_STD;
    elseif L_D==1 % if we decode light
        ST=NEW_STL;
    end
    clear NEW_STD NEW_STL;
end


if Poisson_spikes==1 % replacing by artificial poisson spikes for the analyzed trajectory
    [ST] = Poisson_Spikes_Encoding(Time,X_traj,Y_traj,tr); % generatign Poisson spike trains for given trajectory
end

if dilute_spikes==1 && L_D==1 % do we dilute the spike trains (should be diluting light spike train to match with mean f.r of dark spike trains)   
    Split_to_segments=0; % do we treat all spike train as a single segment or seperate and dilute based on speed
    % for a single segent set to 0, for segments based on speed set to 1
    
    STL=ST; % register the light spike trains   
    clear ST;
    cd(rat_dir);
    load('ST_D.mat'); % loading the dark spike train as a reference for dilution
    cd(ori_dir);
    STD=ST;
    clear ST;
  
    if Split_to_segments==0 % Dilute based on all spike train, a fixed percentage for all the spike train    
        pf=mean(mean(STD))/mean(mean(STL));
        percent_filtered=(1-pf)*100; % percentage tot filter so that mean firing of dilute light will match mean firing of dark
%         percent_filtered=24.2; % optimal percentage to filter to have same mean firing rate as darknes (no speed threshold)
%         percent_filtered=26.7117; % optimal percentage to filter to have same mean firing rate for as darkness (with 2.5 speed threshold)
        [ST] = Filter_spikes(STL,percent_filtered,tr); % diluted spikes trains
        
        clear STL;
        clear STD;
    end    
    
    if Split_to_segments==1 % Dilute based on segment according to velocity
        % when we introduce speed threshold we compare seperatly the mean firing rate
        % between dark and light for segments above the threshold and for segments below the threshold
        
        [ST] = Filter_spikes_speed_thersholded(STL,STD,speed_threshold,tr); % diluted spike train according to speed segements
    end       
end

if permute_spikes==1 % do we permute the spike trains (mean firing rate is kept fixed in this manipulation)
    modules_together=1; % do we perform the same permutation to all neurons that belong to the same module? set to 1
    % for independent permutations across neurons set to 0
    [ST] = permuted_spike_train(Rat,L_D,ST,modules_together,speed_threshold,tr); % permuted spike trains
end

if time_shift==1 % temporal shift to spike trains (for each module independently or for all modules same shift)
%     max_shift=10; % maximal shift each module can get in units of dt
    rng(tr);
    shifts=round(2.*max_shift.*rand(1,3) -max_shift);
    ST(m1_id',:)=circshift(ST(m1_id',:),[0,shifts(1)]);
    ST(m2_id',:)=circshift(ST(m2_id',:),[0,shifts(2)]);
    ST(m3_id',:)=circshift(ST(m3_id',:),[0,shifts(3)]);
end

if identical_time_shift==1
    ST=circshift(ST,[0,identical_shift]);
end

if spatial_shift==1 % sptial shift according to modules
    [Tuning_Curves]=Spatial_shifts(m1_id,m2_id,m3_id,Tuning_Curves,max_shift,tr);
    for k=1:N
        tc=Tuning_Curves(:,:,k);
        tc(idc)=0; % 0 firing rate outside the arena
        Tuning_Curves(:,:,k)=tc; 
    end
end

if identical_spatial_shift==1 % same spatial shift to all neurons
    rng(tr);
    
    shifts_x = round(2.*max_shift.*rand(1) - max_shift); % identical spatial shift to all neurons in x axis
    shifts_y = round(2.*max_shift.*rand(1) - max_shift); % identical spatial shift to all neurons in y axis
    for k=1:N
        tc=Tuning_Curves(:,:,k);
        tc=circshift(tc,[shifts_y(1),shifts_x(1)]);
        tc(idc)=0;
        
        Tuning_Curves(:,:,k)=tc;
    end
end


used_neurons=sort([m1_id',m2_id',m3_id']);
if percentage_of_neurons<100 % meanning if we use less neurons
    rng(tr); % random seed
    AA=rand(1,N);
    AA(AA<=percentage_of_neurons/100)=-1; % tag neurons to keep using
    used_neurons=find(AA==-1); % the neurons we keep for the decoding
end

if groups_of_diff_neuron==1
    [Groups] = Groups_of_neurons(N,Gr,tr); % allocation of the neurons to Gr stranger groups for the tr realization
    used_neurons = Groups(gr,:); % the gr group
end
%% MLE, Abouslute Error and Radi error
MLE=zeros(2,Tt); % will contain the decoder's MLE values
ERR=zeros(1,Tt); % will contain the decoder's MSE
R_distance=zeros(1,Tt); % will contain the difference between estimated radi and true radi

rm_times=zeros(1); % will contain the time points where probabiltiy in the unit cell is 0
rm_index=1; % counter for that
%% Generating movie
generate_movie=-1; % want to generate movie? set to 1
if generate_movie==1
    fr=10; % once every how much dt we save a frame
    movie_stsrt_time=1; % the time we start recording from
    fig=figure(1);
%     set(0, 'DefaultFigureVisible', 'off');
    vidfile=VideoWriter('Likelihood.mp4','MPEG-4');
    vidfile.FrameRate = 120/fr; % data come in 120 Hz
    open(vidfile);
end
%% Decoding likelihood from the spikes
for j=Rel_times+1:Tt % running on time steps  
    
    
    Log_P=zeros(L,L); % will contain the sum of the log of the probability for this time step              

    Tau_Spikes = ST(used_neurons,j-Rel_times:j); % spike trains of all relevant neurons during relevant backward times
    spike_count = sum(Tau_Spikes,2); % sum of spikes during this time interval
    fired = find(spike_count>0); % find which neurons fired in this time interval

    for fi=1:numel(fired) % running only on neurons that did fire during this time interval

        rf = Tuning_Curves(:,:,used_neurons(fired(fi))); % receptive field of the relevant neuron that did fire
        st=Tau_Spikes(fired(fi),:); % spike train during this time interval (of the relevant neuron that did fire)

        spike_timings = find(st>0); % timing of spikes
        num_of_spikes = st(spike_timings); % how many spikes where emitted in these times?
        weighted_spikes = num_of_spikes .* exp(- ( (Rel_times+1) - spike_timings) / Tau_indices); % weighing spikes by their timing and the # of emitted spikes
        weighted_spikes = sum(weighted_spikes); % total effect of the spikes of that neuron

        Log_P = Log_P + weighted_spikes.*log(rf); % log likelihood
    end
    
    if sum(sum(Log_P))==0 % literally 0 probability in the arena (no spikes were fired during Rel_Tau times)
       rm_times(rm_index)=j;
       rm_index=rm_index+1;
   end
    %% MLE and its error with respect to the true location
    [~, max_idx]=max(Log_P(:));
    [r,c]=ind2sub(size(Log_P),max_idx);

    if upper_bound_of_error==1 % if we just guess decoded location from a unifrom posterior
        c = X_guessed(j)+Lim+1+margin; % guessed x location
        r = Y_guessed(j)+Lim+1+margin; % guessed y location
    end

    MLE(:,j)=[c;r]; % [x;y] MLE [pixel]

    Xerror = c - (X_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 
    Yerror = r - (Y_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 

    Pixeldist=sqrt((Xerror^2)+(Yerror^2)); % [Pixels]
    ERR(j)=Pixeldist.*dr; % [cm]

    true_radius = sqrt((X_traj(j))^2 + (Y_traj(j))^2);
    estimated_radius = sqrt((c-Lim-margin)^2 + (r-Lim-margin)^2);
    R_distance(j) = abs(true_radius - estimated_radius); % distance between true radius and estimated radius
    %% Movie
    if generate_movie==1 && mod(j,fr)==0 && j>=movie_stsrt_time
        clf(fig,'reset'); % clearing the figure

        imagesc(Log_P); % posterior likelihood
        set(gca,'YDir','normal');
        hold on;
        plot(X_traj(j)+Lim+margin+1,Y_traj(j)+Lim+margin+1,'o','MarkerEdgeColor','k','MarkerSize',30,'LineWidth',4); % real location
        plot(MLE(1,j),MLE(2,j),'or','MarkerSize',10,'LineWidth',2); % plotting maximum likelihhod estimate
        % +L/2 term since location (0,0) is in the middle of array

        plot(out_y,out_x,'w.','markersize',8);
        ylabel('cm'); xlabel('cm');
        filename=sprintf('Time = %d [sec]',j*dt); % time in ms
        title(filename);
        set(gca,'fontsize',16);     

        Mo1=getframe(gcf); % movie of the first base
        writeVideo(vidfile,Mo1);
    end

end
%% Organizing Results
Simulation_Results=struct;
Simulation_Results.MLE=MLE; % maximum likelihood estimate
Simulation_Results.ERR=ERR; % absolute error
Simulation_Results.R_distance=R_distance;
Simulation_Results.Rel_Tau=Rel_Tau; % kernel time constant [sec]
Simulation_Results.Rel_times=Rel_times; % # of indices we look back at, we didn't have mle fro times 1:Rel_times

Simulation_Results.rm_times=rm_times; % times to be removed, with identity 0 probability at all locations (no spikes occured)
%% Saving
if cluster==1
    cd(Na); % getting into Light/Dark directory
end
save(name,'Simulation_Results');
end