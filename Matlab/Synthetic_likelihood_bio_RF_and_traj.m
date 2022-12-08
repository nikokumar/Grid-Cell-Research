function [ ] = Synthetic_likelihood_bio_RF_and_traj( Rat,L_D, fr_fac, dil_fac, max_shift, tr )
% Generating and then decoding Poisson spikes from measured tuning curves and real trajectory 

cluster=0; % run on cluster? set to 1
rng(tr);
ori_dir=pwd;
% fr_fac = by how much do we multiply the firing rate that generates spikes
% dil_fac = by how much do we dilute the resulted spike train

spatial_shift=1; % do we present a spaital shift to distinct modules? if yes set to 1
identical_spatial_shift=0; % do we present an identical spatial shift to all neurons? if yes set to 1

temporal_shift=0; % do we present a time shift to distinct modules? if yes set to 1
identical_temporal_shift=0; % do we present an identical time shift to all neurons? if yes set to 1


if L_D==0
    Na='Dark';
elseif L_D==1
    Na='Light';
end
%% Check if result already exist on cluster
% cd(Na);
% r_dir=pwd; % where results will be registered
% 
% name = [Na,sprintf(':fr_fac=%d,dilute=%d,tr=%d.mat',round(fr_fac*10),round(dil_fac*10),tr)]; % name of file (factor 10 only for the name)
% if spatial_shift==1 || identical_spatial_shift==1 % if we use spatial shifts on Poisson spikes
%     if spatial_shift==1
%         name = [Na,sprintf(':synthetic_spatial_shift=%d,tr=%d.mat',max_shift,tr)]; % name of file (for spatial shifts case)
%     elseif identical_spatial_shift==1
%         name = [Na,sprintf(':synthetic_identical_spatial_shift=%d,tr=%d.mat',max_shift,tr)]; % name of file (for identical spatial shifts case)
%     end
% end
% if temporal_shift==1 || identical_temporal_shift==1 % if we use spatial shifts on Poisson spikes
%     if temporal_shift==1
%         name = [Na,sprintf(':synthetic_temporal_shift=%d,tr=%d.mat',max_shift,tr)]; % name of file (for temporal shifts case)
%     elseif identical_temporal_shift==1
%         name = [Na,sprintf(':synthetic_identical_temporal_shift=%d,tr=%d.mat',max_shift,tr)]; % name of file (for identical temporal shifts case)
%     end
% end
% if exist(name,'file')~=0 % if this file exists we don't run the script
%     return;
% end
%% Arena
Lim=75; % [cm] radius of the circled arena
dr=1; % resulotion [cm]
margin=5; % [cm] margin addition to the radius for not losing (and reflecting) probability leaks from outside the real arena back in

L = numel(-Lim-margin:dr:Lim+margin);  % size of arena+margin [pixels] = [cm] resulotion

[X,Y] = meshgrid(-Lim-margin:dr:Lim+margin, -Lim-margin:dr:Lim+margin);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena
%% Loading data
if Rat==1 % Bubble
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Rat_Bubble_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data'); % when run locally
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

load('data.mat'); % loading the data file
Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]

load('Tuning_Curves.mat'); % loading the tuning curves
load('ID.mat'); % corrseponding ID's of neurons
load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
N=size(Tuning_Curves,3); % # of total neurons (grid cells)

rounding=0; % [cm], we will work in 1cm unit resulotion

X_traj =roundn( data.task(L_D+1).tracking.x, rounding); % x location of the animal's trajectory [cm]
% rounded to 1cm resulotion
Y_traj = roundn( data.task(L_D+1).tracking.y, rounding); % y location of the animal's trajectory [cm]
% rounded to 1cm resulotion

cd(ori_dir); % go back to original directory
%% Generating Poisson spikes
ST=zeros(N,numel(X_traj)); % will contain the spikes vs time
 
X_traj=X_traj+Lim+1; % set indices to correspond to rf matrix
Y_traj=Y_traj+Lim+1; % set indices to correspond to rf matrix

for j=1:N % running on all neurons in that module
%     rng(1); % same spike train for all trials (which randomize only shifts)
    rf = Tuning_Curves(:,:,j); % receptive field of the j'th neuron
%     rf(idc)=0; % Sanity check - 0 firing rate outside the arena (the trajectory can't also exceed the arena boundaries)

    idx = sub2ind(size(rf),Y_traj,X_traj); % translate trajectory long indices to short index
    fr_vs_time = fr_fac.*rf(idx); % firing rate vs time

    ST(j,:) = poissrnd(fr_vs_time.*dt); % Poisson spikes       
end

X_traj=X_traj-Lim-1; % set indices back
Y_traj=Y_traj-Lim-1; % set indices back

TC=zeros(L,L,size(Tuning_Curves,3)); % embedding inside arena with margins
TC(margin+1:margin+1+2*Lim,margin+1:margin+1+2*Lim,:)=Tuning_Curves; % embedding
Tuning_Curves=TC; clear TC; % renaming
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
%% Dilute spikes?
[ST] = Dilute_Spikes_synthetic( ST,dil_fac,tr ); % Dliuted Dark and light spike trains
%% Spatial shifts?
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
%% Temporal shifts?
if temporal_shift==1 % independent temporal shift to spike trains (for each module independently)

    rng(tr);    
    shifts = 2.*max_shift.*rand(1,3) -max_shift;
    shifts = round(shifts*(1/dt)); % convert sec of shift to indicies
    
    ST(m1_id',:)=circshift(ST(m1_id',:),[0,shifts(1)]);
    ST(m2_id',:)=circshift(ST(m2_id',:),[0,shifts(2)]);
    ST(m3_id',:)=circshift(ST(m3_id',:),[0,shifts(3)]);
end

if identical_temporal_shift==1 % identical temporal shifts to spike train
    
    rng(tr);    
    shifts = 2.*max_shift.*rand(1,1) -max_shift;
    shifts = round(shifts*(1/dt)); % convert sec of shift to indicies
    ST=circshift(ST,[0,shifts(1)]); % all spike trains get the same shift
end
%% Generating Gaussian of the multi module decoder
Mu=[0,0]; % the gaussian is around the origin
D=480; % diffusion coefficient [cm^2/sec]
sigmaXwalk=D*dt; sigmaYwalk=D*dt; % variances of the Gaussian
Sigma=[sigmaXwalk,0; 0,sigmaYwalk]; % Covariance matrix

GaussianPro = mvnpdf([X(:) Y(:)],Mu,Sigma); % Probability function of the Gaussian
GaussianPro = reshape(GaussianPro, length(X), length(Y));
% GaussianPro(idc)=0; % arena is circle so there is 0 probability to be outside of it (shouldn't do anything and keep total probablity=1)

%% Multi module Decoding
Tstart=1; % iteration # we start from
Tlength=size(ST,2); % # of time iterations we analyze

Norms=zeros(1,Tlength); % will contain normalization factors
MLE=zeros(2,Tlength); % will contain the decoder's MLE values
ERR=zeros(1,Tlength); % will contain the decoder's MSE
R_distance=zeros(1,Tlength); % will contain the difference between estimated radi and true radi

probability=ones(L,L); % start with uniform probability
probability(idc)=0; % cant be outside of circled arena
% probability=zeros(L,L); % start with known location
% probability(round(L/2)-2:round(L/2)+2,round(L/2)-2:round(L/2)+2)=1;

probability(:,:,2)=probability(:,:,1); % two likelihoods for just 2 timesteps
n0=1./( (sum(sum(probability(:,:,1)))) .* (dr.^2) ); % initial normalizing factor
probability(:,:,:) = n0.*probability(:,:,:); % normalized initial likelihood
    
Current=2; % parameter that defines which of the two probabilities we update, we start with Current=2, we need to update the 2nd iteration
counter=1; % vectoer index for MSE and movie in the loop in case we dont start from first time step

continue_counter=0; % counting how many times ProSpikeswas ==zeros(L), meaning 0 probability from spikes
im_counter=0; % couring how many times the convulotion generated imaginary output
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
%% Decoding
leak_thresh=10^-10; % probabilites outside the arena smaller than this value won't be reflected to inside the arena and will be set to 0
% S=0; Su=0;
for j=Tstart:(Tstart+Tlength-1) % running on time steps   
    if Current==2
        C=1;
    else
        C=2;
    end

        ProCON=MYconv2(probability(:,:,C),GaussianPro); % The convolution between probability of previous step with the gaussian, the matrices must be squared!
        
        Image=numel(find(imag(ProCON)~=0)); % # of imaginary elements (making sure that we don't have numerical artifacts)
        if Image>0
            im_counter=im_counter+1; % recording it
            ProCON=real(ProCON); % removing imaginary parts
        end
        
        % Reflecting probability back to the arena
        Leak_candidate=ProCON(idc); % convoloved probability values outside the arena        
        
        Leak_candidate(Leak_candidate<=leak_thresh)=0; % setting smaller probabilities from threshold to 0, what left needs to be reflected back into the arena
        reflect_ind=find(Leak_candidate>0); % find indices of terms with actual probabilty needed to be reflected
        out_y_leak=out_y(reflect_ind); % y coordinates of points needs reflection
        out_x_leak=out_x(reflect_ind); % x coordinates of points need reflection
        Leak_candidate=Leak_candidate(reflect_ind); % setting tiny probabilities to 0, what left needs to be reflected back into the arena
        
        
        for l=1:numel(Leak_candidate) % running on points need reflection
            Yl=out_y_leak(l); % y coordinate of point to reflect
            Xl=out_x_leak(l); % x coordinate of point to reflect
            
            dy=Yl-L/2; % y difference between point and arena origin
            dx=Xl-L/2; % x difference between point and arena origin
            h=sqrt(dy^2 +dx^2); % distance of point from origin (should be greater than the radii 'Lim')
            
            D_prime = (Lim/h)*Lim; % desired distance (radii) of point to reflect to from the origin
            
            x_move=1; % defult is moving in the x direction
            slope=dy/dx; % change in y for each integer x step (=slope_X)
            si=sign(L/2-Xl); % direction of movemenet when changing x by (are we left or the right of the center?)
            
            if abs(slope)>1 % if going through x generates big jumps (±) in y then we will alternate x & y 
                slope=1/slope; % change in x for each integer y step (=slope_Y)
                si=sign(L/2-Yl); % direction of movemenet when changing y by (are we above or below  of the center?)
                x_move=0; % we change to move in the y firection now
            end

            while h>=D_prime % as long the distance of the point is greater than that desired
                if x_move==1 % moving 1 integer in the x direction
                    Yl = Yl+si*slope;
                    Xl = Xl+si;
                end
                if x_move==0 % moving 1 integer in the y direction
                    Yl = Yl+si;
                    Xl = Xl+si*slope;
                end
                
                dy=Yl-L/2; % y difference between point and arena origin
                dx=Xl-L/2; % x difference between point and arena origin
                h=sqrt(dy^2 +dx^2); % distance of point from origin (should be greater than the radii 'Lim')
            end
            
            % Reflecting the probabliity into the arena and eliminating what was at the outside
            Yl=round(Yl); % rounding y coordinate
            Xl=round(Xl); % rounding x coordinate
            
            ProCON(Yl,Xl)=ProCON(Yl,Xl)+ProCON(out_y_leak(l),out_x_leak(l)); % reflecting and adding the probabiltiy into the designated arena
            ProCON(out_y_leak(l),out_x_leak(l))=0; % making it zero at the original location outside
        end
               
        %
        ProSpikes=ones(L,L); % we first set the probability because of the spikes to ones
        ProSpikes(idc)=0; % can't be outside of the circled arena
        
        
        for k=1:N % running on all used neurons
            
            rf = Tuning_Curves(:,:,k); % receptive field of the k'th neuron
            spikes=ST(k,j); % # of emitted spikes for the k'th neuron and j'th time step
                        
            ProSpikes=ProSpikes.*(1/factorial(spikes)).*((rf.*dt).^spikes).*(exp(-rf.*dt)); % Poisson dist 
        end
   
        probability(:,:,Current)=ProCON.*ProSpikes; % the likelihood for j time step from all neurons
        
        if sum(sum(probability(:,:,Current)))<10^-45 % if probabilty from spikes is too small numerical errors arise
            continue_counter=continue_counter+1; % register that we skipped a time step
            counter=counter+1;
            continue; % skip this time step
        end
            
        norm = sum(sum(probability(:,:,Current))) .* (dr.^2);
        Norms(1,counter) = log(norm); % normalizing factor
        probability(:,:,Current) = probability(:,:,Current)./norm; % normalizing likelihood

        %% MLE and its error (MSE)
        Pp=probability(:,:,Current);
        [~, max_idx]=max(Pp(:));
        [r,c]=ind2sub(size(Pp),max_idx);
        
%         if upper_bound_of_error==1 % if we just guess decoded location from a unifrom posterior
%             c = X_guessed(j)+Lim+1+margin; % guessed x location
%             r = Y_guessed(j)+Lim+1+margin; % guessed y location
%         end

        MLE(:,counter)=[c;r]; % [x;y] MLE [pixel]
        
        Xerror = c - (X_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 
        Yerror = r - (Y_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 
        
        Pixeldist=sqrt((Xerror^2)+(Yerror^2)); % [Pixels]
        ERR(1,counter)=Pixeldist.*dr; % [cm]
        
        true_radius = sqrt((X_traj(j))^2 + (Y_traj(j))^2);
        estimated_radius = sqrt((c-Lim-margin)^2 + (r-Lim-margin)^2);
        R_distance(counter) = abs(true_radius - estimated_radius); % distance between true radius and estimated radius
        %% Movie
        if generate_movie==1 && mod(j,fr)==0 && j>=movie_stsrt_time
            clf(fig,'reset'); % clearing the figure
            
            imagesc(probability(:,:,Current)); % posterior likelihood
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

        %% update current state #
        if Current==2
            Current=1;
        else % meaning Current=1
            Current=2;
        end
        counter=counter+1;
end

%% Organizing results
Simulation_Results=struct;
Simulation_Results.Norms=Norms; % normalization factor
Simulation_Results.MLE=MLE; % maximum likelihood estimate
Simulation_Results.ERR=ERR; % absolute error
Simulation_Results.R_distance=R_distance;
Simulation_Results.D=D; % diffusion coefficient
Simulation_Results.continue_counter=continue_counter; % # of skipped timestpes
Simulation_Results.im_counter=im_counter; % # of imaginary elements
%% Saving
cd(r_dir);
save(name,'Simulation_Results');
end