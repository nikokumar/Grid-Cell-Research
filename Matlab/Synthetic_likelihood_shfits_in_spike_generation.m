function [] = Synthetic_likelihood_shfits_in_spike_generation(Rat,L_D,max_shift,tr)

% Use idealized neurons 
% Apply spatial shifts to rate maps before generating spike trains and
% then use the original rate maps during decoding.
% (1) Regular trimmed identical spatial shifts and (2)
% untrimmed identical spatial shifts
% and measure the likelihood and error

cluster=0; % run on cluster? set to 1
ori_dir=pwd; % current directory
identical_spatial_shift=0; % identical spatial shifts? set to 1

if L_D==0
    Na='Dark';
elseif L_D==1
    Na='Light';
end
%% Check if result already exist on cluster
if cluster==1    

    cd(Na);
    r_dir=pwd;
    
    name = sprintf('Shifted_spikes_max_shift=%d,tr=%d.mat',max_shift,tr); % name of file
    if identical_spatial_shift==1
        name = sprintf('Identical_Shifted_spikes_max_shift=%d,tr=%d.mat',max_shift,tr); % name of file
    end
    
    if exist(name,'file')~=0 % if this file exists we don't run the script
        return;
    end
    cd(ori_dir); % going back to the directory we came from
end
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
end
% rat_dir=pwd;

load('data.mat'); % loading the data file
Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
Tt = round(numel(Time)); % # of time points

rounding=0; % [cm], we will work in 1cm unit resulotion

X_traj =roundn( data.task(L_D+1).tracking.x, rounding); % x location of the animal's trajectory [cm]
% rounded to 1cm resulotion
Y_traj = roundn( data.task(L_D+1).tracking.y, rounding); % y location of the animal's trajectory [cm]
% rounded to 1cm resulotion


load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
FN=fieldnames(Neurons_sorted_according_to_modules); % names of fields
Modules = numel(FN); % # of modules

N=zeros(1,Modules); % will contain the # of neurons in each module

for nn=1:Modules
    N(nn) = numel(Neurons_sorted_according_to_modules.(FN{nn})); % # of neurons in the nn'th module
end

cd(ori_dir); % go back to original directory where we came from
%% Allocation to modules
% Corresponding allocation to modules of the neurons which their spikes are
% in ST file and tuning curve 
m1_id = 1:N(1);
m2_id = m1_id(end)+1:m1_id(end)+N(2);
m3_id = m1_id(end)+N(2)+1:sum(N);
%% Generate synthetic tuning curves and spike trains for shifted trajectory
% same tuning curves and spikes for all trials (which are only for spatial shifts)
[TCs,STs] = Sythetic_TC_and_poiss_spikes(Rat,X_traj,Y_traj,Modules,N,dt,Lim, identical_spatial_shift,max_shift,tr);

% TCs are the tuning curves used for decoding
% STs are the simulated spike trinas generated from shifted tuning curves
% spike trains

ST = [STs{1} ; STs{2}; STs{3}]; % spike trains of all neurons
clear STs;

Tuning_Curves = zeros(2*Lim+1,2*Lim+1,sum(N));
Tuning_Curves(:,:,m1_id)=TCs{1};
Tuning_Curves(:,:,m2_id)=TCs{2};
Tuning_Curves(:,:,m3_id)=TCs{3};

TC=zeros(L,L,size(Tuning_Curves,3)); % embedding inside arena with margins
TC(margin+1:margin+1+2*Lim,margin+1:margin+1+2*Lim,:)=Tuning_Curves; % embedding
Tuning_Curves=TC; clear TC; % renaming

clear TCs;
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
Tlength=Tt; % # of time iterations we analyze

Norms=zeros(1,Tlength); % will contain normalization factors
MLE=zeros(2,Tlength); % will contain the decoder's MLE values
ERR=zeros(1,Tlength); % will contain the decoder's MSE
R_distance=zeros(1,Tlength); % will contain the difference between estimated radi and true radi

probability=ones(L,L);
probability(idc)=0; % cant be outside of circled arena
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
        
        
        for k=1:sum(N) % running on all used neurons
            
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
            plot(X_R_traj(j)+Lim+margin+1,Y_R_traj(j)+Lim+margin+1,'o','MarkerEdgeColor','k','MarkerSize',30,'LineWidth',4); % real location
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
Simulation_Results.Norms=Norms; % log normalization factor
Simulation_Results.MLE=MLE; % maximum likelihood estimate
Simulation_Results.ERR=ERR; % absolute error
Simulation_Results.R_distance=R_distance;
Simulation_Results.continue_counter=continue_counter; % # of skipped timestpes
Simulation_Results.im_counter=im_counter; % # of imaginary elements
%% Saving
cd(r_dir);
save(name,'Simulation_Results');
end