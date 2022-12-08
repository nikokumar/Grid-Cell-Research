function [ ] = Decoding_Markov_Movie( Rat,L_D,tr ) %
% Decoding Tergeir data set (Light/Dark) using estimated tuning curves and
% recorded spikes
% Use Biological estimated tuning curves and Biological spikes to decode

ori_dir=pwd;

speed_threshold=3; % diluting spikes trains  - but according to speed thresohld
cluster=0; % run on cluster? set to 1
%% Video properties
seg=2; % the segment to simulate and make movie of
if L_D==0 % dark
    Na = 'Dark';
elseif L_D==1 % light
    Na = 'Light';
end
vid_name = [Na,sprintf(' Markov decoding %d.mp4',seg)];
b_width = .75; % width of plot box axes

prior_gap=2400; % # of time points to simulate before the begining of the segment (to approximate the markov probablity)
%% Set conditions
Dilute_dark_AND_light=1; % do we dilute both dark and light spikes ? if yes set to 1
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

speed=data.task(L_D+1).tracking.MY_speed; % recorded speed

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
margin=5; % [cm] margin addition to the radius for not losing (and reflecting) probability leaks from outside the real arena back in

L = numel(-Lim-margin:dr:Lim+margin);  % size of arena+margin [pixels] = [cm] resulotion

[X,Y] = meshgrid(-Lim-margin:dr:Lim+margin, -Lim-margin:dr:Lim+margin);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
Tt = round(numel(Time)); % # of time points
%% Load tuning curves and spikes
load('Tuning_Curves.mat'); % loading light measured tuning curves
% if d_TC==1
%     clear Tuning_Curves;
%     load('Tuning_Curves_d.mat'); % loading dark measured tuning curves
% end
TC=zeros(L,L,size(Tuning_Curves,3)); % embedding inside arena with margins
TC(margin+1:margin+1+2*Lim,margin+1:margin+1+2*Lim,:)=Tuning_Curves; % embedding
Tuning_Curves=TC; clear TC; % renaming

load('ID.mat'); % corrseponding ID's of neurons
load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
N=size(Tuning_Curves,3); % # of neurons (grid cells)

if L_D==1 % loading Light spike trains
    load('ST_L.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory
elseif L_D==0 % loading Dark spike trains
    load('ST_D.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory   
end
%% Allocation to modules
m1_id=Neurons_sorted_according_to_modules.M1; % id of neuorons from first module
m2_id=Neurons_sorted_according_to_modules.M2; % id of neuorons from second module
m3_id=Neurons_sorted_according_to_modules.M3; % id of neuorons from third module

for m=1:numel(m1_id) % running on 1st module neurons
    m1_id(m)=find(m1_id(m)==ID); % translating unit number to serial neuron # (from units # to 1 to N)
end
for m=1:numel(m2_id) % running on 1st module neurons
    m2_id(m)=find(m2_id(m)==ID);
end
for m=1:numel(m3_id) % running on 1st module neurons
    m3_id(m)=find(m3_id(m)==ID);
end

used_neurons=sort([m1_id',m2_id',m3_id']); % neurons to use for deocding

cd(ori_dir); % go back to original directory where we came from
%% Spike dilution
if Dilute_dark_AND_light==1 % we use new spike trains diluted based on light and dark, and on speed
     
    [ NEW_STD, NEW_STL ] = Dilute_Spikes_neuron_wise( Rat, speed_threshold,tr ,cluster); % Dliuted Dark and light spike trains
    if L_D==0 % if we decode dark
        ST=NEW_STD;
    elseif L_D==1 % if we decode light
        ST=NEW_STL;
    end
    clear NEW_STD NEW_STL;
end
%% Getting continuous segments which are more than speed threshold
% Use continuous segments and analyze them seperatly

min_con_time=30; % minimal consecutive time of segments [sec]
min_segment_time = round((1/dt)*min_con_time); % minimal consecutive indices


TP = ones(1,numel(speed)); % Time points template to generate continouous segments
TP(speed<speed_threshold)=0; % flag points to filter based on speed

Keep=find(TP==1); % points of interset not to be removed
Change=find(diff(Keep)~=1); % changing segment point, signals where segments ends and start

n_segments = numel(Change)+1; % total # of segments
start_points=zeros(n_segments,1); % will contain segments time starting points
end_points=zeros(n_segments,1); % will contain segments time ending points

start_points(1) = Keep(1); % starting point of first segment
for n_seg=2:n_segments % running on segments
    end_points(n_seg-1) = Keep(Change(n_seg-1)); % end point of segment
    start_points(n_seg) = Keep(Change(n_seg-1) +1); % start point of segment
end
end_points(end) = Keep(end); % end point of final segment

Below_min_consecutive = find( (end_points-start_points)<=min_segment_time); % segments shorter than threshold, to be filtered out

start_points(Below_min_consecutive)=[]; % filter starting points
end_points(Below_min_consecutive)=[]; % filter ending points
n_segments=numel(start_points); % new filtered total # of segments
%% Generating Gaussian of the multi module decoder
Mu=[0,0]; % the gaussian is around the origin
D=480; % diffusion coefficient [cm^2/sec]
sigmaXwalk=D*dt; sigmaYwalk=D*dt; % variances of the Gaussian [cm^2]
Sigma=[sigmaXwalk,0; 0,sigmaYwalk]; % Covariance matrix

GaussianPro = mvnpdf([X(:) Y(:)],Mu,Sigma); % Probability function of the Gaussian
GaussianPro = reshape(GaussianPro, length(X), length(Y));
% GaussianPro(idc)=0; % arena is circle so there is 0 probability to be outside of it (shouldn't do anything and keep total probablity=1)
%% Multi module Decoding
Tstart=1; % iteration # we start from
Tlength=Tt; % # of time iterations we analyze

Norms=zeros(1,Tlength); % will contain normalization factors
MLE=zeros(2,Tlength); % will contain the decoder's MLE values
ERR=zeros(1,Tlength); % will contain the decoder's absoule error
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
generate_movie=1; % want to generate movie? set to 1
if generate_movie==1
    fr=4; % once every how much dt we save a frame (3 correspond to 25ms, 4 correspond to 33ms and 6 correspond to 50ms)
%     movie_start_time = 1+round(st/dt); % start time for decoding and recording [sec]
%     movie_length = round(lt/dt); % total time of decoding and recording [sec]
    fig=figure(1);
%     set(0, 'DefaultFigureVisible', 'off');
    vidfile=VideoWriter(vid_name,'MPEG-4');
    vidfile.FrameRate = round(1/dt)/fr; % how many frames per sec? data comes in 120 Hz
    open(vidfile);
end

Movie_frames = start_points(seg) :fr: start_points(seg) + min_segment_time; % frames to record
%% Color in white outside the arena
imAlpha=ones(L-2*margin);

[X_movie,Y_movie] = meshgrid(-Lim:dr:Lim, -Lim:dr:Lim);

R_sq_movie=X_movie.^2 + Y_movie.^2; % squared distance from origin
[out_y_movie,out_x_movie]=find(R_sq_movie>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc_movie = sub2ind(size(X_movie), out_y_movie,out_x_movie); % linear indices of points outside of circled arena

imAlpha(idc_movie)=0; % making outisde values complete transperant

ticks = [1,25:25:125,151];
%% Decoding
leak_thresh=10^-10; % probabilites outside the arena smaller than this value won't be reflected to inside the arena and will be set to 0

for j = start_points(seg) - prior_gap :start_points(seg) + min_segment_time % running on all time steps relevant for the movie
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
               
        
%         ProCON(idc)=0; % can't be outside of circled arena
%         ProCON(ProCON<0)=0; % compensate for numerical error


%         S=S+sum(sum(ProCON));
%         disp(sum(sum(ProCON)));
%         Su=Su+sum(ProCON(idc));
        %
        ProSpikes=ones(L,L); % we first set the probability because of the spikes to ones
        ProSpikes(idc)=0; % can't be outside of the circled arena
        
        
        for k=used_neurons % running on all used neurons
            
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
        Norms(1,counter) = log(norm); % log normalization factor vs time
        probability(:,:,Current) = probability(:,:,Current)./norm; % normalizing likelihood

        %% MLE and its error (MSE)
        Pp=probability(:,:,Current);
        [~, max_idx]=max(Pp(:));
        [r,c]=ind2sub(size(Pp),max_idx);

        MLE(:,counter)=[c;r]; % [x;y] MLE [pixel]
        
        Xerror = c - (X_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 
        Yerror = r - (Y_traj(j)+Lim+1+margin); % adding Lim+1 and the margin to account for coordinates location 
        
        Pixeldist = sqrt((Xerror^2)+(Yerror^2)); % [Pixels]
        ERR(1,counter) = Pixeldist.*dr; % [cm]
        
        true_radius = sqrt((X_traj(j))^2 + (Y_traj(j))^2);
        estimated_radius = sqrt((c-Lim-margin)^2 + (r-Lim-margin)^2);
        R_distance(counter) = abs(true_radius - estimated_radius); % distance between true radius and estimated radius
        %% Movie
        if generate_movie==1 && sum(Movie_frames==j)==1 % if we record and it's a frame to record

            clf(fig,'reset'); % clearing the figure
            set(gcf, 'Position',  [500, 400, 500, 400]);
            hold on;

            tit = title(sprintf('Time = %.3f [s]',j*dt)); % title time in sec, in ms resultion
            tit.FontSize = 20;
            tit.FontWeight = 'normal'; % normal title font

            imagesc(probability(margin+1:end-margin,margin+1:end-margin,Current),'AlphaData',imAlpha); % posterior likelihood
            xlim([1,size(X_movie,1)]);
            ylim([1,size(X_movie,1)]);
            set(gca,'YDir','normal');
            set(gca,'linewidth',b_width);
            
            plot(c-margin,r-margin,'or','MarkerSize',15,'LineWidth',2.5); % plotting maximum likelihood estimate
            plot(X_traj(j)+Lim+1,Y_traj(j)+Lim+1,'xk','MarkerSize',15,'LineWidth',2); % real location

            xticks(ticks); yticks(ticks);
            xticklabels({'0','25','50','75','100','125','150'});
            yticklabels({'0','25','50','75','100','125','150'});

            ylabel('cm'); xlabel('cm');           
            set(gca,'fontsize',14); box on;
            axis square; % make x and y axis have same length (since the arena is a circle)

            Co = colorbar;
            Co.LineWidth=b_width;

            Mo1=getframe(gcf); % get movie frame
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

% Simulation_Results.Spike_train=ST; % saving the filtered spike train
%% Saving
if cluster==1
    cd(Na); % getting into Light/Dark directory
end
save(name,'Simulation_Results');
end