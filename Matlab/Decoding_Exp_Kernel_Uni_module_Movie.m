function [ ] = Decoding_Exp_Kernel_Uni_module_Movie( Rat,L_D,Tau ) %
% Decoding Tergeir data set (Light/Dark) using estimated tuning curves and
% recorded spikes
% Use Biological estimated tuning curves and Biological spikes to decode

%% Video properties
seg=10; % the segment to simulate and make movie of
if L_D==0 % dark
    Na = 'Dark';
elseif L_D==1 % light
    Na = 'Light';
end
vid_name = [Na,sprintf(' Log Likelihood %d.mp4',seg)];
b_width = .75; % width of plot box axes
%% Decoding properties
ori_dir=pwd;
speed_threshold=3; % diluting spikes trains  - but according to speed thresohld
Modules = 3; % plotting movie only for cell with 3 modules
%%
cluster=0; % run on cluster? set to 1

Rel_Tau=10; % how many Tau's back is it relevant to look at on the spike trains? (should be the same as in multi module code)
% Beyond that spikes are weighted by a 0 instead of the kernal
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
margin=0; % [cm] margin addition to the radius for not losing (and reflecting) probability leaks from outside the real arena back in
% Don't need margin in this case since there is no leak of probability (Must be as in multi exp deocded)

L = numel(-Lim-margin:dr:Lim+margin);  % size of arena+margin [pixels] = [cm] resulotion

[X,Y] = meshgrid(-Lim-margin:dr:Lim+margin, -Lim-margin:dr:Lim+margin);

R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>=Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

[in_y,in_x]=find(R_sq<Lim^2); % coordinates inside the circled arena
idc_in = sub2ind(size(X), in_y,in_x); % linear indices of points inside of circled arena

dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]
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

if L_D==1 % loading Light spike trains
    load('ST_L.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory
elseif L_D==0 % loading Dark spike trains
    load('ST_D.mat'); % loading emitted spikes of all neurons binned in time corresponding to trajectory   
end

cd(ori_dir); % go back to original directory where we came from
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
%% Uni - module neurons
used_neurons_modules=cell(1,4); % will contain the neurons from each module
used_neurons_modules{1} = [m1_id',m2_id',m3_id']; % all the neurons for multi module decoding as reference
used_neurons_modules{2} = m1_id'; % neurons that belong to the 1st module
used_neurons_modules{3} = m2_id'; % neurons that belong to the 2nd module
used_neurons_modules{4} = m3_id'; % neurons that belong to the 3rd module

% SAME VALUES AS IN USED IN UNI MODULE DECODING
if Rat==1 || Rat==11 % Bubble and Bubble 2nd dataset
    R_Real_Loc = [21,28,40];
elseif Rat==2 % Roger
    R_Real_Loc = [25,40,52];
elseif Rat==4 % Rolf
    R_Real_Loc = [22,28,43];
end

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
%% Color in white outside the arena
imAlpha=ones(L);
imAlpha(idc)=0; % making outisde values complete transperant

ticks = [1,25:25:125,151];
%% Decoding likelihood from the spikes
clf(fig,'reset'); % clearing the figure
set(gcf, 'Position',  [200, 500, 1000, 700]);
image(100); % initial frame to start on
       
for j=start_points(seg) :fr: start_points(seg) + min_segment_time % running on segment

   Log_P=zeros(L,L,Modules+1); % will contain the sum of the log of the probability for all modules and for multi module

   Tau_Spikes = ST(:,j-Rel_times:j); % spike trains of all neurons during relevant backward times
   spike_count = sum(Tau_Spikes,2); % sum of spikes during this time interval


   for mo=1:Modules+1 % running on the modules + multi module

       mo_neurons=used_neurons_modules{mo};
       fired = find(spike_count(mo_neurons)>0); % find which neurons from the mo module fired in this time interval

       for fi=1:numel(fired) % running only on neurons that did fire during this time interval

           rf = Tuning_Curves(:,:,mo_neurons(fired(fi))); % receptive field of the relevant neuron that did fire
           st=Tau_Spikes(mo_neurons(fired(fi)),:); % spike train during this time interval (of the relevant neuron that did fire)

           spike_timings = find(st>0); % timing of spikes
           num_of_spikes = st(spike_timings); % how many spikes where emitted in these times?
           weighted_spikes = num_of_spikes .* exp(- ( (Rel_times+1) - spike_timings) / Tau_indices); % weighing spikes by their timing and the # of emitted spikes
           weighted_spikes = sum(weighted_spikes); % total effect of the spikes of that neuron

           Log_P(:,:,mo) = Log_P(:,:,mo) + weighted_spikes.*log(rf); % log likelihood of the mo'th module
       end
   end
   %% Normalize the probabilites (for colorbar)
%        P = exp(Log_P); % convert to probability
% 
%        for mo = 1:Modules+1 % running on the modules + multi module
%            p = P(:,:,mo); % probability for the mo'th module
% 
%            l_p = p(idc_in); % points inside the arena
% 
%            norm = sum(l_p) .* (dr.^2);
%            l_p = l_p./norm; % normalized probabilty inside the arena
%            
%            p(idc_in) = l_p; % put back the inside arena points in the whole arena matrix
%            p(idc) = min(l_p); % set outside point with minimal value of inside point for imagesc color scaling
%        
%            Log_P(:,:,mo) = log(p); % normalized log likelihood
%        end
   %% Multi module log likelihood
   Log_P_Multi=Log_P(:,:,1); % Multi-module log likelihood
   [~, max_idx]=max(Log_P_Multi(:));
   [r,c]=ind2sub(size(Log_P_Multi),max_idx);

   %% Extract uni module positions
   XX = X - (c-Lim-margin); % shifted x coordinates
   YY = Y - (r-Lim-margin); % shifted y coordinates
   radi = sqrt( XX.^2 + YY.^2); % corresponding radius

   r_locs = zeros(Modules,1); % will contain the y local positions for each module within its unit cell
   c_locs = zeros(Modules,1); % will contain the x local positions for each module within its unit cell
   for um=1:3 % run on the three uni-modules
   
       r_real_loc = R_Real_Loc(um); % radi of the um module
       [local_out_y,local_out_x]=find(radi>=r_real_loc); % coordinates to set to 0 firing rate since they are outside of circle
       local_idc = sub2ind(size(X), local_out_y,local_out_x);
       
       log_p = Log_P(:,:,um+1); % log likelihood of the um'th module
       log_p(local_idc)=min(min(log_p)); % leaving just a circle around desired location with higest log(probability) for MLE and ERR readout

       [~, max_idx]=max(log_p(:));
       [r_locs(um),c_locs(um)]=ind2sub(size(log_p),max_idx);
   end
   %% Movie
   if generate_movie==1
       
       clf(fig,'reset'); % clearing the figure
       set(gcf, 'Position',  [400, 300, 700, 600]);
       sgt = sgtitle(sprintf('Time = %.3f [s]',j*dt)); % title time in sec, in ms resultion
       sgt.FontSize = 20;

       for gm=1:Modules+1 % runnning on all modules (+Multi module)

           subplot(2,2,gm);

           l_P = Log_P(:,:,gm); % log likelihood of the gm group
           l_P(idc) = min(min(l_P(idc_in))); % set outside point with minimal value of inside point for imagesc color scaling

           % Normalizing the log likelihood so colorbar won't move as much
           l_P = l_P - max(max(l_P)); % setting max value of log likelihood to 0


           imagesc(l_P,'AlphaData',imAlpha); % posterior likelihood
           set(gca,'YDir','normal');
           set(gca,'linewidth',b_width);
           Co = colorbar;
           Co.LineWidth=b_width;

           % Specifing colorbar ticks
           d_L = abs(min(min(l_P))); % the maximal difference between log(L) values as the maximum was set to 0
           Co_Ticks = -100:10:0;
           if d_L>100
               Co_Ticks = -300:50:0;
               if d_L>300
                   Co_Ticks = -1000:100:0;
               end
           end
           Co.Ticks=Co_Ticks;

           hold on;
           plot(c,r,'or','MarkerSize',15,'LineWidth',2.5); % plotting maximum likelihood estimate
           plot(X_traj(j)+Lim+margin+1,Y_traj(j)+Lim+margin+1,'xk','MarkerSize',15,'LineWidth',2); % real location

           xticks(ticks); yticks(ticks);
           xticklabels({'0','25','50','75','100','125','150'});
           yticklabels({'0','25','50','75','100','125','150'});

           ylabel('cm'); xlabel('cm');
           title_name=sprintf('Module %d',gm-1); % module # title
           if gm==1 % if it's multi module
              title_name = 'Multi module'; % mutli module title
           end
           ti = title(title_name); % title
           ti.FontWeight = 'normal'; % normal title font
           
           set(gca,'fontsize',14);     
           axis square; % make x and y axis have same length (since the arena is a circle)

           %% Plot uni-module unit cell circles
           if gm>1 % if it's uni module decoding
               Ci = viscircles([c,r],R_Real_Loc(gm-1)); % plot unit circle for uni modules
               
               Ci.Children(2).LineStyle = 'none'; % remove background line

               Ci.Children(1).LineStyle = '--'; % define circle line style
               Ci.Children(1).Color = [.4,.4,.4]; % define circle color
               Ci.Children(1).LineWidth = 1; % define line width

               Ci_X = Ci.Children(1).XData; % circle x points
               Ci_Y = Ci.Children(1).YData; % circle y points

               Ci_X(isnan(Ci_X))=[]; % remove nan values x
               Ci_Y(isnan(Ci_Y))=[]; % remove nan values x

               Ci_Xx = Ci_X - Lim; % center around the origin
               Ci_Yy = Ci_Y - Lim; % center around the origin

               Ci_r = sqrt( Ci_Xx.^2 + Ci_Yy.^2);
               Ci_out = find(Ci_r>=Lim); % find points outside the arena
              
               xi_o = Ci_X(Ci_out); % the outside x points of the circle
               yi_o = Ci_Y(Ci_out); % the outside y points of the circle
               
               % Might need to circshift points so superimposed line won't cross the arena
               dxx = xi_o(2:end)-xi_o(1:end-1); % diff x between consecutive points
               dyy = yi_o(2:end)-yi_o(1:end-1); % diff y between consecutive points

               Di = sqrt( dxx.^2  + dyy.^2); % distnace between line's consecutive points
               Di = roundn(Di,-4); % rounding decimal before comparing

               U = unique(Di); % unique values of the distances
               if numel(U)>1 % if there is more than a single value than the circle is crossing the arena we need the top one
                   [~,l]=max(Di); % point to circhshift

                   xi_o = circshift(xi_o,-l); % circshifting x coordinate
                   yi_o = circshift(yi_o,-l); % circshifting y coordinate
               end

               plot(xi_o,yi_o,'w','LineWidth',Ci.Children(1).LineWidth); % superimposed outside lines with white

               %% Plot uni module estimate
               plot(c_locs(gm-1),r_locs(gm-1),'p','Color',[0,0.6,0],'MarkerSize',16,'LineWidth',1.5);
           end

       end
       Mo1=getframe(gcf); % get movie frame
       writeVideo(vidfile,Mo1);
   end
end
  
end