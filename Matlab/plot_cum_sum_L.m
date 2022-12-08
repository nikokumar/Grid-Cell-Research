clear;

ori_dir=pwd;
dt=1/120;
Rat=3; % 1 for Bubble, 11 for Bubble 2, 2 gor Roger, 3 for Coconut and 4 fo Rolf

dir = '/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Dark and Light TC no dilution'; % directory for non adjusted rate results
cd(dir);

if Rat==1
    cd('Bubble');
elseif Rat==11
    cd('Bubble 2');
elseif Rat==2
    cd('Roger');
elseif Rat==3
    cd('Coconut');
elseif Rat==4
    cd('Rolf');
end

cd('Markov results');

load('Dark_spikes,Dark_TC=0'); % loading dark results
Nd = Simulation_Results.Norms; % log(norms) vs time from dark

load('Light_spikes,Dark_TC=0'); % loading light results
Nl = Simulation_Results.Norms; % log(norms) vs time from dark

td=0:dt:(numel(Nd)-1)*dt; % dark time
tl=0:dt:(numel(Nl)-1)*dt; % light time
td=td./60; % min units
tl=tl./60; % min units

Cd = cumsum(Nd); % cummulative sum dark
Cl = cumsum(Nl); % cummulative sum light
%% plot
close all;
subplot(1,2,1); % dark
plot(td,cumsum(Nd));
xlim([0,td(end)]);
xlabel('Time [min]');
ylabel('\Sigma_{\it{i}=0}^{\it{t}}log(\it{Z}_{\it{i}})');
title('Dark');
set(gca,'fontsize',16); grid on;

subplot(1,2,2); % light
plot(tl,cumsum(Nl));
xlim([0,tl(end)]);
xlabel('Time [min]');
ylabel('\Sigma_{\it{i}=0}^{\it{t}}log(\it{Z}_{\it{i}})');
title('Light');
set(gca,'fontsize',16); grid on;

set(gcf, 'Position',  [100, 300, 800, 250]);

cd(ori_dir);