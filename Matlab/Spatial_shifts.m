function [TC] = Spatial_shifts(m1_id,m2_id,m3_id,Tuning_Curves,max_shift,tr,identical)
% Inputs: allocation of cells to modules and their tuning curves
%         max_shift is the maximal spatial shifts to introduce to the modules [cm]=[pixel]
%         identical specifies if shifts are identical to all cells
% Output: TC - shifted tunning curves
rng(tr);

shifts_x = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in x axis
shifts_y = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in y axis

if identical==1 % if we actually perform identical spatial shifts
    rng(tr*1000); % new seed for identical shifts
    shifts_x = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in x axis
    shifts_y = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in y axis

    shifts_x(2:3) = shifts_x(1); % same x shifts for all
    shifts_y(2:3) = shifts_y(1); % same y shifts for all
end


L=size(Tuning_Curves,1); % size of length of a tuning curve (including margin for leak)
LL=L+2*max_shift; % size of extended tuning curve
N=size(Tuning_Curves,3); % # of neurons

TC=zeros(LL,LL,N); % embedding inside a greater arena
TC(max_shift+1:max_shift+L,max_shift+1:max_shift+L,:)=Tuning_Curves; % embedding


for k=1:1:numel(m1_id) % running on first module neurons
    TC(:,:,m1_id(k)) = circshift(TC(:,:,m1_id(k)),[shifts_y(1),shifts_x(1)]);
end
for k=1:1:numel(m2_id) % running on first module neurons
    TC(:,:,m2_id(k)) = circshift(TC(:,:,m2_id(k)),[shifts_y(2),shifts_x(2)]);
end
for k=1:1:numel(m3_id) % running on first module neurons
    TC(:,:,m3_id(k)) = circshift(TC(:,:,m3_id(k)),[shifts_y(3),shifts_x(3)]);
end

TC = TC(max_shift+1:max_shift+L,max_shift+1:max_shift+L,:);
end