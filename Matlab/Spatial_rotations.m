function [TC] = Spatial_rotations(m1_id,m2_id,m3_id,Tuning_Curves,max_rot,tr,identical)
% Inputs: allocation of cells to modules and their tuning curves
%         max_rot is the maximal spatial rotation to introduce to the modules [degree]
%         identical specifies if shifts are identical to all cells
% Output: TC - rotated tunning curves
% This works because the tuning curve is centered around the origin
rng(tr);

rot = round(2.*max_rot.*rand(1,3) - max_rot); % rotation [degree]

if identical==1 % if we actually perform identical rotations
    rng(tr*1000); % new seed for identical shifts

    rot = round(2.*max_rot.*rand(1,3) - max_rot); % rotation [degree]
    rot(2:3) = rot(1); % same rotation for all [degree]
end

%% Rotating rate maps
TC = zeros(size(Tuning_Curves)); % will contain rotated rate maps

for k=1:1:numel(m1_id) % running on first module neurons
    TC(:,:,m1_id(k)) = imrotate(Tuning_Curves(:,:,m1_id(k)),rot(1),'bilinear','crop');
end
for k=1:1:numel(m2_id) % running on first module neurons
    TC(:,:,m2_id(k)) = imrotate(Tuning_Curves(:,:,m2_id(k)),rot(2),'bilinear','crop');
end
for k=1:1:numel(m3_id) % running on first module neurons
    TC(:,:,m3_id(k)) = imrotate(Tuning_Curves(:,:,m3_id(k)),rot(3),'bilinear','crop');
end

end