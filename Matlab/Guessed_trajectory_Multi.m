function [X_guessed, Y_guessed] = Guessed_trajectory_Multi(Tt,Lim,rounding,tr)
% Provides guessed trajectory when using uniform posterior: gives higher
% bound of error
%% Generate random points within the circle
rng(tr);

X_guessed  = (Lim-(-Lim)).*rand(Tt,1) + (-Lim);
Y_guessed  = (Lim-(-Lim)).*rand(Tt,1) + (-Lim);

R = sqrt(X_guessed.^2 + Y_guessed.^2); % the radii of randomly distibuted point

over=find(R>Lim); % find points outside the arena
n_over=numel(over); % how many points?
while n_over>0 % while we have points outside the circle
    
    X_guessed(over) = (Lim-(-Lim)).*rand(n_over,1) + (-Lim); % new x points
    Y_guessed(over) = (Lim-(-Lim)).*rand(n_over,1) + (-Lim); % new y points
    
    
    R(over) = sqrt(X_guessed(over).^2 + Y_guessed(over).^2); % new radii
    over=find(R>Lim); % find points outside the arena
    n_over=numel(over); % how many points?
end

X_guessed = roundn(X_guessed,rounding); % same rounding as x actual tracjectory
Y_guessed = roundn(Y_guessed,rounding); % same rounding as y actual tracjectory


end