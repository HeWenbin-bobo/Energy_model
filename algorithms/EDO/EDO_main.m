function [ bestind, bestfit ] = EDO_main(lb,ub,dim,fobj)
% Define heuristic function options (optional)
NP = 5;
Max_iter = 5;

% Randomize random seed
% rng('shuffle'); % We don't want repeatability in the heuristic
rng(0); % We don't want repeatability in the heuristic

% Execute EDO
[bestfit, bestind]=EDO(NP,Max_iter,lb,ub,dim,fobj);

end