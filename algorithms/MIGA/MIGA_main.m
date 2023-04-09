function [ bestind, bestfit, nite, lastpop, lastfit, history ] = MIGA_main(lb,ub,dim,fobj)
% Define heuristic function options (optional)
opts.ninfo = 0; % Verbosity level (0=none, 1=minimal, 2=extended)
opts.label = 10; % Label (identification purposes)
opts.dopar = 0; % Parallel execution of fitness function
opts.nhist = 2; % Save history (0=none, 1=fitness, 2=all{pop,fit})

% Define Islands Model parameters
ni = 3; % Number of islands
ngg = 6; % Number of global iterations
nemi = 5; % Number of emigrants
goal = -Inf; % Target fitness value
heufun = @aga; % Heuristic function to use with the Islands Model

% Define AGA parameters
ng = 6; % Number of local generations
np = 40; % population of each island
N = [3,... % Number of elites
    floor(np*0.1),... % Number of mutants
    floor(np*0.05),...% Number of newcomers
    floor(np*0.2)]; % Number of parents

% Auxiliary function
ranrange = @(a,b,n) a + (b-a)*rand(n,1); % n random values between a i b

% Define AGA functions
unifun = @(x,f) deal(x,f); % Discard identical individuals: currently not in use
fitfun = @(x) fobj(x); % Fitness function (to be minimized)
mutfun = @(x,f) x + ranrange(-0.1,0.1,dim); % Mutation: small random mov
repfun = @(x,y,fx,fy) (x+y)/2; % Reproduction: average
ranfun = @() ranrange(lb,ub,dim); % Random individual
prifun = @(x) fprintf('%f',x(dim)); % Print an individual

% Assemble AGA data structure
DATA{1} = ng;
DATA{2} = N;
DATA{3} = unifun;
DATA{4} = fitfun;
DATA{5} = mutfun;
DATA{6} = repfun;
DATA{7} = ranfun;
DATA{8} = []; % Set as empty so that the program will not print information each iteration

% Randomize random seed
% rng('shuffle'); % We don't want repeatability in the heuristic
rng(0); % We don't want repeatability in the heuristic

% We can just give the number of islands and individuals (pops = [ni np]),
% then aim generates the populations calling our ranfun. Or we can
% generate our own initial populations:
pops = cell(1,ni);
for isl=1:ni
    for j=1:np
        pops{isl}{j} = ranfun(); % Create random individual
    end
end

% Execute Islands Model + Genetic Algorithm (GA)
[ bestind, bestfit, nite, lastpop, lastfit, history ] = aim ( ...
    opts, pops, ngg, nemi, goal, heufun, DATA );

end