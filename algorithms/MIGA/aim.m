function [ bestind, bestfit, nite, lastpop, lastfit, history ] = aim ( ...
    opts, pops, ngg, nemi, goal, heufun, DATA )
%AIM finds minimum of a function using the Islands Model
%
%Programmers:   David de la Torre   (UPC/ETSEIAT)
%               Manel Soria         (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          10/05/2018
%Revision:      2
%
%Usage:         [ bestind, bestfit, nite, lastpop, lastfit, history ] = ...
%                   AIM( opts, pops, ngg, nemi, goal, fheu, DATA )
%
%Inputs:
%   opts:       function control parameters [struct] (optional)
%       ninfo:  verbosity level (0=none, 1=minimal, 2=extended)
%       label:  integer number that precedes the prints in case output is
%               to be filtered
%       dopar:  parallel execution of fitness function [1,0]
%       nhist:  save history (0=none, 1=fitness, 2=all{pop,fit})
%                   0: history = []
%                   1: history(ngg,ni) = bestfit(i)
%                   2: history{ngg,ni} = (hist, nite)
%                      history{ngg,ni+1} = best individual of each gg
%                      history{ngg,ni+2} = best individual fitness
%   pops:       list with initial population elements
%   ngg:        number of global generations
%   nemi:       number of emigrations between islands per global generation
%   goal:       if function value is below goal, iterations are stopped
%   heufun:     handle to the heuristic function (@aga, @ade, etc.)
%   DATA:       cell array with the specific parameters and callback
%               functions of the heuristic function.
%               The heuristic is called in the following way:
%               [outputs] = fheu(opts,pop,goal,DATA);
%               Do ensure the DATA array has all the variables in the order
%               required by the heuristic function.
%               Note that last variable in DATA should be prifun(), and
%               last-1 variable in DATA should be ranfun().
%
%Outputs:
%   bestind:    best individual (among all the islands)
%   bestfit:    fitness value of best individual
%   nite:       number of global iterations (generations) performed
%   lastpop:    list with last populations of each island
%   lastfit:    best fitness values of last population of each island
%   history:    array with saved global history array

% Get configuration options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else, ninfo = 0; end
if isfield(opts,'label'), label = opts.label; else, label = 0; end
if isfield(opts,'dopar'), dopar = opts.dopar; else, dopar = 0; end
if isfield(opts,'nhist'), nhist = opts.nhist; else, nhist = 0; end

% Create history array
history = [];

% Build population
if (~iscell(pops)) % Only population size is given                                   
    ni = pops(1); % Number of islands
    np = pops(2); % Population of an island
    pops = cell(1,ni); % Preallocate var
    for island=1:ni % Fill islands with population
        for i=1:np % Fill population with individuals
            pops{island}{i} = DATA{end-1}(); % Create random individual
        end
    end
else % Initial population is given. Check structural consistency
    ss = size(pops); % Size (m,n) of input population
    if (ss(1)~=1) % Error in pops shape
        error('AIM global population shape must be: pops{1,ni}');
    end
    ni = ss(2); % Number of islands
    for i=1:ni
        ss = size(pops{i}); % Size
        if (ss(1)~=1)
            error('AIM population shape must be: pop{1,ni}');
        end
        if i==1, np = ss(2); % Population length
        else, if np~=ss(2), error('AIM population length mismatch'); end
        end
    end
end

% Show info
if ninfo>0
    fprintf('AIM begin ngg=%d ni=%d np=%2d\n',ngg,ni,np);
end

% Preallocate top best fitness
topbestfit = 0;

% Iterate through generations
for gg=1:ngg
    
    % Bests of each island
    ibestinds = cell(ni,1);
    ibestfits = zeros(ni,1);
    
    % Evolve each island separately
    for island=1:ni
        
        % Island options
        iopts.ninfo = ninfo-1;
        iopts.label = label + gg*1000 + island;
        iopts.dopar = dopar;
        iopts.nhist = nhist;
        
        % Execute heuristic algorithm
        [ibestind, ibestfit, initer, ilastpop, ~, hist] = heufun ( ...
            iopts, pops{island}, goal, DATA{:} );
        
        % Save values at the end of local island iteration
        pops{island} = ilastpop; % Save last population
        ibestinds{island} = ibestind; % Save best individual
        ibestfits(island) = ibestfit; % Save best fitness value

        % Save history
        if nhist>1 % Save full history for each island
            history{gg,island} = {hist,initer}; %#ok
        elseif nhist>0 % Save best fitness only
            history(gg,island) = ibestfit; %#ok
        end

        % Show extended info
        if ninfo>1
            fprintf('AIM label=%d gg=%d island=%d fitbest=%e',...
                label,gg,island,ibestfit);
            if ~isempty(DATA{end}) % Print best individual
                fprintf(' best='); DATA{end}(ibestind);
            end
            fprintf('\n');
        end
        
    end
    
    % Find best individual of all islands and island where it lives
    [fbest, ibest] = min(ibestfits);
    
    % Save top best individual and its fitness
    if fbest<topbestfit || gg==1 % Better individual found, or first gg
        topbestind = ibestinds{ibest};
        topbestfit = fbest;
    end
    
    % Save history
    if nhist>1 % Save full history for each island
        history{gg,ni+1} = ibestinds{ibest}; %#ok % Best individual
        history{gg,ni+2} = fbest; %#ok % Fitness of best individual
    end

    % Check if reached target fitness or max generations 
    if fbest<=goal || gg>=ngg % Target achieved; end simulation
        
        % Save output values
        bestind = topbestind;
        bestfit = topbestfit;
        nite = gg;
        lastpop = pops;
        lastfit = ibestfits;
        
        % Show info if required
        if ninfo>0
            fprintf('AIM label=%d nite=%2d fitbest=%8.3e',label,gg,bestfit);
            if ~isempty(DATA{end}) % Print best individual
                fprintf(' best='); DATA{end}(bestind);
            end
            if gg==ngg % Maximum number of iterations reached
                fprintf(' max. iterations reached, leaving\n');
            else % Goal has been reached
                fprintf(' goal=%e achieved, leaving\n',goal);
            end
        end
        
        % Stop iterating
        break;
        
    end
    
    % Show info
    if ninfo>0
        fprintf('AIM label=%d gg=%d fitbest=%8.3e islandbest=%d',...
            label,gg,fbest,ibest);
        if ~isempty(DATA{end})
            fprintf(' best='); DATA{end}(ibestinds{ibest});
        end
        fprintf('\n');
    end
    
    % Emigration
    % Copy nemi individuals (among the best) on origin island, and then
    % replace nemi individuals (among the worst) on destination island
    for island=1:ni
        
        % Select origin and destination islands (ring topology)
        dest = island; % Destination island
        orig = island+1; if orig>ni, orig=1; end % Origin island
        
        % Migrate individuals
        for migrant=1:nemi
            
            % Individual that will be killed by the migrant
            killed = length(pops{island}) - migrant + 1;
            
            % Migrate individual
            pops{dest}{killed} = pops{orig}{migrant};
            
            % Extended info
            if ninfo>1
                fprintf(['AIM label=%d gg=%d migrant from ',...
                    'island=%d ind=%2d replaces island=%d ind=%2d\n'],...
                    label,gg,orig,migrant,dest,killed);
            end
            
        end
        
    end

end

end

