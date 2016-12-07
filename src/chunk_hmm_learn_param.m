function [P,Model,gamma] = chunk_hmm_learn_param(Data, Model,initialParam,varargin)
% function [rho, self_t, log_like, fm, T, rho_er, v, v_er, ...
%     initial_dist, mean_pause, mean_inchunk, ...
%     mean_pause_er, mean_inchunk_er, ...
%     chunks, cor_chunks, gamma] = ...
%     chunk_hmm_learn_param(mt_seq, er_seq, varargin)
% Finds best parameters for chunks using Baum-Welch algorithm
% INPUT:
%  Data: Structure of the data to be fitted
%    mt_seq:        n x k matrix of reaction times (n: number of trials, k: sequence length)
%    er_seq:        n x k matrix of error data (optional)
%    seqType:       n x 1 vector of sequence type if multiple different
%                   sequences are fitted simulatenously  (optional)
%    trialNumber    n x 1 trial number per sequence   (optional)
%  Model: Structure of model options (see Model Structure)
%  initialParam: Structure of intial parameters (see Model Parameters)
%  VARAGIN:
%    'verbose',0/1  printing fit?
% OUTPUT:
%    Params:         Structure of model parameters
%    Model:          Structure of model options
%    gamma:          NxQ Posterior probability on different ways of chunking
% MODEL STRUCTURE:
%    chunks:         Types of chuning structures that are considered
%    fit_<paramatername>: Fit or not fit the parameter to the data?
%
% MODEL PARAMETERS:
%    initial_dist:    initial distribution of chunks (one per seuqence type)
%    self_t:          self-transition parameter
%    T:               transition matrix
%    mean_pause:      mean reaction time pause at beginning of chunk
%    mean_inchunk:    mean reaction time within a chunk
%    v:               variance in reaction time
%    rho:             correlation parameter for reaction time
%    mean_pause_er:   Error frequency in the beginning of a chunk
%    mean_inchunk_er: Error frequency in the beginning of a chunk
%    v_er:            variance in errors
%    rho_er :         correlation parmater for errors within a chunk
%
% OUTPUT parameters:
%    log_like: log-likelihood on training data
%    fm: forward messages for hmm learning
%    chunks: QxK matrix indicating chunk indetities
%    cor_chunks:
% some parameters can be given
p = inputParser;
p.KeepUnmatched = true;
addParamValue(p, 'verbose', false, @islogical);
addParamValue(p, 'maxIter', 30, @isnumeric);
parse(p, varargin{:});
verbose = p.Results.verbose;
maxIter = p.Results.maxIter; 
tic;

P = initialParam;

% Determine if we are fitting multiple sequence types 
% at the same time 
if (isfield(Data,'seqType')) 
    seq = unique(Data.seqType); 
    numSeq = length(seq); 
else 
    Data.seqType = ones(size(Data.mt_seq,1),1); 
    numSeq = 1; 
end; 


if (Model.fit_rt && ~isfield(P,'rho')) 
    P.rho = 0;
end

if (Model.fit_er && ~isfield(P,'rho_er')) 
    P.rho_er = 0;
end

log_like = nan;

% Make chunking structure if missing
if ~isfield(Model,'chunks') || isempty(Model.chunks)
    Model.chunks = create_chunks_nospace('n_seqlen', size(Data.mt_seq, 2));
end

% Generate raw correlation matrix for each chunk structure
Model.cor_chunks = to_corr_chunks(Model.chunks);

% Number of chunking structures
n_chunks = size(Model.chunks, 1);
nTrials  = size(Data.mt_seq,1); 
% ---------------------------------------
% Set the intial values for all parameters
% Self transition parameters

if ~isfield(P,'self_t')
    P.self_t = 0.6;
end;

% Set up initial guess of transition matrix
P.T = ((ones(n_chunks, n_chunks)-eye(n_chunks))*...
    (1-P.self_t))/(n_chunks-1) + ...
    eye(n_chunks)*P.self_t;

% Initial guesses of variance parameters: Exclude missing observations
if (Model.fit_rt) 
    if (~isfield(P,'v'))
        P.v    = var(Data.mt_seq(~isnan(Data.mt_seq)));
    end; 
    if (~isfield(P,'mean_pause'))
        P.mean_pause = prctile(Data.mt_seq(~isnan(Data.mt_seq)),75);
        P.mean_inchunk = prctile(Data.mt_seq(~isnan(Data.mt_seq)),25);
    end; 
end;

if (Model.fit_er) 
    if (~isfield(P,'v_er'))
        P.v_er = var(Data.er_seq(~isnan(Data.er_seq)));
    end; 
    if (~isfield(P,'mean_pause_er'))
        P.mean_pause_er = prctile(Data.mt_seq(~isnan(Data.mt_seq)),75);
        P.mean_inchunk_er = prctile(Data.mt_seq(~isnan(Data.mt_seq)),25);
    end; 

end;

% initial guess for the intial distribution: uniform
if (~isfield('initial_dist',P))
    P.initial_dist = ones(numSeq, n_chunks)/n_chunks;
end;


% Initialize
n_iteration = 0;
not_exit = true;

% alpha = 1;
% beta = 1/(n_chunks-1);
while not_exit
    current_log_like = sum(log_like);
    % COMPUTATION OF EXPECTATION
    p_obs_un = create_emission_for_chunks(P,Data,Model);
    
    gamma  = nan(nTrials,n_chunks); 
     marg_epsilon = zeros(n_chunks,n_chunks); 
    for i=1:numSeq 
        indx = find (Data.seqType==i); 
        [~, ~, g, log_like(indx), me] = hmm_inference(p_obs_un(indx,:), P.T, 'initial_dist', P.initial_dist(i,:));
        % Collect sufficient statistics 
        numTrials(i)=length(indx);
        numStay(i)  = trace(g(1:end-1, :)' * g(2:end, :)); 
        gammaInit(i,:)   = g(1,:); 
        marg_epsilon = marg_epsilon+numTrials(i)*me; 
        gamma(indx,:)    = g(2:end,:); % Posterior probabilities 
    end; 

    n_iteration = n_iteration + 1;
    % Feedback if required 
    if verbose
        
        fprintf('Iteration: %.0f (elapsed time: %.1f s)\n', n_iteration, ...
            toc);
        fprintf(['\t log_like        = %f\n'], sum(log_like));
        fprintf(['\t self_t        = %f\n'], P.self_t);
        fprintf(['\t mean_pause    = %f\n'], P.mean_pause);
        fprintf(['\t mean_inchunk  = %f\n'], P.mean_inchunk);
        fprintf(['\t rho           = %f\n'], P.rho);
        
    end
    
    % Update likelihood
    delta_fval = (sum(log_like) - current_log_like);
    
    if delta_fval<0 
        fprintf('likelihood decreased... exiting\n'); 
        P     = Pbest;             % Parameters
        gamma = gamma_best;        % Expectation
        log_like = log_like_best;  % Loglikelihood
        not_exit = false;
    elseif delta_fval<10e-3
        fprintf('likelihood changed less than 10e-3... exiting\n'); 
        not_exit = false;
    elseif n_iteration > maxIter 
        fprintf('maximal number of iterations reached... exiting\n'); 
        not_exit = false;
    else
        % save previous results
        Pbest         = P;
        gamma_best    = gamma;
        log_like_best = log_like;
        
        % MAXIMIZATION OF EXPECTATION
        % Restimate transition maximization
        if Model.fit_T
            T = marg_epsilon/sum(numTrials);
            
            P.self_t = sum(numStay)/sum(numTrials); 
        end
        if Model.fit_T && Model.diagonal_T
            % Transition probability probability of self-transition
            T = ((ones(n_chunks, n_chunks)-eye(n_chunks))*...
                (1-P.self_t))/(n_chunks-1) + ...
                eye(n_chunks)*P.self_t;
        end
        
        P.initial_dist = gammaInit; 
        
        % reestimate MT parameters
        if (isfield(Data,'mt_seq'))
            % find mean_pause
            switch (Model.fit_rt)
                case 1          % Learn constant Pause model
                    [P.mean_pause, P.mean_inchunk] = learn_pause(Model.chunks, Data.mt_seq, gamma);
                case 2
                    % Learn detrneding pause model
            end;
            
            P.v = learn_variance(Model.chunks, Data.mt_seq, gamma, P.mean_pause, P.mean_inchunk,Model.fit_rt);
            
            % reestimate covariance for movement time
            if Model.fit_rho
                P.rho = learn_cor(Model.chunks, Model.cor_chunks, Data.mt_seq, gamma, ...
                    P.mean_pause, P.mean_inchunk, P.v, Model.fit_rt);
            end
        end
        
        if (isfield(Data,'er_seq'))
            % find mean_pause
            if Model.fit_er
                [P.mean_pause_er, P.mean_inchunk_er] = learn_pause(chunks, er_seq, gamma);
            end
            
            P.v_er = learn_variance(chunks, er_seq, gamma, ...
                mean_pause_er, mean_inchunk_er, fit_er);
            
            if P.fit_rho_er
                P.rho_er = learn_cor(chunks, cor_chunks, er_seq, gamma, ...
                    mean_pause_er, mean_inchunk_er, v_er, fit_er);
            end
        end
    end
end

