function [P,Model,gamma] = chunk_hmm_learn_param(Data, Model,initialParam,varargin)
% function [rho, self_t, log_like, fm, T, rho_er, v, v_er, ...
%     initial_dist, mean_pause, mean_inchunk, ...
%     mean_pause_er, mean_inchunk_er, ...
%     chunks, cor_chunks, gamma] = ...
%     chunk_hmm_learn_param(mt_seq, er_seq, varargin)
% Finds best parameters for chunks using Baum-Welch algorithm
% INPUT:
%  Data: Structure of the data to be fitted
%    mt_seq:       n x k matrix of reaction times (n: number of trials, k: sequence length)
%    er_seq:       n x k matrix of error data (optional)
%    seqType:      n x 1 vector of sequence type if multiple different
%                   sequences are fitted simulatenously  (optional)
%    trialNum      n x 1 trial number per sequence   (optional)
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
Model.cor_chunks = to_corr_chunks(Model.chunks,Model.cov_structure);

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
switch (Model.fit_rt)
    case 1 % This is fitting the means
        if (~isfield(P,'mean_pause'))
            P.mean_pause = prctile(Data.mt_seq(~isnan(Data.mt_seq)),75);
            P.mean_inchunk = prctile(Data.mt_seq(~isnan(Data.mt_seq)),25);
        end;
    case 2 % This fit between and within with seperate exponentials
        start = prctile(Data.trialNum,3); 
        ende = prctile(Data.trialNum,95); 
        P.pause_tau =  ende/3; 
        P.inchunk_tau= ende/3; 
        dataStart = Data.mt_seq(Data.trialNum<=start,:); 
        dataEnd   = Data.mt_seq(Data.trialNum>=ende,:); 
        P.pause_B   = prctile(dataStart(~isnan(dataStart)),75);
        P.inchunk_B = prctile(dataStart(~isnan(dataStart)),25);
        P.pause_A   = prctile(dataEnd(~isnan(dataEnd)),75);
        P.inchunk_A = prctile(dataEnd(~isnan(dataEnd)),25);
end;

if (~isfield(P,'v'))
    P.v    = var(Data.mt_seq(~isnan(Data.mt_seq)));
    if (Model.fit_rt_var==2) 
        P.v(1,2)=P.v(1,1)/2; 
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
        switch(Model.fit_rt)
            case 1
                fprintf(['\t mean_pause    = %f\n'], P.mean_pause);
                fprintf(['\t mean_inchunk  = %f\n'], P.mean_inchunk);
            case 2 
                fprintf(['\t exp_pause    = %f %f %f\n'], P.pause_B,P.pause_A,P.pause_tau);
                fprintf(['\t exp_inchunk  = %f %f %f\n'], P.inchunk_B,P.inchunk_A,P.inchunk_tau);
        end; 
        fprintf(['\t rho           = %f\n'], P.rho);
        
        chunk_plot_fit(Data,Model,gamma,P); 
        drawnow;
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
        
        % Now calcualte for each observation, what the likelihood is that
        % it is a chunk boundary or within-chunk interval 
        boundary = diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0;
        p_boundary = bsxfun(@times,permute(gamma,[1 3 2]),permute(boundary,[3 2 1])); 
        p_boundary  = sum(p_boundary,3); 

        % reestimate MT parameters
        if (isfield(Data,'mt_seq'))
            % find mean_pause
            switch (Model.fit_rt)
                case 1          % Learn constant Pause model
                    [P.mean_pause, P.mean_inchunk,res1,res2] = learn_pause(Data.mt_seq, p_boundary);
                case 2   % Learn detrneding pause model
                    [P,res1,res2] = learn_exponential(Data,p_boundary,P);
            end;        
            P.v = learn_variance(res1,res2, p_boundary,Model.fit_rt_var);
            
            % reestimate covariance for movement time
            if Model.fit_rho
                P.rho = learn_cor(Model.chunks, Model.cor_chunks, res1,res2, gamma); 
            end
        end
        
        if (isfield(Data,'er_seq'))
            % find mean_pause
            if Model.fit_er
                [P.mean_pause_er, P.mean_inchunk_er,res1,res2] = learn_pause(Model.chunks, er_seq, gamma);
            end
            
            P.v_er = learn_variance(res1,res2, p_boundary);
            
            if P.fit_rho_er
                P.rho_er = learn_cor(Model.chunks, Model.cor_chunks, res1,res2, gamma); 
            end
        end
    end
end
P.log_like = log_like; 

