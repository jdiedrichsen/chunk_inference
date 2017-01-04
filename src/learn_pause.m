function [pause, inchunk,res1,res2] = learn_pause(data, p_boundary)
% function [start_pause, nonstart_pause] = learn_pause(data, p_boundary)
% learn pause at the start of chunk and pause not at the start
% based on probability of chunks in matrix gamma for the data `data`
% indicator variable for when chunk starts
% Ignores NaNs in the data 
% Parameters: 
%   data:        Data to learn from (nTrials x nMovements) 
%   p_boundary:  Posterior probability of the hidden state (nTrials x nChunks) 
% Returns: 
%   pause:      Mean Value for beginning of chunk 
%   inchunk:    Mean Value for within chunk 
%   res1:       Residual, assuming pause
%   res2:       REsidual, assuming withinchunk 

% Mean for the chunk having started 
indx=~isnan(data);
d=data(indx); 
w=p_boundary(indx); 
pause = sum(d.*w)/sum(w);
res1 = data-pause; 

% Mean for within chunk 
w=1-p_boundary(indx); 
inchunk = sum(d.*w)/sum(w);
res2 = data-inchunk; 
