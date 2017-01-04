function [P,res1,res2] = learn_exponential(Data, p_boundary,P); 
% function [P,res1,res2] = learn_exponential(Data, p_boundary,P); 
% learn pause at the start of chunk and pause not at the start
% based on probability of chunks in matrix gamma for the data `data`
% indicator variable for when chunk starts
% Ignores NaNs in the data 
% Parameters: 
%   data:   Data to learn from, field mt_seq and trialNum  
%   p_boundary:  Posterior probability of whether it is a boundary or not (nTrials x nChunks) 
% Returns: 
%      P: parameters for exponentials for within and between chunks intervals
%      res1: Residuals under chunk boundary assumption 
%      res2: Resdiual under within chunk assumption 

indx = ~isnan(Data.mt_seq); 
y = Data.mt_seq(indx); 
t = repmat(Data.trialNum,1,size(Data.mt_seq,2)); 
t = t(indx);
w = p_boundary(indx); 

[P.pause_B,P.pause_A,P.pause_tau]=fit_exponential(t,y,w,P.pause_B,P.pause_A,P.pause_tau);
[P.inchunk_B,P.inchunk_A,P.inchunk_tau]=fit_exponential(t,y,1-w,P.inchunk_B,P.inchunk_A,P.inchunk_tau);

yp1 = (P.pause_B-P.pause_A)*exp(-t/P.pause_tau)+P.pause_A;
yp2 = (P.inchunk_B-P.inchunk_A)*exp(-t/P.inchunk_tau)+P.inchunk_A;

res1 = nan(size(Data.mt_seq)); 
res1(indx) = y-yp1;  
res2 = nan(size(Data.mt_seq)); 
res2(indx) = y-yp2; 



function [B,A,tau]=fit_exponential(t,y,w,B,A,tau); 
x0=[B;A;tau]; 
OPT=optimset(@fmincon); 
OPT.Display='off'; 
fcn=@(x) weighterror(x,t,y,w); 

x = fmincon(fcn,x0,[],[],[],[],[-inf -inf 0.1],[],[],OPT);
B=x(1); 
A=x(2); 
tau=x(3); 

function err=weighterror(x,t,y,w); 
yp  = (x(1)-x(2))*exp(-t/x(3))+x(2); 
err=sum((y-yp).*(y-yp).*w); 
