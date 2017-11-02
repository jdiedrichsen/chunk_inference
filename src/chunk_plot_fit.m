function chunk_plot_fit(Data,Model,gamma,P,varargin); 
% Plots the fit of the current
% Model as a function of trialNumber and sequenceType 

colordensity = 1; % Color the density plot by within / between ? 

vararginoptions(varargin,{'colordensity'}); 

% Determine whether transition is more likely within or between
% (p_boundary) 
boundary = diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0;
p_boundary = bsxfun(@times,permute(gamma,[1 3 2]),permute(boundary,[3 2 1])); 
p_boundary  = sum(p_boundary,3); 

% Get densities in bins of 25ms 
t = repmat(Data.trialNum,1,size(Data.mt_seq,2)); 
MTd  = floor(Data.mt_seq(:)/25)*25; 
row = unique(MTd(~isnan(MTd))); 
col  = unique(Data.trialNum(:)); 
Data1 = pivottable(MTd,t(:),p_boundary(:),'nansum','forcerow',row,'subset',~isnan(MTd)); 
Data2 = pivottable(MTd,t(:),1-p_boundary(:),'nansum','forcerow',row,'subset',~isnan(MTd)); 
Data1(isnan(Data1))=0; 
Data2(isnan(Data2))=0; 
Data1(Data1>0)=Data1(Data1>0)+3; % Add intercept for better visibility  
Data2(Data2>0)=Data2(Data2>0)+3; 


% Now determine the color of each bin 
countd =(Data1 + Data2);  % Overall density if the data 
if (colordensity) 
    sc=max([Data1(:);Data2(:)]); 
    C(:,:,1)=1-Data2/sc;   % red for within 
    C(:,:,2)=1-max(Data1/sc,Data2/sc); 
    C(:,:,3)=1-Data1/sc; 
    C(C>1)=1; 
else 
    sc = max(countd(:));   % Gray scale for count density 
    C(:,:,1)=1-countd/sc;   
    C(:,:,2)=1-countd/sc; 
    C(:,:,3)=1-countd/sc; 
end; 
imagesc_rectangle(col,row,C); 


set(gca,'YDir','normal'); 
xlabel('trial Number'); 
ylabel('MT'); 

switch(Model.fit_rt)
    case 1
        pred_pause  =  ones(length(col),1).*P.mean_pause; 
        pred_inchunk  =  ones(length(col),1).*P.mean_inchunk; 
    case 2 
        pred_pause  =  (P.pause_B-P.pause_A)*exp(-col/P.pause_tau)+P.pause_A;
        pred_inchunk = (P.inchunk_B-P.inchunk_A)*exp(-col/P.inchunk_tau)+P.inchunk_A;
end; 

hold on; 
plot(col,pred_pause,'r',col,pred_inchunk,'b');
hold off; 