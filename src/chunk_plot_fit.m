function chunk_plot_fit(Data,Model,gamma,P); 
% Plots the fit of the current
% Model as a function of trialNumber and sequenceType 


% Determine whether data point is more likely classsified as 
% Between chunk or within chunk itnerval 
boundary = diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0;
p_boundary = bsxfun(@times,permute(gamma,[1 3 2]),permute(boundary,[3 2 1])); 
p_boundary  = sum(p_boundary,3); 

t = repmat(Data.trialNum,1,size(Data.mt_seq,2)); 
MTd  = floor(Data.mt_seq(:)/25)*25; 
row = unique(MTd(~isnan(MTd))); 
col  = unique(Data.trialNum(:)); 
Data1 = pivottable(MTd,t(:),p_boundary(:),'nansum','forcerow',row,'subset',~isnan(MTd)); 
Data2 = pivottable(MTd,t(:),1-p_boundary(:),'nansum','forcerow',row,'subset',~isnan(MTd)); 
Data1(isnan(Data1))=0; 
Data2(isnan(Data2))=0; 
sc=max([Data1(:);Data2(:)]); 
C(:,:,3)=1-Data1/sc; 
C(:,:,1)=1-Data2/sc; 
C(:,:,2)=1-(Data1/sc + Data2/sc)/2; 
image(col,row,C); 
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