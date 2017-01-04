function v = learn_variance(res1,res2, p_boundary,fit_rt_var)
% function v = learn_variance(res1,res2, p_boundary)
% learn variance similar
% INPUT:
%      res1:   Residuals under assumption that it is chunk boudnary
%      res2:   Resdiaul under assumption that is is within chunk
%      p_boundary: probability of being on the chunk boudnary
%      fit_rt_var: 1: fir one single variance 2: fit seperate variances

indx=~isnan(res1);
d1=res1(indx);
d2=res2(indx);
w=p_boundary(indx);
switch(fit_rt_var)
    case 1
        v = sum((d1.^2).*w+(d2.^2).*(1-w))/sum(indx(:));
    case 2
        v(1,1) = sum((d1.^2).*w)/sum(w);
        v(1,2) = sum((d2.^2).*(1-w))/sum(1-w);
end;
