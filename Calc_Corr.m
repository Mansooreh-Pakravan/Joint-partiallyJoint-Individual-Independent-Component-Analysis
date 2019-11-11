function [corr_out,sir] = Calc_Corr(Sr,Ss,C,verbose)

corr_out = zeros(C,1);
sir = zeros(C,1);
[Snew,~] = sort_Corr_twoSources(Sr,Ss);
for c=1:C
    vec = Snew(c,:);
    vecr = Sr(c,:);   
    vec = vec./norm(vec);
%     vec = vec./max(vec);
    vecr = vecr./norm(vecr);     
%     vecr = vecr./max(vecr);    
    corr_out(c,1) = corr(vec',vecr');
    corr_out(c,1) = norm(vec'-vecr');
    sir(c,1)=CalcSIR(vec',vecr');
end
corr_out = mean(corr_out);
sir = mean(sir);
V = sqrt(length(vec));
if verbose
    figure
    
    for c=1:C
        subplot(C,1,c)
        vec = Snew(c,:);
        vecr = Sr(c,:);   
        vec = vec./norm(vec);
        vec = vec./max(vec);
        vecr = vecr./norm(vecr);     
        vecr = vecr./max(vecr); 
%         imagesc(reshape(vec,V,V))
        plot(vecr,'g','LineWidth',2)
        hold on
        plot(vec)

    end
end