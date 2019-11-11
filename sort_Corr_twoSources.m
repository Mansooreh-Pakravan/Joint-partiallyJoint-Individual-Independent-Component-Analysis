function [Snew,indsort] = sort_Corr_twoSources(S1,S2)

X = abs(S2);
Snew = zeros(size(S1));
indsort = zeros(size(S1,1),1);
 for n=1:size(S1,1)
        R = corr((S1(n,:))',X');
        mR = max(R);
        imR = find(R == mR);
        ind=imR(1);
        Snew(n,:) = X(ind,:);
        X(ind,:) = zeros(1,size(S1,2));
        indsort(n) = ind;
 end
% 
% for n=1:size(S1,1)
%     res0 = zeros(1,size(S1,1));
%     for c=1:size(S1,1)
%         vec = S1(n,:);
%         vecr = X(c,:);   
% %         vec = vec./norm(vec);
% %         vecr = vecr./norm(vecr);     
%         res0(c) = norm(vec-vecr);
%     end
%     n
%     res0
% 
%     
%     mR = min(res0);
%     imR = find(res0 == mR);
%     ind=imR(1);
%     Snew(n,:) = X(ind,:);
%     X(ind,:) = inf(1,size(S1,2));
%     indsort(n) = ind;
% end