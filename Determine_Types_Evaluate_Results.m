function [Sigma_opt,sirs,CostJ,CostpJ,CostI,EstJ,EstPJ,EstI,Acc]=Determine_Types_Evaluate_Results(Dim,flag,JpJIF,S,Sjtica)

for k=1:Dim.K
    [corrcn(k), sirsn(k)] = Calc_Corr(abs(S{k}),abs(Sjtica{k}),Dim.Cest(k),0);
end       
sirs = mean(sirsn);
corrr = mean(corrcn);

for k=1:Dim.K
    for c=1:min(Dim.Cest)
        f0 = flag(c,k);
        flag0(c,k) = abs(f0-(Dim.K-1));
    end
end

c01 = 1;
Cj=[];
EstJ=0;
for c=1:min(Dim.Cest)
    idf = length(find(flag0(c,:)<2));
    if idf>=Dim.K-2
        EstJ = EstJ+1;
        Cj(c01) = c;
        c01 = c01+1;
    end
end

a = JpJIF(Cj,:);
CostJ = mean(a(:));

JpJIF0 = JpJIF;
JpJIF0(Cj,:)=[];

c0=1;
for c=EstJ+1:max(Dim.Cest)
    f = JpJIF(c,:);
    f(f==-1)=[];
    if mean(f)==0
        jf(c0) = 0;
    else
        jf(c0) = CostJ./mean(f);
    end
    c0=c0+1;
end

label = kmeans(jf',2);

ind1 = find(label==1);
ind2 = find(label==2);

ratioclass1 = jf(ind1);
ratioclass2 = jf(ind2);


featureclass1 = JpJIF0(ind1,:);
featureclass2 = JpJIF0(ind2,:);

if min(ratioclass1)>max(ratioclass2)
    % class 1 is individual
    EstPJ = length(ratioclass2);
    lowBestTHR = max(featureclass1(:));
    highBestTHR = min(featureclass2(:));
    CostpJ = mean(featureclass2(:));
    CostI = mean(featureclass1(:));
else
    % class 2 is individual
    EstPJ = length(ratioclass1);
    lowBestTHR = max(featureclass2(:));
    highBestTHR = min(featureclass1(:));
    CostpJ = mean(featureclass1(:));
    CostI = mean(featureclass2(:));
end

Sigma_opt = (lowBestTHR+highBestTHR)/2;

for k=1:Dim.K
    EstI(k) = Dim.Cest(k)-(EstJ+EstPJ);
end

