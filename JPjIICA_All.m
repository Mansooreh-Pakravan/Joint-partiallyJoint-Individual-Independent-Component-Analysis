function [Sjtica,JpJIFAll,Fratio] = JPjIICA_All(Xin,Dim,Statistics,verbose,MaxIter,mythr)

%% ------------Joint/partially-joint/Individual Independent Component Analysis---------- 
% Mansooreh Pakravan & Mohammad Bagher Shamsollahi----------------------

%Please cite the paper:
%"Joint, Partially-joint, and Individual Independent Component Analysis in Multi-Subject fMRI Data"
% published in Transaction on biomedical engineering

%% The code is inspired by:
%S. Cruces, A. Cichocki, "Combining blind source extraction and 
%     joint approximate diagonalization", Proc. fourth symp. on ICA/BSS,
%     pp. 463--468, Japan, 2003.

%% Dimension reduction using PCA
% Suppose we have estimated correctly the number of all of sources for each
% subject using a model order selection method-> stored in Dim.Cest \n R^(1\times Dim.K)
Ytensor = Data_ItsPCA_DiffC(Dim,Xin,Dim.Cest);  

%% Initialization : Identity matrix
Sjtica = cell(1,Dim.K);
for k=1:Dim.K
    Sjtica{k} = Ytensor{k}(1:Dim.Cest(k),:);
end

JpJIFAll = zeros(max(Dim.Cest),Dim.K)-1;
Fratio = zeros(max(Dim.Cest),Dim.K)-1;

%% Outer iteration for each source and each subject
for iter = 1:MaxIter 
    Ytensor0 = Ytensor;
    for k = 1:Dim.K
        for c = 1:Dim.Cest(k)
            
            %% The main Cost function
            % Note that in the stop rule of the inner iterations, the parameter \epsilon
            % is set to 10^(-7), so it maybe alittle slow.

            [B,costAll,Ratio] = JPjIICA(Ytensor0,Sjtica,k,c,Dim,Statistics);

            JpJIFAll(c,k) = costAll(end); 
            Fratio(c,k) = Ratio(end);
            
            %% If the source is indiviudal, then:
            if JpJIFAll(c,k)<mythr
                [B,~,~] = JPjIICA(Ytensor0,Sjtica,k,c,Dim,Statistics); 
            end
            %% Compute the desired source from estimated demixing matrix 
            Sjtica{k}(c,:) = B*Ytensor0{k};
            
            %% Deflation : remove the extracted source from the contents of
            % the observation
            Sk = Sjtica{k}(c,:);
            for t = 1:size(Ytensor0{k},1)
                y = Ytensor0{k}(t,:);
                y = y-mean(y.*Sk)*(1/mean(Sk.*Sk))*Sk;
                Ytensor0{k}(t,:) = y;
            end
        end                               
    end
 
end
%% plot the extracted sources
if verbose
    strN = ['JpJI-ICA Results'];
    figure('Name',strN,'Visible','On','NumberTitle','off') 
   for k=1:Dim.K
        i=(k-1)*max(Dim.Cest)+1;
        for c=1:Dim.Cest(k)
            if i-(k-1)*max(Dim.Cest)<=Dim.Cest(k)+k*max(Dim.Cest)
                subplot(Dim.K,max(Dim.Cest),i)
                temp = abs(Sjtica{k}(c,:));
                temp = temp./max(temp);
                imagesc(reshape(temp,Dim.V,Dim.V))
                 axis off
                drawnow
                i=i+1;
            end
        end
    end

end 


    

