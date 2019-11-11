function Main_function(Xin,Sclear,Dim)

%% Inputs are the source and the obsercation of a multi-subject data according to the JpJI-MDU source model 
% 'Sclear': the original clear source: the sources are stored in a cell[1
% \in Dim.K] and the $k$th cell belogs to $k$th subject. Sclear{k,1} has
% the dimension Dim.C (number of sources) \times Dim.V3 (nnumber of all voxels: for example the vectorized version of fmri spatial source in the 3D volumn )

% Note that we need Sclear as the gold stundard only for validating the
% algorithm (in function Determine_Types_Evaluate_Results.m).
%If you do not have the Sclear, you can call only JPjIICA_All.m function.

% 'Dim' : a struct :the dimension of the simulated data
% 'Dim.K' : the number of subjects
% 'Dim.C1' : the number of joint sources
% 'Dim.C2' : the number of partially-joint sources (in this version of published code C_2 is equal for all of the subjects)
% 'Dim.C3k' : it is a vector that its 'k'th element is equal to C_3^{(k)}

% 'Xin': Observations are stored in a cell[1 \in Dim.K] and each cell
% has the dimension [Dim.T \times Dim.V3] (as the observations of each subject)

%% plot original sources (If you have the gold standard)
verbose = 1;
if verbose
    strN = 'Original Sources';
    figure('Name',strN,'Visible','On','NumberTitle','off')
    for k=1:Dim.K
        i=(k-1)*max(Dim.Cest)+1;
        for c=1:Dim.Cest(k)
            if i-(k-1)*max(Dim.Cest)<=Dim.Cest(k)+k*max(Dim.Cest)
                subplot(Dim.K,max(Dim.Cest),i)
                temp = Sclear{k}(c,:);
                temp = temp./max(temp);
                imagesc(reshape(temp,Dim.V,Dim.V))
                 axis off
                i=i+1;
            else        
            end
        end
    end
end

%% set the parameters of the JpJI-ICA algorithm
w = [0 2 3 4];
Statistics=w/max(w);   % Weights of each cumulant order.
MaxIter = 10;
 
%% The Joint/partially-Joint/Individual Independent Component Analysis algorithm
mythr = 1; % you can set the threshold with the function "Determine_Types_Evaluate_Results" automatically using your train data.
[Sjtica,JpJIF,flag] = JPjIICA_All(Xin,Dim,Statistics,verbose,MaxIter,mythr);  

%% Determining the \sigma_opt automatically, Estimating the number of joint, partiall-joint,...
%% and individual sources, computing jSIR and clustering...
%% the subjetcs based on their partially-joint sources
[Sigma_opt,sirs,CostJ,CostpJ,CostI,EstC1,EstC2,EstC3] = Determine_Types_Evaluate_Results(Dim,flag,JpJIF,Sclear,Sjtica);
