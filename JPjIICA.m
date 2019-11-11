function [B,costAll,Ratio] = JPjIICA(Ytensor0,S,kkk,ccc,Dim,Statistics)

x = Ytensor0{kkk};
[N,T]=size(x);

MaxInnerIter     =50;
StopRule        =1e-7;  

w=Statistics/max(Statistics);   % Weights of each statistic.

%% PREWHITENING 
x=x-mean(x')'*ones(1,T);
Rxx=x*x'/T;
[Q,D]=eig(Rxx);
[v,i]=sort(diag(D));
vx=flipud([v,i]);
vx(vx==0)=(10^-20);
W = diag(1./sqrt(vx(1:N,1)))*Q(:,vx(1:N,2))';
z=W*x;                 

%% INITIALIZATION the raws of the demixing matrix
U=eye(size(z,1),1);
yy=U'*z;

%% Defining Ktild
Kx0 = 1:Dim.K;
Kx0(kkk) = [];
Kx0 = Kx0(randperm(Dim.K-1,Dim.K-1));

NotConvergence=1;
costAll = [];
cost0 = [];

inneriter=0;
while (inneriter<MaxInnerIter && NotConvergence)
    inneriter=inneriter+1;
    U_old=U;
    if ccc<=min(Dim.Cest)  % because in this version of published code we assume C_1+C_2 is equal for all of the subjects.
        [U,FjointALL,Fjoint0] = Cost_JPJI(Kx0,z,ccc,N,w,U_old,S,T);            
        costAll(inneriter) = FjointALL;
        cost0(inneriter) = Fjoint0;
        Ratio(inneriter) = FjointALL/Fjoint0;
    else
        [U,~,~] = Cost_JPJI(kkk,z,ccc,N,w,U_old,S,T);  
        costAll(inneriter)=0;
        cost0(inneriter)=0;
        Ratio(inneriter)=0;
        Mu=[];
    end

    %% Stop rule
    if inneriter>3
        NotConvergence=max(abs(U(:)-U_old(:)))>StopRule;
    end
end 

%% Demixing matrix.
B=U'*W;  
