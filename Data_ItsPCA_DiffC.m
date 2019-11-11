function Ytensor = Data_ItsPCA_DiffC(Dim,Xnoise,DimCk)


for k=1:Dim.K
    C = DimCk(k) ;
    Data0 = Xnoise{k};
    [Coef,score,~,~] = pca(Data0','Centered','off','NumComponents',C);
    Ytensor{k} = score';
         
end