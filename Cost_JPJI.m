function [U,cost,cost0] = Cost_JPJI(Kx0,z,c,N,w,U_old,S,T)
C2zy = [];
C3zy = [];
C4zy = [];
MM = [];


for i=1:length(Kx0)
    k2 = Kx0(i);
    if (length(Kx0)==1)
        k3 = [k2,k2];
        k4 = [k2,k2,k2];
    elseif (length(Kx0)==2)
        k3 = [k2,Kx0(2)];
        k4 = [k2,Kx0(2),k2];
    else
        if i==length(Kx0)-1
            k3 = [k2,Kx0(i+1)];
            k4 = [k3,Kx0(1)];
        elseif i==length(Kx0)
            k3 = [k2,Kx0(1)];
            k4 = [k3,Kx0(2)];
        else
            k3 = [k2,Kx0(i+1)];
            k4 = [k3,Kx0(i+2)];
        end

    end
    %ind = ind+1;

    C2zy(:,i)=(z*S{k2}(c,:).')/(T);
    C3zy(:,i)=z*(S{k3(1)}(c,:) .*S{k3(2)}(c,:) ).'/T;
    C4zy(:,i)=...
        z*(S{k4(1)}(c,:) .*S{k4(2)}(c,:) .*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*S{k4(2)}(c,:)).')/T)*z(c,:) ...
        -diag(sum((S{k4(2)}(c,:).*z(c,:)).')/T)*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*z(c,:)).')/T)*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*S{k4(2)}(c,:)).')/T)*S{k4(3)}(c,:) ).'/T;   
    
end

MM(1:N,1:N,2)=w(2)*(C2zy*C2zy');
MM(1:N,1:N,3)=w(3)*(C3zy*C3zy');
MM(1:N,1:N,4)=w(4)*(C4zy*C4zy');

Mi = squeeze(sum(MM,3));

Czy(:,1)=Mi*U_old; 
cost=abs(U_old'*Czy);

[Q_L,~,Q_R]=svd(Czy,0);          % Optimize using thin-SVD.

U=Q_L*Q_R';
yy=U'*z;
yy = abs(yy);

%%
C2zy = [];
C3zy = [];
C4zy = [];
MM = [];
ind = 0;
for i=1
    k2 = Kx0(i);
    if (length(Kx0)==1)
        k2 = k2;
        k3 = [k2,k2];
        k4 = [k2,k2,k2];
    elseif (length(Kx0)==2)
        k2 = k2;
        k3 = [k2,Kx0(2)];
        k4 = [k2,Kx0(2),k2];
    else
        k2 = Kx0(i);
        if i==length(Kx0)-1
            k3 = [k2,Kx0(i+1)];
            k4 = [k3,Kx0(1)];
        elseif i==length(Kx0)
            k3 = [k2,Kx0(1)];
            k4 = [k3,Kx0(2)];
        else
            k3 = [k2,Kx0(i+1)];
            k4 = [k3,Kx0(i+2)];
        end

    end
    ind = ind+1;

    C2zy(:,ind)=(z*S{k2}(c,:).')/(T);
    C3zy(:,ind)=z*(S{k3(1)}(c,:) .*S{k3(2)}(c,:) ).'/T;
    C4zy(:,ind)=...
        z*(S{k4(1)}(c,:) .*S{k4(2)}(c,:) .*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*S{k4(2)}(c,:)).')/T)*z(c,:) ...
        -diag(sum((S{k4(2)}(c,:).*z(c,:)).')/T)*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*z(c,:)).')/T)*S{k4(3)}(c,:) ...
        -diag(sum((S{k4(1)}(c,:).*S{k4(2)}(c,:)).')/T)*S{k4(3)}(c,:) ).'/T;    
end

MM(1:N,1:N,2)=w(2)*(C2zy*C2zy');
MM(1:N,1:N,3)=w(3)*(C3zy*C3zy');
MM(1:N,1:N,4)=w(4)*(C4zy*C4zy');
Mi=squeeze(sum(MM,3));

Czy(:,1)=Mi*U_old; 
cost0=abs(U_old'*Czy);


