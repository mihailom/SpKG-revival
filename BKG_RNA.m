function KG = BKG_RNA(theta_n, C, MVar, X, zeta, previousChoices, sampleSize)

%theta_n: column vector of estimated means at the n-th iteration
%C: covariance matrix at the n-th iteration
%MVar: variance of measurement noise
%X: design matrix M*K, with M the number of alternatives and K the number
%of features
%zeta: group structure
%previousChoises: previous choises with the batch, initially all zeros
%sampleSize: sample batch size

% obtain active alternatives
a_theta=theta_n(zeta==1);
a_C=C(zeta==1,zeta==1);
a_X=X(:,zeta==1);

%bs: the bs-th decision with a batch
bs=length(previousChoices(previousChoices~=0))+1;

[M,K]=size(a_X);
sigmatilde=zeros(M,bs);
%calcultate the maximal knowledge gradient value of bs alternatives under the i-th choice
%the first choice is found using exact KG fomular
%recuresively find the following alternatives using MC simulatoin
for j=1:bs-1
    % update statistics within the batch
    sigmatilde(:,j)=sigmaTilde(a_X,a_C, previousChoices(j), MVar);
    gammaU=MVar+a_X(previousChoices(j),:)*a_C*a_X(previousChoices(j),:)';
    a_C=a_C-(1./gammaU)*a_C*a_X(previousChoices(j),:)'*a_X(previousChoices(j),:)*a_C;
end
KG=FND(bs-1, a_X, sigmatilde, a_theta, a_C, MVar,sampleSize);    

end



function KG = FND(nbd, X, sigmatilde, theta, C, MVar,sampleSize)
%nbd: number of decisions made in this batch
%i: current choice in the other dimensions 
%M2: number of choices in the batch dimension
  [M,K]=size(X);
  KG=zeros(M,1);
  a=X*theta;
  for iter=1:M
       b=sigmaTilde(X, C, iter, MVar);
       if nbd==0      %find the fist alternative using KG formular
           [KGf,LogKG] = LogEmaxAffine(a',b');      
           KG(iter)=KGf;
       else
           %Monte Carlo Sampling
             sigmatilde(:,nbd+1)=b;
             KGtemp=sum(max(sigmatilde(:,1:nbd+1)*randn(nbd+1,sampleSize)+a*ones(1,sampleSize)))/sampleSize;
             KG(iter)=KGtemp;      
       end
  end
end


function sigmatilde=sigmaTilde(X, C, i, MVar)
  B=X*C;
  Sigma=B*(X(i,:)');
  sigmatilde=Sigma./sqrt(MVar+Sigma(i));
end

        
        
        
        
        
        
        
