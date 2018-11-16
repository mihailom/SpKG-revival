%{
Knowledge Gradient for Sparse Linear Regression
This is used to provide suggestions to cover 100% of the target molecule 

1/2 \sum(Loss)+\lambda*penalty.

notation for the following:
p:   number of features.
number of nucleotides
M:   number of alternatives.
number of regions


This function takes in
T:        target molecule
X:        the design matrix for linear regression model (M x p)  
IN:       index set for all alternatives (M x 2)
theta_0:  prior for the linear regression parameters (p x 1)
MVar:     measurement variance (scalar)
a good estimate of the variance -- variance of energies
C:        prior covariance matrix for theta
along diagonals are variances
sq rt diagonal will tell us the confidence on the theta_0
lambda:   regularization parameter
calculated via cross validation or validation??
units are squared energy

options:  a structure containing
 epsilon:  calculation precision
 gamma:    the forgetting factor for least squares, default is 1 
 xi,eta:   the parameters for beta priors
 flipN:    the number of fliping groups
 cases:    number of MC simulation to estimate lasso estimator covariance matrix

And returns
choices:  alternatives picked at each iteration 
KGVals:   KG values at each iteration
kg_vals:  KG value of the best alternative
%}

function [choices, KGVals, kg_vals]= KGSpCover(T,X,IN,theta_0,MVar,C,lambda,options)


% check and initialize
p      = length(T);
G      = eye(p);
w      = 1;

if nargin<8
  options=[];
end
if(~isfield(options, 'epsilon'))
  options.epsilon = 1e-6;
end
if(~isfield(options, 'xi'))
  for i = 1:p
    options.xi(i) = w*logical(any(theta_0(find(G(i,:)))));
  end
  options.xi = options.xi+ones(1,p);
end
if(~isfield(options, 'flipN'))
  options.flipN = min(4,p);
end
if(~isfield(options, 'sampleSize'))
  options.sampleSize = 1;
end

% check if flipN is smaller than p
if options.flipN > max(p,32) 
    error('The number of flip groups is larger than number of overall groups or it may overflow!');
end

% check if lambda is a vector of length p
if length(lambda)==1
    lambda = lambda.*ones(p,1);
else if length(lambda)~=p
        error('Input dimension of lambda vector may be wrong!');
    end
end



xi                = options.xi;
flipN             = options.flipN;
sampleSize        = options.sampleSize;
eta               = w+2-logical(xi);

% theta_0, xi,eta has to be a column vector
theta_0      = theta_0(:);
xi           = xi(:);
eta          = eta(:);

theta        = theta_0;
KGVals       = [];
kg_vals      = [];


%initialize
previousChoices = [];
cover           = zeros(p,1);
M               = size(X,1);
count           = 0;
 while(any(cover == 0))
     count = count+1;
     % choose using sparse KG, get batch decisions
     unnorm_KGVal = KGSpLin_batch(theta,C,X,G,MVar,xi,eta,flipN,previousChoices,sampleSize);
     KGVal        = (unnorm_KGVal - min(unnorm_KGVal))./(max(unnorm_KGVal)-min(unnorm_KGVal)); 
      % compute cost coefficient for each alternative
      cost         = zeros(M,1);
      for i = 1:M
          cost(i)         = (lambda - KGVal(i))/length(setdiff(IN(i,1):IN(i,2),find(cover)));
      end
      cost(previousChoices) = Inf;
      cost(cost==-Inf) = Inf; 
      [~,bestX]    = min(cost);

     % update
     previousChoices = [previousChoices,bestX];
     KGVals          = [KGVals,KGVal]; 
     kg_vals         = [kg_vals,KGVal(bestX)];
     cover(IN(bestX,1):IN(bestX,2)) = 1;
 end
choices = previousChoices; 

end

%KGSpLin_batch returns the KGVals for all alternatives
%X:     matrix for alternatives (M x p)
%theta: belief about parameters (p x 1)
%C:     covariance for parameters  
function KGVal = KGSpLin_batch(theta,C,X,G,MVar,xi,eta,flipN,previousChoices,sampleSize)

    [M, p]           = size(X);
    KGMat            = zeros(M,2^flipN-1);
    probs            = zeros(1,2^flipN-1);
    allindices       = 1:p;
    Mean             = xi./(xi+eta);
    MeanC            = 1-Mean;
    AbsDev           = abs(Mean-0.5);
    [~, ind]         = sort(AbsDev);
    flipind          = ind(1:flipN);
    ins              = setdiff(find(Mean>=0.5),flipind);  % current index of "in" features
    
    for j = 1:2^flipN-1
        bits         =  dec2bin(j,flipN)-'0';
        curr_ins     =  union(ins, flipind(find(bits)));
        curr_outs    =  setdiff(allindices,curr_ins);
        zeta         =  sum(G(curr_ins,:),1);
        KGMat(:,j)   =  BKG_RNA(theta,C,MVar,X,zeta,previousChoices, sampleSize);
        
        probs(j)     =  prod([Mean(curr_ins);MeanC(curr_outs)]);
        
    end
    probs           = probs./sum(probs);
    KGVal           = sum(repmat(probs,M,1).*KGMat,2);

end

