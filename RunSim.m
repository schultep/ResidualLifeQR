%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates the results for the simulations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expected Runtime: 5.3 hours - reduce MM to run smaller number of
%iterations;

N=[100,200,400,800,1600];    % N: sample sizes;
MM=1000;          % MM: Number of replications per sample size;
                  % Change MM to 1000 to mimic the simulations in the
                  % mansucript;
theta=[-.5,-.5,.01,-.5,-.5];    % Model parameters, assumed coefficients;
                                % For b1, intercept, slope, Z0, Z1 ;
rndEf=[0,1,1,1,.2];             % Model parameters, random effects.;
                                % Mean and variance of b0, mean and variance of b1, anc covariance;
tau=[.05,.10,.20,.8];           % tau - the quantiles of interest;
bw=.025;        % bw: a precision parameter

% Length of quantiles and sample sizes for looping ;
LL=length(tau);     NN=length(N);

% Declare matrices to store the simulation results ;
WQQR=zeros([LL,MM,NN]);  % results from training data;
QQR=zeros([LL,MM,NN]);   % results from test data ;



tic
for j=1:NN
    j
    for i=1:MM
        QR=simDat2(N(j),theta,rndEf);   % Training set
        QR2=simDat2(N(j),theta,rndEf);  % Test set of identical size
        MLE1=MLEfit_sim1(QR);   % Generate starting value for quantile regression:
                                % Ad hoc MLE fit to start the QR fit for each value of tau
        for k=1:LL
            QRF=QRfit_sim1(QR,MLE1,tau(k),bw); %Fit QR model for given tau
            [QQR(k,i,j)]=QRfit_sim1.testSet(QRF,QR2,MLE1); % See how the model did with the Test set
            [WQQR(k,i,j)]=QRfit_sim1.testSet(QRF,QR,MLE1); % See how the model did with the training set
        end
    end
end
toc 
   



% Define the i to get the results for the desired tau
i=2
tau(i)
"Quantile Regression results"
"Training Set"
N
"mean bias in quantile"
mean(squeeze(WQQR(i,:,:))-tau(i))
"root mean squared error"
sqrt(mean((squeeze(WQQR(i,:,:))-tau(i)).^2))
"mean absolute error"
mean(abs(squeeze(WQQR(i,:,:))-tau(i)))


"Test Set"
"mean bias in quantile"
mean(squeeze(QQR(i,:,:))-tau(i))
"mean absolute error"
mean(abs(squeeze(QQR(i,:,:))-tau(i)))
"root mean squared error"
sqrt(mean((squeeze(QQR(i,:,:))-tau(i)).^2))



