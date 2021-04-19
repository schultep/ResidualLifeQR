classdef QRfit_sim1
    properties   
        theta_Z  % static covariate value
        theta_X  % longitudinal covariate value
        theta_t  % constant, and linear time variable.
        val      % Value of the likelihood cost function
        FLAG     % Did the optimization complete    
        MLE_Cost
        QR_Cost
        MLE_Perc 
        QR_Perc
        QR_Perc2
        wtP
        MLE_Quant
        QR_Quant
        tau
        bw
    end
    methods
       function fit=QRfit_sim1(DAT,MLEfit,tau,bw)
           if nargin<3
               return
           end
           fit.tau=tau;
           if nargin<5
               bw=.05;
           end
      
           fit.bw=bw;
           %starting values
           th=[MLEfit.theta_X, MLEfit.theta_t,MLEfit.theta_Z];
           [fit.MLE_Cost,fit.MLE_Quant]=QRfit_sim1.gee(th,tau,DAT,bw);
 %           fit.MLE_Perc=QRfit_sim1.actPerc(fit.MLE_Quant,DAT,DAT.trueB,DAT.trueTheta);
          
           cst=@(TTT) sum(abs(QRfit_sim1.gee(TTT,tau,DAT,bw)));
           for i=1:5     
               [th,fit.val,fit.FLAG]=fminsearch(cst,th,optimset('Display','off'));
           end
           
           [fit.QR_Cost,fit.QR_Quant]=QRfit_sim1.gee(th,tau,DAT,bw);
    %       [fit.QR_Perc,fit.wtP]=QRfit_sim1.actPerc(fit.QR_Quant,DAT,DAT.trueTheta);

           
           fit.theta_X=th(1); 
           fit.theta_t=th(2:3);
           fit.theta_Z=th(4:end);
           return
       end
       %methods end
    end
    methods (Static)
        function [ZZ,Q] = gee(theta,tau,DAT,bw)
            Q=QRfit_sim1.gorfQ(theta,tau,DAT);
            u=(DAT.Y-Q);            
            a=1+exp(-u/bw);
            a=a.^(-1);
            a(u<-(100*bw))=0;
            a(u>(100*bw))=1;
            tmp=DAT.wt.*(a-1+tau);  
            
            ZZ=[mean(tmp)/mean(DAT.wt);(DAT.Z'*tmp)./sum(DAT.Z,1)';
                mean(DAT.X.*tmp)/mean(DAT.X);mean(DAT.T.*tmp)/mean(DAT.T)];%.2*mean(a(unique(DAT.id,'first'))-1+tau)];
        end
        
        function [Q] = gorfQ(theta,tau,DAT)
           cT= -log(1-tau);
           MX=5+DAT.T;    %2*max(DAT.Y);
           const=exp(theta(2)+DAT.Z*theta(4:end)'+theta(1)*DAT.b_QR);  
           inner=theta(3)*cT./const+exp(theta(3)*DAT.T);
           ind=inner<=0;
           Q=log(inner)/theta(3);
           Q(ind)=MX(ind);
           Q(Q>MX)=MX(Q>MX);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % actPerc
        % Return - Actual percentiles
        % Arguments - 
        %   Fit - quantile fits through GorfQ
        %   DAT - QRDat object
        %   b - Actual values of b for each id
        %   thR - Real theta values
        %%%%%%%%%%%%%%%%%%%%%%%
        function [P,wtP] = actPerc(Fit,DAT,thR)
          % Fit=QRF.QR_Quant;
           b=DAT.trueB(DAT.id);
           const=exp(thR(2)+DAT.Z*thR(4:end)'+thR(1)*b);
           if thR(3)==0
                P= 1-exp(const.*(DAT.T-Fit));
           else
               P= 1-exp(const/thR(3).*(exp(thR(3)*DAT.T)-exp(thR(3)*Fit)));
           end
            tmp=DAT.wt.*P;  
            wtP=mean(tmp)/mean(DAT.wt);
            
            
            
        end
                
        function [p_vals] = boot(QRfit,DAT,NN)
            lTH=3+size(DAT.Z,2);
            mQR=zeros([lTH,NN]);
            mMLE=zeros(lTH);
            tau=QRfit.tau; bw=QRfit.bw;
            for i=1:NN
                DT1=QRdata.Boot(DAT);
                MLE1=MLEfit_sim1(DT1);
                th=[MLE1.theta_X, MLE1.theta_t,MLE1.theta_Z];
                cst=@(TTT) sum(abs(QRfit_sim1.gee(TTT,tau,DAT,bw)));
                for j=1:5     
                    [th,fit.val,fit.FLAG]=fminsearch(cst,th,optimset('Display','off'));
                end
                 mQR(:,i)=th;
            end
            %95% CI based on Normal approximation
            resA=mQR<=0;
            resB=mQR>=0;
            
            p_vals=2*min(mean(resA'),mean(resB'));
%             sQR=std(mQR'); 
%             mnQR=mean(mQR');
%             tmp =    abs(mnQR./sQR);
%             p_vals=2*(1-normcdf(tmp));
            
            
            
        end
        
        
        function [QRquant,MLEquant,qQR] = testSet(QRfit,DAT,MLEfit)           
            thMLE=[MLEfit.theta_X, MLEfit.theta_t,MLEfit.theta_Z];
            thQR=[QRfit.theta_X, QRfit.theta_t,QRfit.theta_Z];
            qMLE=QRfit_sim1.gorfQ(thMLE,QRfit.tau,DAT);
            qQR=QRfit_sim1.gorfQ(thQR,QRfit.tau,DAT);
            MLEquant=mean(DAT.wt.*(qMLE>=DAT.Y))/mean(DAT.wt);
            QRquant=mean(DAT.wt.*(qQR>=DAT.Y))/mean(DAT.wt);
            return
        end
        
    end
    %classdef ends
end