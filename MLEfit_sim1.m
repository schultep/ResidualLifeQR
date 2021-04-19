classdef MLEfit_sim1
    properties   
        theta_Z  % static covariate value
        theta_X  % longitudinal covariate value
        theta_t  % constant, and linear time variable.
        val      % Value of the likelihood cost function
        FLAG     % Did the optimization complete   
    end
    methods
       function fit=MLEfit_sim1(DAT)
           %starting values
           vls=coxphfit(DAT.Z,DAT.Y, ... 
               'Censoring',DAT.delta);
           th=zeros([1,(3+size(DAT.Z,2))]);            
           th(4:end)=vls;
           MLE=@(TTT) -MLEfit_sim1.MLE_cost(TTT,DAT);
           for i=1:3  
               [th,fit.val,fit.FLAG]=fminsearch(MLE,th,optimset('Display','off'));
           end
           
           fit.theta_X=th(1); 
           fit.theta_t=th(2:3);
           fit.theta_Z=th(4:end);
           return
           
         
       end
       %methods end
    end
    methods (Static)
       function [LogLik] = MLE_cost(theta,DAT)
           %change so that all vectors are the size of the number of id's
           const=exp(theta(2)+DAT.Z_MLE*theta(4:end)'+theta(1)*DAT.b_MLE');
           Y=DAT.Y_MLE; delta=DAT.delta_MLE;
           h=theta(2)+DAT.Z_MLE*theta(4:end)'+theta(1)*DAT.b_MLE'+theta(3)*Y;
           
           if theta(3)==0
               H= -const.*Y; 
           else
               H= -const.*(exp(theta(3)*Y)-1)  /theta(3); 
           end
           
           ind=find(delta==1);
           H(ind)=H(ind)+h(ind);
           LogLik=sum(H,'omitnan');
       end
       
   end
    %classdef ends
end