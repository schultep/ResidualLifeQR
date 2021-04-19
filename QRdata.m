classdef QRdata
    properties
        id
        X
        X0
        T
        Y
        delta
        Z
        Y_MLE
        delta_MLE
        Z_MLE
        b_MLE
        b_QR
        b0_MLE
        b0_QR
        typeB
        wt
        trueTheta
        trueB
        remZ
        remXT
        remb0
        remb
    end
    
    methods
        function QRobj=QRdata(id,X,T,Y,delta,Z,typeB,theta,b)
            if nargin>=8
                QRobj.trueTheta=theta;
            end
            if nargin==9
                QRobj.trueB=b;
            end
            if nargin<7
                typeB=0;
            end
            QRobj.typeB=typeB;
            id=id-id(1)+1;
            csID=diff(id)-1;
            csID(csID<0)=0;
            csID=cumsum(csID);
            id(2:end)=id(2:end)-csID;
            
            QRobj.id=id; %NN=length(unique(id));
            QRobj.X=X;  %lnX=length(X);
            QRobj.T=T;
            QRobj.Y=Y;
            QRobj.delta=delta;
            QRobj.Z=Z;
            
            [~,begs]=unique(QRobj.id,'first');ends=[(begs(2:end)-1);length(T)];
            [f,t]=ecdf(Y(begs),'function','survivor',...
                'censoring',delta(begs),'bounds','on');
            if f(end)==0
                f(end)=f(end-1);  
            end
            QRobj.wt= zeros(size(Y));
            for i = 1:length(QRobj.wt)
                QRobj.wt(i)= (sum(id==id(i))-1)*f( sum(t<=Y(i))); 
            end
            QRobj.wt=QRobj.delta.*QRobj.wt.^-1;
            QRobj.wt(begs)=0;
            QRobj.Y_MLE=Y(begs);
            QRobj.delta_MLE=delta(begs);
            QRobj.Z_MLE=Z(begs,:);
            indB=find(ends>begs);
            QRobj.X0=X(begs);
            
            if typeB==0
                QRobj.b_MLE=0*X(begs);
                QRobj.b0_MLE=X(begs)';
                QRobj.b_MLE(indB)=(X(ends(indB))-X(begs(indB)))./(T(ends(indB))- T(begs(indB)));
                QRobj.b_MLE=QRdata.windsor(QRobj.b_MLE,1)';     
                QRobj.b_MLE(QRobj.b_MLE==0)=mean(QRobj.b_MLE(indB));
                b_QR1=(X- X(begs(id)))./(T- T(begs(id)));           
                b_QR1(begs)=0;
                b_QR2=diff(X)./diff(T);             
                b_QR2=[0;b_QR2];         
                b_QR2(begs)=0;
                QRobj.b_QR=(b_QR1+b_QR2)./2;
                QRobj.b0_QR=X;
            end
            if typeB==1               
                QRobj.b_QR=zeros([size(Y,1),1]);
                QRobj.b0_QR=zeros([size(Y,1),1]);
                for i=1:size(indB,1)              
                    g=indB(i);              
                    for j=(begs(g)+1):(ends(g))                    
                        md1=fitlm(T(begs(g):j) ,X(begs(g):j));
                        tX0=md1.Coefficients.Estimate(1,1);
                        tBX=md1.Coefficients.Estimate(2,1); 
                        QRobj.b_QR(j)=tBX;
                        QRobj.b0_QR(j)=tX0+tBX*T(j);                 
                    end                    
                end                           
                QRobj.b_QR=QRdata.windsor(QRobj.b_QR,.5);                
                QRobj.b_MLE=QRobj.b_QR(ends)';           
                QRobj.b0_MLE=QRobj.b0_QR(ends)';
            end
        end
    end
    methods (Static)
           
        function QRobj=removeInd(DAT,i)
            QRobj=DAT;
            ind=find(DAT.id==i);
            id=QRobj.id;
            id(ind)=[];
            id=id-id(1)+1;
            csID=diff(id)-1;
            csID(csID<0)=0;
            csID=cumsum(csID);
            id(2:end)=id(2:end)-csID;
            QRobj.id=id;
            QRobj.remXT=[QRobj.X(ind),QRobj.T(ind)];
            QRobj.remZ=QRobj.Z(ind,:);
            QRobj.remb=QRobj.b_QR(ind);
            QRobj.remb0=QRobj.b0_QR(ind);
            QRobj.X(ind)=[];
            QRobj.T(ind)=[];
            QRobj.Y(ind)=[];
            QRobj.delta(ind)=[];
            QRobj.Z(ind,:)=[];
            QRobj.b_QR(ind)=[];
            QRobj.b0_QR(ind)=[];
            QRobj.wt(ind)=[];
            
            QRobj.X0(i)=[];
            QRobj.Y_MLE(i)=[];
            QRobj.delta_MLE(i)=[];
            QRobj.Z_MLE(i,:)=[];
            QRobj.b_MLE(i)=[];
            QRobj.b0_MLE(i)=[];
        end

        function QRobj2=Boot(QRobj)
            NID=max(QRobj.id);
            ind=randi(NID,[NID,1]);
            idTmp=[]; whLine=[];
          
            for j=1:NID
                ind2= find(QRobj.id==ind(j)); 
                idTmp=[idTmp;j*ones([length(ind2),1])];
                whLine=[whLine;ind2];                
            end
           QRobj2=QRdata(idTmp,QRobj.X(whLine),QRobj.T(whLine),QRobj.Y(whLine),QRobj.delta(whLine),QRobj.Z(whLine,:));
        end

        function Y=windsor(Y,p)
            low=prctile(Y,p);
            high=prctile(Y,100-p);
            Y(Y<low)=low;
            Y(Y>high)=high;
        end
    end
end
