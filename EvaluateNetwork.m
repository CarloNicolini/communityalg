function Values=EvaluateNetwork(Res)

m=Res.m_orig;
uM=unique(m);
Values=nan(3,6,size(Res.pacoSurpriseWeight,1));
for i=1:3
    switch i
        case 1          
            s=Res.pacoSurpriseWeight;       
        case 2
            s=Res.pacoInfoMapWeight;
        case 3
            s=Res.NewmanWeight;
    end
    if not(size(s,2)==size(m,2))
        s=s';
        if not(size(s,2)==size(m,2))
            error('Problem with the size of your results')
        end
    end
    
    for itry=1:size(s,1)
        uS=unique(s(itry,:));
        finalMat=zeros(length(uS),length(uM));
        for iM=1:length(uM)
            mtmp=(m==uM(iM));
            for iS=1:length(uS)
                stmp=(s(itry,:)==uS(iS));
                finalMat(iS,iM)=sum(double(stmp.*mtmp));
            end
        end
        
        MccCommunity=zeros(1,length(uM));
        AccCommunity=zeros(1,length(uM));
        SensCommunity=zeros(1,length(uM));
        SpecCommunity=zeros(1,length(uM));
        PPVCommunity=zeros(1,length(uM));
        NPVCommunity=zeros(1,length(uM));
        
        for iM=1:length(uM)
            finalMattmp=finalMat;
            tab=zeros(1,4);
            [val,ind] =max(finalMattmp(:,iM));
            tab(1)=val;
            finalMattmp(ind,iM)=0;
            tab(2)=sum(finalMattmp(:,iM));
            tab(3)=sum(finalMattmp(ind,:));
            finalMattmp(ind,:)=0;finalMattmp(:,iM)=0;
            tab(4)=sum(finalMattmp(:));
            MccCommunity(iM)=MccVal(tab);
            AccCommunity(iM)=AccVal(tab);
            SensCommunity(iM)=SensitivityVal(tab);
            SpecCommunity(iM)=SpecificityVal(tab);
            PPVCommunity(iM)=PPVVal(tab);
            NPVCommunity(iM)=NPVVal(tab);
            
        end
        
        Values(i,1,itry)=mean(MccCommunity);
        Values(i,2,itry)=mean(AccCommunity);
        Values(i,3,itry)=mean(SensCommunity);
        Values(i,4,itry)=mean(SpecCommunity);
        Values(i,5,itry)=mean(PPVCommunity);
        Values(i,6,itry)=mean(NPVCommunity);
      %  Values(i,7:10,itry) =valid_RandIndex(s(itry,:),m+1);
    end
end



function coef=MccVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN

coef=(tab(1)*tab(4)+tab(2)*tab(3))/ ...
    sqrt((tab(1)+tab(2))*(tab(1)+tab(3))*(tab(4)+tab(2))*(tab(4)+tab(3)));

function coef=AccVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN
coef=(tab(1)+tab(4))/sum(tab)  ;

function coef=SensitivityVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN
coef=tab(1)/(tab(1)+tab(3));


function coef=SpecificityVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN
coef=tab(4)/(tab(4)+tab(2));


function coef=PPVVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN
coef=tab(1)/(tab(1)+tab(2));


function coef=NPVVal(tab)
%tab(1)=TP; tab(2)=FP; tab(3)=FN; tab(4)=TN
coef=tab(4)/(tab(4)+tab(3));