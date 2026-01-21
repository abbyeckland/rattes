%% sediment stored in terminal dams
function [calc1all_wrivers,calc6_delta] = DeltaTerm(calc1all_wrivers,calc6_delta,calc5all,data5all,data1wrivers,data6_delta)

calc5=calc5all;
data=data5all;
USdams=find(data.SID<500000 & data.PermStorag==1);

calc1=calc1all_wrivers;
calc6=calc6_delta;

deldata=data6_delta;

riverdata=data1wrivers;
rivers=find(data1wrivers.isriver==1);

numt=numel(calc5.t);
numdelt=numel(deldata.deltaID);
numall=numel(data1wrivers.SID);

%% get number of dams above rivers 

calc1.DamsAboveRivers=zeros(numall,numt);
calc1.TermSedAboveRivers=zeros(numall,numt);
calc1.TermDamsAboveRivers=zeros(numall,numt);
calc1.maxTermDamsAboveRivers=zeros(numall,numt);

%get the number of dams built above rivers through time
for j=rivers'
    riverSID=riverdata.SID(j);
    idx=[];
    termdams=[];
    for iA=USdams' % identify USdams that are above this river
        tmpA = data.rivertags{iA};
        if ismember([riverSID],[tmpA])==1
            idx= [idx iA]; %index location for a dam that flows to target river
        end
    end

    if isempty(idx)~=1
        idx=idx';
        termdams=find(data.term(idx)==1);
        calc1.maxTermDamsAboveRivers(j,1)=numel(termdams);
        yrbuilt=calc5.builtidx(idx);
        yrremoved=calc5.removedidx(idx);
        yrremoved = yrremoved(yrremoved > 0);
        for i=1:numt
            numbuilt=numel(find(yrbuilt<=i));
            numremoved=numel(find(yrremoved<=i));
            calc1.DamsAboveRivers(j,i)=numbuilt-numremoved;
            if isempty(termdams)~=1
                termidx=idx(termdams);
                yrbuiltterm=calc5.builtidx(termidx);
                yrremovedterm=calc5.removedidx(termidx);
                yrremovedterm = yrremovedterm(yrremovedterm > 0);
                numbuiltterm=numel(find(yrbuiltterm<=i));
                numremovedterm=numel(find(yrremovedterm<=i));
                calc1.TermDamsAboveRiversUS(j,i)=numbuiltterm-numremovedterm;
                calc1.TermSedAboveRivers(j,i)=sum(calc5.sed(termidx,i));
            end
            clear numbuilt numremoved termidx numbuiltterm numremovedterm
        end
    end
end
clear j iA i idx tmpA

%% get number of dams above deltas 

calc6.DamsAboveDeltas=zeros(numdelt,numt);
calc6.TermSedAboveDeltas=zeros(numdelt,numt);
calc6.maxTermDamsAboveDeltas=zeros(numdelt,1);
% calc6.MaxPathxTermDelta=zeros(numdelt,1);
% calc6.MinPathXTermDelta=zeros(numdelt,1);
% calc6.MedPathXTermDelta=zeros(numdelt,1);
% calc6.MeanPathXTermDelta=zeros(numdelt,1);

%get the number of dams built above deltas through time
for j=1:numdelt
    deltaID=deldata.deltaID(j);
    idx=[];
    termdams=[];
    termidx=[];
    for iA=USdams' % identify USdams that are above this river
        tmpA = data.deltatag(iA);
        if ismember([deltaID],[tmpA])==1
            idx= [idx iA]; %index location for a dam that flows to target delta
        end
    end

    if isempty(idx)~=1
        idx=idx';
        termdams=find(data.term(idx)==1);
        calc6.maxTermDamsAboveDeltas(j,1)=numel(termdams);
        yrbuilt=calc5.builtidx(idx);
        yrremoved=calc5.removedidx(idx);
        yrremoved = yrremoved(yrremoved > 0);
        if isempty(termdams)~=1
                termidx=idx(termdams);
                % calc6.MaxPathxTermDelta(j,1)=max(data.pathx(termidx));
                % calc6.MinPathXTermDelta(j,1)=min(data.pathx(termidx));
                % calc6.MedPathXTermDelta(j,1)=median(data.pathx(termidx));
                % calc6.MeanPathXTermDelta(j,1)=mean(data.pathx(termidx));
        end
        
        for i=1:numt
            numbuilt=numel(find(yrbuilt<=i));
            numremoved=numel(find(yrremoved<=i));
            calc6.DamsAboveDeltas(j,i)=numbuilt-numremoved;
            if isempty(termdams)~=1
                calc6.TermSedAboveDeltas(j,i)=sum(calc5.sed(termidx,i));
            end
            clear numbuilt numremoved
        end
    end
clear numbuilt numremoved 
end
clear j iA i idx tmpA

%now=find(calc5.t==2025);

calc6.PercentSedInTermDams=(calc6.TermSedAboveDeltas./calc6.sedabovedeltas)*100;
%calc6.Sed2025TermDeltaDams=calc6.TermSedAboveDeltas(:,now);
%calc6.PercentofTotalSedinTermDams2025=calc6.PercentSedInTermDams(:,now);

    
%% Save stuff

calc6_delta=calc6;
calc1all_wrivers=calc1;

save('deltaterm.mat','calc1all_wrivers','calc6_delta','-v7.3');
disp('*')
disp('Done with DeltaTerm function')
end