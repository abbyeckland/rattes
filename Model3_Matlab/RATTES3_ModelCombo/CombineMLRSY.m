% Load MLR data
function [data5all,calc5all] = CombineMLRSY(filenameMLRtrap,filenameMLRcap,filenameMLRcapup,filenameMLRcaplow,filenameMLRsv,filenameMLRsvup,filenameMLRsvlow,data2norivers,calc2_MLRinput,data4SYmodel,calc4_SY2,convert1,convert2)%

%% read files
% read files, does not have site data included

%trap, time, and SID
MLRtrap=readtable(filenameMLRtrap,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
trapMLR=table2array(MLRtrap);

%tMLR=str2double(MLRtrap.Properties.VariableNames);
SIDMLR=str2double(MLRtrap.Properties.RowNames);

%cap, with error
MLRcap=readtable(filenameMLRcap,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
capMLR=table2array(MLRcap);

MLRcapup=readtable(filenameMLRcapup,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
capMLRup=table2array(MLRcapup);

MLRcaplo=readtable(filenameMLRcaplow,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
capMLRlo=table2array(MLRcaplo);

%sv, with error
MLRsv=readtable(filenameMLRsv,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
svMLR=table2array(MLRsv);

MLRsvup=readtable(filenameMLRsvup,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
svMLRup=table2array(MLRsvup);

MLRsvlo=readtable(filenameMLRsvlow,'ReadRowNames', true,'ReadVariableNames',true,'VariableNamingRule','preserve');
svMLRlo=table2array(MLRsvlo);

clear MLRcap MLRcapup MLRcaplow MLRtrap MLRsv MLRsvup MLRsvlow 

%% define data and calc fields
data=data2norivers;
dataSY=data4SYmodel;
SIDall=data.SID;
SIDSY=dataSY.SID;

calc=calc2_MLRinput; % this is the dataset we will fill in to create a new combo
calcSY=calc4_SY2;
t=calc.t;

numdam=height(data);
numt=numel(t);

%identify dams in the MLR predictor list and flag identifier. Fo now, some dams will be missing until I get new input from abby.
data.inMLR=ismember(SIDall,SIDMLR);
data.inSY=ismember(SIDall,SIDSY); 

data.capcMLR=NaN(numdam,1);
data.capcMLRup=NaN(numdam,1);
data.capcMLRlo=NaN(numdam,1);
data.capcALL=NaN(numdam,1); %this will combine national model output for sites and MLR output for non-sites

% List of field names to initialize
fieldNames = {'trapMLR','capMLR', 'capMLRlo', 'capMLRup','sedMLR','sedMLRup','sedMLRlo','seddtMLR' ...
              'trapSY', 'capSY', 'sedSY','seddtSY'};

% Preallocate NaN arrays for all fields'
for i = 1:numel(fieldNames)
    calc.(fieldNames{i}) = NaN(numdam, numt);
end
clear i

%fields with zeros
calc.trapMLRup=zeros(numdam,numt);

calc.wSDRyieldSY=NaN(numdam,1);
calc.SAatdam_MLRinput=calc.SAatdam;
calc.sedshed_MLRinput=calc.sedshed;

% List of fields to remove
fieldsToRemove = {'trap', 'cap', 'AveTrap', 'SAatdam', 'SedDAtoDS', ...
                  'sedshed', 'wSDR', 'wSDRyield', 'wsedshed', 'wSAatdam', 'wseddel'};

% Remove the fields
calc = rmfield(calc, fieldsToRemove);

%sites and not sites
notsite=find(data.issite==0);
site=find(data.issite==1);

%% rename fields to combine, SIDs are all ascending, so can infill data

%create a subset of SY data that we will want in the final table.
dataSYsub=dataSY(:,{'SID','maxstor'});
dataSYsub=renamevars(dataSYsub,"maxstor","capcSY"); %rename to capc and infill with data.maxstor

%now infill calc datasets SY model
inSY=find(data.inSY==1);
calc.capSY(inSY,:)=calcSY.cap; %fill in from SY output
calc.trapSY(inSY,:)=calcSY.trap; %fill in from SY output
calc.sedSY(inSY,:)=calcSY.sed; %fill in from SY output
calc.seddtSY(inSY,:)=calcSY.sed_dt; %fill in from SY output
calc.wSDRyieldSY(inSY,:)=calcSY.wSDRyield2; %fill in from SY output
calc.SAatdamSY(inSY,:)=calcSY.SAatdam;
calc.sedshedSY(inSY,:)=calcSY.sedshed;

%now infill calc datasets MLR model
inMLR=find(data.inMLR==1); 
calc.capMLR(inMLR,:)=capMLR; %fill in from MLR output
calc.capMLRup(inMLR,:)=capMLRup; %fill in from MLR output
calc.capMLRlo(inMLR,:)=capMLRlo; %fill in from MLR output

calc.trapMLR(inMLR,:)=trapMLR; %fill in from MLR output
for i=inMLR'
    calc.trapMLRup(i,calc.predictidx(i))=calc.trapp(i);
end
calc.trapMLRlo=calc.trapMLRup;

calc.sedMLR(inMLR,:)=svMLR; %fill in from MLR output
calc.sedMLRup(inMLR,:)=svMLRup; %fill in from MLR output
calc.sedMLRlo(inMLR,:)=svMLRlo; %fill in from MLR output

%fill in begining of ALL dataset; this is the SY model at the sites and the
%MLR model everywhere else
calc.capALL=calc.capMLR; %create dataset 
calc.capALLup=calc.capMLRup;
calc.capALLlo=calc.capMLRlo;

calc.trapALL=calc.trapMLR; %create dataset
calc.trapALLup=calc.trapMLRup;
calc.trapALLlo=calc.trapMLRlo;

calc.sedALL=calc.sedMLR; %create dataset
calc.sedALLup=calc.sedMLRup; %create dataset
calc.sedALLlo=calc.sedMLRlo; %create dataset

%SY model is treated as a known in the model combination, with no error
calc.capALL(site,:)=calc.capSY(site,:); %fill in from other struct for sites
calc.capALLup(site,:)=calc.capSY(site,:); %fill in from other struct for sites
calc.capALLlo(site,:)=calc.capSY(site,:); %fill in from other struct for sites

%SY model is treated as a known in the model combination, with no error
calc.trapALL(site,:)=calc.trapSY(site,:); %fill in from other struct for sites
calc.trapALLup(site,:)=calc.trapSY(site,:); %fill in from other struct for sites
calc.trapALLlo(site,:)=calc.trapSY(site,:); %fill in from other struct for sites

%SY model is treated as a known in the model combination, with no error
calc.sedALL(site,:)=calc.sedSY(site,:); %fill in from other struct for sites
calc.sedALLup(site,:)=calc.sedSY(site,:); %fill in from other struct for sites
calc.sedALLlo(site,:)=calc.sedSY(site,:); %fill in from other struct for sites

%create sed dt timeseries for MLR model and then create ALL dataset for
%seddt
calc.seddtMLR(inMLR,:)=0;
calc.seddtMLRup(inMLR,:)=calc.seddtMLR(inMLR,:);
calc.seddtMLRlo(inMLR,:)=calc.seddtMLR(inMLR,:);

%get sed_dt timeseries for the MLR model, fill in low and high trap from
%predictidx + 1 to end of time or removal
for i= inMLR'
    %set the time for this to run until
    if ~isnan(calc.removedidx(i))
        stop=(calc.removedidx(i))-1;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    start=calc.predictidx(i);
    start2=start+1;

    calc.trapMLRup(i,start)=1-1./(1+calc.kappa(i)*(((calc.capMLRup(i,start))/convert1)./((data.DA(i))/convert2))); %trap efficiency low and high between predict year and end
    calc.trapMLRlo(i,start)=1-1./(1+calc.kappa(i)*(((calc.capMLRlo(i,start))/convert1)./((data.DA(i))/convert2))); %trap efficiency low and high between predict year and end
    
    for n=start2:stop
        calc.seddtMLR(i,n)=calc.capMLR(i,n-1)-calc.capMLR(i,n);
        calc.seddtMLRup(i,n)=calc.capMLRlo(i,n-1)-calc.capMLRlo(i,n); %low capacity means more sediment, so sedhi goes with caplo
        calc.seddtMLRlo(i,n)=calc.capMLRup(i,n-1)-calc.capMLRup(i,n); % high capacity means less sediment, so sedlo goes with capup
        
        % calc trapup and trapdown for MLR model
        calc.trapMLRup(i,n)=1-1./(1+calc.kappa(i)*(((calc.capMLRup(i,n))/convert1)./((data.DA(i))/convert2))); %trap efficiency low and high between predict year and end
        calc.trapMLRlo(i,n)=1-1./(1+calc.kappa(i)*(((calc.capMLRlo(i,n))/convert1)./((data.DA(i))/convert2))); %trap efficiency low and high between predict year and end

    end
    clear start; clear stop; 
end
clear i

calc.seddtALL=calc.seddtMLR;
calc.seddtALLup=calc.seddtMLRup;
calc.seddtALLlo=calc.seddtMLRlo;

%%%% For the national sedimentation model, sites are treated as a known (no error estimate).
calc.seddtALL(site,:)=calc.seddtSY(site,:);
calc.seddtALLup(site,:)=calc.seddtSY(site,:);
calc.seddtALLlo(site,:)=calc.seddtSY(site,:);

%% Find RESNET reservoirs without perm storage MLR model didn't predict at.  keep storage constant
% find dams in old MLR output that have been deleted from RESNET
noperm=find(data.PermStorag==0 & data.issite==0);

%We just want to keep the storage static from built to removed.
for i=noperm'
    start=calc.builtidx(i);
     %set the time for this to run until
    if ~isnan(calc.removedidx(i))
        stop=(calc.removedidx(i))-1;

        %fill in zeros for capacity and trap after removed
        calc.capALL(i,(stop+1):numt)=0;
        calc.capALLlo(i,(stop+1):numt)=0;
        calc.capALLhi(i,(stop+1):numt)=0;

        calc.trapALL(i,(stop+1):numt)=0;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    
    %infill from beginning of time to built
    calc.capALL(i,1:(start-1))=0;
    calc.capALLup(i,1:(start-1))=0;
    calc.capALLlo(i,1:(start-1))=0;
    
    calc.trapALL(i,1:(start-1))=0;
    
    %will never have sed b/c holding cap steady. from from 1:numt sed=0
    calc.sedALL(i,:)=0;
    calc.sedALLup(i,:)=0;
    calc.sedALLlo(i,:)=0;
    
    calc.seddtALL(i,:)=0;
    calc.seddtALLup(i,:)=0;
    calc.seddtALLlo(i,:)=0;
    
    %fill in from built to removed or end of time.
    calc.capALL(i,start:stop)=data.maxstor(i);
    calc.capALLup(i,start:stop)=data.maxstor(i);
    calc.capALLlo(i,start:stop)=data.maxstor(i);
    
    calc.trapALL(i,start:stop)=calc.trapp(i);
    calc.trapALLup(i,start:stop)=calc.trapp(i);
    calc.trapALLlo(i,start:stop)=calc.trapp(i);
end
clear i

%% BackFill MLR data to get total sed volume and sed_dt.  
%The MLR predictions are based on the "design Sed contributing DA." Get sed. yield for Design Sed contribute DA and project mean back through time. 

goback=find((data.yrc<data.yrp) & data.inMLR==1);

for idx=goback'
    tstart = (find(calc.seddtMLR(idx, :) > 0, 1)); % Find the first index with sedimentation
    tend = find(calc.seddtMLR(idx, :) > 0, 1, 'last'); % Find the last index with sedimentation
    
    %model result
    SYtemp=(calc.seddtMLR(idx,tstart:tend))./(calc.SAatdam_MLRinput(idx,((tstart-1):(tend-1))).*calc.trapMLR(idx,((tstart-1):(tend-1)))); % because the first year sed yield predicts the second year sedimentation.  diff is X(2)-X(1)
    SYmean=(sum(SYtemp))/(numel(tstart:tend));
    
    %up err
    SYtempup=(calc.seddtMLRup(idx,tstart:tend))./(calc.SAatdam_MLRinput(idx,((tstart-1):(tend-1))).*calc.trapMLR(idx,((tstart-1):(tend-1)))); % because the first year sed yield predicts the second year sedimentation.  diff is X(2)-X(1)
    SYmeanup=(sum(SYtempup))/(numel(tstart:tend));
    
    %lo err
    SYtemplo=(calc.seddtMLRlo(idx,tstart:tend))./(calc.SAatdam_MLRinput(idx,((tstart-1):(tend-1))).*calc.trapMLR(idx,((tstart-1):(tend-1)))); % because the first year sed yield predicts the second year sedimentation.  diff is X(2)-X(1)
    SYmeanlo=(sum(SYtemplo))/(numel(tstart:tend));
    
    nmax=calc.predictidx(idx,1)-1;
    nmin=calc.builtidx(idx,1);

    for n = (nmax:-1:nmin)
        calc.capMLR(idx,n)=calc.capMLR(idx,(n+1))+(SYmean*(calc.SAatdam_MLRinput(idx,(n+1)).*calc.trapMLR(idx,(n+1)))); % This is based on the design SA used in the MLR model. 
        calc.capMLRup(idx,n)=calc.capMLRup(idx,(n+1))+(SYmeanup*(calc.SAatdam_MLRinput(idx,(n+1)).*calc.trapMLRup(idx,(n+1)))); %for projecting backward, hi sediment yield goes with hi capacity  
        calc.capMLRlo(idx,n)=calc.capMLRlo(idx,(n+1))+(SYmeanlo*(calc.SAatdam_MLRinput(idx,(n+1)).*calc.trapMLRlo(idx,(n+1)))); %  high sediment yield goes wtih lo capacity
                    
        calc.trapMLR(idx,n)=1-1/(1+(calc.kappa(idx))*(((calc.capMLR(idx,n))/convert1)/((data.DA(idx))/convert2)));
        calc.trapMLRup(idx,n)=1-1/(1+(calc.kappa(idx))*(((calc.capMLRup(idx,n))/convert1)/((data.DA(idx))/convert2)));
        calc.trapMLRlo(idx,n)=1-1/(1+(calc.kappa(idx))*(((calc.capMLRlo(idx,n))/convert1)/((data.DA(idx))/convert2)));      
    end

    data.maxstor(idx)=calc.capMLR(idx,nmin);
    maxstoruptemp=calc.capMLRup(idx,nmin);
    maxstorlotemp=calc.capMLRlo(idx,nmin);

    %update total sed volume until dam removal or end of time
     %set the time for this to run until
    if ~isnan(calc.removedidx(idx))
        stop=(calc.removedidx(idx))-1;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    start=nmin+1;%this is yrc + 1
    pause=nmax+1; %this is yrp
    
    %update calc.sedMLR using sedimentation since construction
    calc.sedMLR(idx,start:stop)=data.maxstor(idx)-calc.capMLR(idx,start:stop);
    calc.sedMLRup(idx,start:pause)=maxstoruptemp-calc.capMLRup(idx,start:pause); % between yrc and yrp, hi capacity back projection accumulates sediment faster to intersect at yrp
    calc.sedMLRup(idx,(pause+1):stop)=maxstoruptemp-calc.capMLRlo(idx,(pause+1):stop); % between yrp and end, the low capacity estimate is accumulating sediment faster
    
    calc.sedMLRlo(idx,start:pause)=maxstorlotemp-calc.capMLRlo(idx,start:pause); % between yrc and yrp, low capacity back projection accumulates sediment slower to intersect at yrp
    calc.sedMLRlo(idx,(pause+1):stop)=maxstorlotemp-calc.capMLRup(idx,(pause+1):stop); % between yrp and end, high capacity back projection accumulates sediment slower 
    

    for n=start:stop
        calc.seddtMLR(idx,n)=calc.sedMLR(idx,n)-calc.sedMLR(idx,n-1);
        calc.seddtMLRup(idx,n)=calc.sedMLRup(idx,n)-calc.sedMLRup(idx,n-1); % hi sed low cap
        calc.seddtMLRlo(idx,n)=calc.sedMLRlo(idx,n)-calc.sedMLRlo(idx,n-1); %low sed high cap
    end
    clear SYmean SYtemp tstart tend nmax nmin maxstoruptemp maxstorlotemp
end
clear idx

calc.sedALL(goback,:)=calc.sedMLR(goback,:);
calc.sedALLup(goback,:)=calc.sedMLRup(goback,:);
calc.sedALLlo(goback,:)=calc.sedMLRlo(goback,:);

calc.seddtALL(goback,:)=calc.seddtMLR(goback,:);
calc.seddtALLup(goback,:)=calc.seddtMLRup(goback,:);
calc.seddtALLlo(goback,:)=calc.seddtMLRlo(goback,:);

calc.capALL(goback,:)=calc.capMLR(goback,:);
calc.capALLup(goback,:)=calc.capMLRup(goback,:);
calc.capALLlo(goback,:)=calc.capMLRlo(goback,:);

calc.trapALL(goback,:)=calc.trapMLR(goback,:);
calc.trapALLup(goback,:)=calc.trapMLRup(goback,:);
calc.trapALLlo(goback,:)=calc.trapMLRlo(goback,:);

%% combine fields from SY dataset into here
data=outerjoin(data,dataSYsub,'Keys','SID','MergeKeys',true,'Type','left');
data.capcMLR(notsite,1)=data.maxstor(notsite,1);
data.capcALL=data.capcMLR;
data.capcALL(site,1)=data.capcSY(site,1); %overwrites sites with SY capc for those sites

data.maxstor=[];

%% save mat file and output
%decided to get rid of the ALL and have the unspecified variables just be cap trap sed and seddt b/c easier.  kept above for clarity
calc.cap=calc.capALL;
calc.capup=calc.capALLup;
calc.caplo=calc.capALLlo;

calc.trap=calc.trapALL;
calc.trapup=calc.trapALLup;
calc.traplo=calc.trapALLlo;

calc.sed=calc.sedALL;
calc.sedup=calc.sedALLup;
calc.sedlo=calc.sedALLlo;

calc.seddt=calc.seddtALL;
calc.seddtup=calc.seddtALLup;
calc.seddtlo=calc.seddtALLlo;

% List of fields to remove
fieldsToRemove = {'trapALL','trapALLup','trapALLlo','capALL','capALLup','capALLlo','sedALL','sedALLup','sedALLlo','seddtALL','seddtALLup','seddtALLlo'};
% Remove the fields
calc = rmfield(calc, fieldsToRemove);

data5all=data;
calc5all=calc;

save('ComboModel.mat','data5all','calc5all','-v7.3');
disp('*')
disp('Done with Model Combo')


end







