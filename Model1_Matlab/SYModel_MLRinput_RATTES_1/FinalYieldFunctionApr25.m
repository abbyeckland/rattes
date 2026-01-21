%function [data4SYmodel,calc4_SY2] = FinalYieldFunctionApr25(data3SYmodel,calc4_SY2,calc3_SY1,convert1,convert2,OutputCapName,OutputSedName, OutputTrapName)% 
function [data4SYmodel,calc4_SY2] = FinalYieldFunctionApr25(data3SYmodel,calc4_SY2,calc3_SY1,convert1,convert2)% 

%% pre allocate
data=data3SYmodel; 
calc=calc4_SY2;

% setup
t=calc.t;

numdam=numel(data.SID);
numt=numel(t);

calc.SAatdam=calc.origDA;
site=find(data.issite==1);
notsite=find(data.issite==0);
calc.wSDRyield2=NaN(numdam,1);
calc=rmfield(calc,'wSDRyield');


%head=find(data.head==1);
calc.trapc=NaN(numdam,1);


%% fill in start data; make a field of all dams that flow to each site, find the min yrc for all dams upstream from a site. fill in cap and trap1
siteidx=cell(numdam,1);

for j=site'
    idx=[];
    idx=find(data.sitetag1==data.SID(j));
    siteidx{j,1}=[idx];
    if data.yrc(j)<data.yr1(j)
        data.maxstor(j)=NaN;
    end
    calc.trap(j,(calc.sur1idx(j)))=calc.trap1(j);
    calc.cap(j,(calc.sur1idx(j)))=data.cap1(j);
end

%% fill in between yrc and yrp for dams with only 1 survey that postdates the completion year, but aren't sites. Treated as static for now but later back projected with solution.
for i=notsite'    
        calc.trapc(i)=1-1./(1+calc.kappa(i)*(((data.maxstor(i))/convert1)./((data.DA(i))/convert2)));
        calc.trap(i,calc.builtidx(i))=calc.trapc(i);
        calc.cap(i,calc.builtidx(i))=data.maxstor(i);
        if (data.yrp(i)~=data.yrc(i)) % sets data constant for non-sites before dam closure.  will back calculate later.
            calc.trap(i,(calc.builtidx(i)):(calc.predictidx(i)))=calc.trapc(i);
            calc.cap(i,(calc.builtidx(i)):(calc.predictidx(i)))=data.maxstor(i); %keeps constant trap efficiency and capacity between year built and prediction year. For when 1st survey after built or survey showed capacity growth.
        end
end

calc.sedshed=calc.trap.*calc.SAatdam;
calc.SedDAtoDS=calc.SAatdam-calc.sedshed;

clear i; clear idx; clear j; clear jmax; clear a; clear idxtemp; clear siteSID; clear val; clear head;

%% Go through by rank, find the final sediment yield rate

iterationfinal=NaN(numdam,1);
leng=3;
imax=max(data.Rank);


for i=1:imax %looping through Ranks from low to high.
    tic
    disp(['Loop on Rank i: ', num2str(i)])
    val=find(data.Rank==i & data.issite==1);
    for j=val'
        ready=0;
        iteration=0;
        hi=calc3_SY1.wSDRyield1(j)*2; % high value
        lo=calc3_SY1.wSDRyield1(j)*.05; % low value
        wSDRyield=calc3_SY1.wSDRyield1(j); %starting SDR yield test, taken from the site
        testvalue = wSDRyield;

        %identify dams that flow to this site
        idx=siteidx{j,1};

        %if wSDRyield is already finalized it has already been removed from SA at the site. So it is already accounted for.
        %happens when there is a site above a site.
        idx(~isnan(calc.wSDRyield2(idx))) = [];

        % now we need to go from start to survey 2 to update SDA at dam for the target dam.
        while ready==0
            %moved this here
            SAatdamtemp=calc.SAatdam;
            capcalctemp=calc.cap;
            calctraptemp=calc.trap;
            SedDAtoDStemp=calc.SedDAtoDS;
            sedshedtemp=calc.sedshed;
            iteration=iteration+1;
            if iteration>100
                error('Lots of iterations')
            end

            % flush out SA at dam with upstream non-site dams
            if ~isempty(idx)%
                %order the upstream dams by rank so we start at low rank and then move to high rank, so sediment moves downstream
                sortme = [data.Rank(idx), idx];
                D = sortrows(sortme);  % Sort rows based on the first column
                idx = D(:, 2);         % Extract sorted indices
                for a=idx'  %
                    nmin=calc.predictidx(a)+1;
                    for n=nmin:1:numt %this fills in sediment to the end of the timeloop
                        if t(n)<data.yrr(a)% but will only fill sediment if dam has not been removed
                            capcalctemp(a,n)=capcalctemp(a,n-1)-(testvalue*sedshedtemp(a,n-1));
                            %if sedimentation made exceeds available capacity set capacity at 0
                            if(capcalctemp(a,n)<=0) % If there isn't any capacity there isn't any trap efficiency.
                                capcalctemp(a,n)=0;
                                calctraptemp(a,n)=0;
                                sedshedtemp(a,n)=0;
                                SedDAtoDStemp(a,n)=SAatdamtemp(a,n);
                            else % if there is capacity, update TE and sedshed and sedDADS.
                                calctraptemp(a,n)=1-1/(1+calc.kappa(a)*(((capcalctemp(a,n))/convert1)/((data.DA(a))/convert2)));
                                sedshedtemp(a,n)=SAatdamtemp(a,n).*calctraptemp(a,n);
                                SedDAtoDStemp(a,n)=SAatdamtemp(a,n)-sedshedtemp(a,n);
                            end
                        else %dam has been removed
                            SedDAtoDStemp(a,n)=SAatdamtemp(a,n);
                        end
                    end
                    %Move that drainage area downstream !
                    goesto=data.todam(a); %this is the dam SID that Dam a goes to (could be another dam between Dam a and target site)
                    thatdam=find(data.SID==goesto); %this is the index location of that dam
                    SAatdamtemp(thatdam,:)=SAatdamtemp(thatdam,:)-(calc.origDA(a,:)-SedDAtoDStemp(a,:)); %Adjusts SA at dam for next time loop. Rank of this dam must be higher than dam it comes from
                    sedshedtemp(thatdam,:)=SAatdamtemp(thatdam,:).*calctraptemp(thatdam,:); %update sedshed to reflect changes in SAatdam
                    SedDAtoDStemp(thatdam,:)=SAatdamtemp(thatdam,:)-sedshedtemp(thatdam,:);

                    %check for errors
                    if any(SAatdamtemp(thatdam,:)<0)
                        error('Drainage area errors downstream from idx(k)')
                    end
                    clear thatdam; clear goesto;
                end
            end

            % OKAY now that SDA at dam is set up, test testvalue at the site
            nmin=calc.sur1idx(j)+1; nmax=calc.sur2idx(j); %
            sedshedtemp(j,:)=SAatdamtemp(j,:).*calctraptemp(j,:); %

            for n=nmin:nmax
                capcalctemp(j,n)=capcalctemp(j,n-1)-(testvalue*sedshedtemp(j,n-1));
                if(capcalctemp(j,n)<=0)
                    capcalctemp(j,n)=0;
                    calctraptemp(j,n)=0;
                    sedshedtemp(j,n)=0;
                    SedDAtoDStemp(j,n)=SAatdamtemp(j,n);
                else
                    calctraptemp(j,n)=1-1/(1+calc.kappa(j)*(((capcalctemp(j,n))/convert1)/((data.DA(j))/convert2)));
                    sedshedtemp(j,n)=SAatdamtemp(j,n)*calctraptemp(j,n);
                    SedDAtoDStemp(j,n)=SAatdamtemp(j,n)-sedshedtemp(j,n);
                end
            end

            if capcalctemp(j,nmax)>(data.cap2(j)*.9999) && capcalctemp(j,nmax)<(data.cap2(j)*1.0001) % means you found a solution
            %if capcalctemp(j,nmax)>(data.cap2(j)*.98) && capcalctemp(j,nmax)<(data.cap2(j)*1.02) %faster loop for troubleshooting
                ready=1;
                calc.wSDRyield2(j)=testvalue;
                calc.wSDRyield2(idx)=testvalue;
                iterationfinal(j)=iteration;
                iterationfinal(idx)=iteration;

                %let temp overwrite the final set, because it worked
                calc.SAatdam=SAatdamtemp;
                calc.cap=capcalctemp;
                calc.trap=calctraptemp;
                calc.SedDAtoDS=SedDAtoDStemp;
                calc.sedshed=sedshedtemp;
                clear sedshedtemp; clear SedDAtoDStemp; clear calctraptemp; clear capcalctemp; clear SAatdamtemp;
            else
                %reset temp datasets to try again
                if capcalctemp(j,nmax)>data.cap2(j)
                    lo=testvalue;
                else
                    hi=testvalue;
                end
                if lo==hi
                    error('Lo equals Hi')
                end
                clear sedshedtemp; clear SedDAtoDStemp; clear calctraptemp; clear capcalctemp; clear SAatdamtemp;
                testrange=linspace(lo,hi,leng);
                testvalue=median(testrange);
            end
        end

        clear idx; clear sortme; clear D;

        %Once wsDR is solved for, need to fill in sedimentation at site from survey 2 to end of time loop, if it has permanent storage
        nmin=calc.sur2idx(j)+1; nmax=numt;
        for n=nmin:nmax
            if data.PermStorag(j)==1 % if the site has permanent storage
                if t(n)<data.yrr(j)
                    calc.cap(j,n)=calc.cap(j,n-1)-(calc.wSDRyield2(j)*calc.sedshed(j,n-1));
                end
                %don't let capacity go below 0
                if(calc.cap(j,n)<=0)
                    calc.cap(j,n)=0;
                    calc.trap(j,n)=0;
                    calc.sedshed(j,n)=0;
                    calc.SedDAtoDS(j,n)=calc.SAatdam(j,n);
                else
                    calc.trap(j,n)=1-1/(1+calc.kappa(j)*(((calc.cap(j,n))/convert1)/((data.DA(j))/convert2)));
                    calc.sedshed(j,n)=calc.SAatdam(j,n).*calc.trap(j,n);
                    calc.SedDAtoDS(j,n)=calc.SAatdam(j,n)-calc.sedshed(j,n);
                end
            else % so if there isn't permanent storage leave capacity after survey 2 constant
                if t(n)<data.yrr(j)
                    calc.cap(j,n)=calc.cap(j,(calc.sur2idx(j)));
                end
                if(calc.cap(j,n)<=0)
                    calc.cap(j,n)=0;
                    calc.trap(j,n)=0;
                    calc.sedshed(j,n)=0;
                    calc.SedDAtoDS(j,n)=calc.SAatdam(j,n);
                else
                    calc.trap(j,n)=1-1/(1+calc.kappa(j)*(((calc.cap(j,n))/convert1)/((data.DA(j))/convert2)));
                    calc.sedshed(j,n)=calc.SAatdam(j,n).*calc.trap(j,n);
                    calc.SedDAtoDS(j,n)=calc.SAatdam(j,n)-calc.sedshed(j,n);
                end
            end
        end
        clear n; clear nmin; clear nmax;

        %now fill in from yrcomp to yr1 Between yr comp and sur1.
        if data.yrc(j)<data.yr1(j)
            nmin=calc.builtidx(j); nmax=calc.sur1idx(j)-1;
            for n=nmax:-1:nmin
                if data.PermStorag(j)==1 %only backproject storage at reservoirs with permanent storage
                    calc.cap(j,n)=calc.cap(j,(n+1))+(calc.wSDRyield2(j)*calc.sedshed(j,n+1)); %calcs capacity between completion and first survey. this would be the weighted DA at that timestep. (so DA-DA of upstream dams)- so if there were fewer dams the prod would be bigger.
                else
                    calc.cap(j,n)=calc.cap(j,(calc.sur1idx(j)));
                end
                calc.trap(j,n)=1-1/(1+calc.kappa(j)*(((calc.cap(j,n))/convert1)/((data.DA(j))/convert2)));
                calc.sedshed(j,n)=calc.SAatdam(j,n).*calc.trap(j,n);
                calc.SedDAtoDS(j,n)=calc.SAatdam(j,n)-calc.sedshed(j,n);
            end
            clear n nmin nmax;
        end
        %now if this isn't a terminal dam, move that drainage area downstream
        if data.term(j)==0 %if it isn't a term dam, move that drainage area downstream
            goesto=data.todam(j);
            thatdam=find(data.SID==goesto);
            calc.SAatdam(thatdam,:)=calc.SAatdam(thatdam,:)-(calc.origDA(j,:)-calc.SedDAtoDS(j,:)); %Adjusts SA at dam for next time loop. Rank of this dam must be higher than dam it comes from
            calc.sedshed(thatdam,:)=calc.SAatdam(thatdam,:).*calc.trap(thatdam,:); %update sedshed to reflect changes in SAatdam
            calc.SedDAtoDS(thatdam,:)=calc.SAatdam(thatdam,:)-calc.sedshed(thatdam,:);

            %double check for any drainage area errors, can be caused by flow diversions
            if any(calc.SAatdam(thatdam,:)<0)
                error('Drainage area errors downstream from target site')
            end
            clear goesto; clear thatdam; 
        end
    end
    toc
end

%% backproject sedimentation at sites where yrp> yrc using the wSDR yield 
%sediment contributing drainage area will be off at these sites but the total sedimentation will be more accurate.
%this will also yield a better comparison with the MLR model, which does back project. 
% not updating SAatdam or sedshed

predictdam=find(data.issite==0 & data.PermStorag==1 &(data.yrp>data.yrc));
for j=predictdam'
    nmin=calc.builtidx(j); nmax=calc.predictidx(j)-1;
    for n=nmax:-1:nmin
        calc.cap(j,n)=calc.cap(j,(n+1))+(calc.wSDRyield2(j)*calc.sedshed(j,n+1)); %calcs capacity between completion and first survey. this would be the weighted DA at that timestep. (so DA-DA of upstream dams)- so if there were fewer dams the prod would be bigger.
        calc.trap(j,n)=1-1/(1+calc.kappa(j)*(((calc.cap(j,n))/convert1)/((data.DA(j))/convert2)));
    end
    calc.trapc(j)=calc.trap(j,nmin);
    data.maxstor(j)=calc.cap(j,nmin);
    clear n nmin nmax;
end



%% for sites we want the maxstor data to match the capacity at yrc 

for i=site'
    data.maxstor(i)=max(calc.cap(i,:));
end
data.capp(site)=data.maxstor(site);
calc.trapc(site)=1-1./(1+calc.kappa(site).*(((data.maxstor(site))/convert1)./((data.DA(site))/convert2)));        
clear i; 

%calc sed timeseries.
calcsed=zeros(numdam,numt);
calcsed_dt=zeros(numdam,numt);
for i= 1:numdam
    %set the time for this to run until
    if ~isnan(calc.removedidx(i))
        stop=(calc.removedidx(i))-1;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    start=calc.builtidx(i)+1;
    calcsed(i,start:stop)=data.maxstor(i)-calc.cap(i,start:stop);
    for n=start:stop
        calcsed_dt(i,n)=calc.cap(i,n-1)-calc.cap(i,n);
    end
    clear start stop; 
end


%% Create output save mat
calc4_SY2=calc;
calc4_SY2.sed=calcsed;
calc4_SY2.sed_dt=calcsed_dt;
calc4_SY2.iterations=iterationfinal;

data4SYmodel=data; 

% rows=numdam+1;
% columns=numt+1;
% SedOutput_Compare_wMLR=NaN(rows,columns);
% SedOutput_Compare_wMLR(2:rows,1)=data.SID;
% SedOutput_Compare_wMLR(1,2:columns)=calc.t(1,:);
% CapOutput_Compare_wMLR =SedOutput_Compare_wMLR;
% TrapOutput_Compare_wMLR=SedOutput_Compare_wMLR;
% SedOutput_Compare_wMLR(2:rows,2:columns)=calcsed;
% CapOutput_Compare_wMLR(2:end,2:end)=calc.cap;
% TrapOutput_Compare_wMLR(2:end,2:end)=calc.trap;
% 
% writematrix(CapOutput_Compare_wMLR,OutputCapName,'Delimiter',',');
% writematrix(SedOutput_Compare_wMLR,OutputSedName,'Delimiter',',');
% writematrix(TrapOutput_Compare_wMLR,OutputTrapName,'Delimiter',',');

save('YieldFinal.mat','data4SYmodel','calc4_SY2','-v7.3');
disp('Done with Final Yield Function')


end