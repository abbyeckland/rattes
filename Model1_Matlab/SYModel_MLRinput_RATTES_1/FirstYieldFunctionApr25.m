function [data3SYmodel,calc3_SY1] = FirstYieldFunctionApr25(data3SYmodel,calc3_SY1,convert1,convert2)
data=data3SYmodel;
calc=calc3_SY1;

%% Trap efficiency

%trap efficiency at non-sites (static with time)
t=calc.t;
numt=numel(t);
numdam=height(data);
notsite = find(data.issite==0); 
site = find(data.issite == 1);

b=NaN(numdam,1);
m=NaN(numdam,1);
rise=NaN(numdam,1);
run=NaN(numdam,1);
AveTrap=NaN(numdam,1);

%Trap efficiency at sites, linear interp for now
rise(site)=calc.trap2(site)-calc.trap1(site);
run(site)=(data.yr2(site)-data.yr1(site));%
m(site)=rise(site)./run(site); %
b(site)=calc.trap1(site)-(m(site).*data.yr1(site));

for j=site'
       
        %fill in year for survey 1, capacity 
        calc.cap(j,calc.sur1idx(j))=data.cap1(j);

        %trap including and between between survey 1 and 2 (matches calc.trap1 & calc.trap2 exactly)
        btwn=(calc.sur1idx(j):1:calc.sur2idx(j));
        calc.trap(j,btwn)=(m(j)*(t(btwn)))+b(j); %%%straight line
        
        %fill in year for survey 2, capacity 
        calc.cap(j,calc.sur2idx(j))=data.cap2(j);

        %calc the ave trap efficiency over the period between surveys for sites. linear for now
        temp=calc.trap(j,calc.sur1idx(j):calc.sur2idx(j)-1);
        temp2=numel(temp);
        AveTrap(j,1)=(sum(temp))/(temp2); % this is the average trap effiency btwn surveys
        clear temp temp2
 end
clear j

% Non-sites, trap efficiency at non-sites w/ capacity static for now. 
for j=notsite'
    % find pre& post dam values in time loop
    if ~isnan(calc.removedidx(j))
        postdam=(calc.builtidx(j):(calc.removedidx(j)-1));
    else
        postdam=(calc.builtidx(j):numt);
    end
    %fill in static non site data, because this fills in maxstor post completion, it already covers the gap between yrc and yrp
    calc.trap(j,postdam)=1-1./(1+calc.kappa(j)*(((data.maxstor(j))/convert1)./((data.DA(j))/convert2))); %stationary non-site trap efficiency
    calc.cap(j,postdam)=data.maxstor(j);
    calc.AveTrap(j,1)=calc.trap(j,(postdam(1))); %since static just fills in the first trap efficiency
    clear postdam
      
end
clear j

%% Unravel DA changes and calc cap and trap

% this is the sediment contributing drainage area upstream from a reservoir X (does not include trap efficiency at reservoir X). 
%starts at orig value and decreases as upstream dams are built
calc.SAatdam=(calc.origDA);

imin=1;
imax=max(data.Rank);
for i=imin:imax
    ranknum=find(data.Rank==i);
    calc.sedshed(ranknum,:)=calc.SAatdam(ranknum,:).*calc.trap(ranknum,:); %km2, this will calc calc.sedshed except beyond second survey for the rank we are on
    calc.SedDAtoDS(ranknum,:)=calc.SAatdam(ranknum,:)-calc.sedshed(ranknum,:); %km2, this is the volume moving downstream past dam site in any given year, up to second survey.
    for val=ranknum'
        if data.issite(val)==1
            built=calc.builtidx(val);
            sur1=calc.sur1idx(val);
            sur2=calc.sur2idx(val);
            calc.wsedshed(val)=(sum(calc.sedshed(val,sur1:(sur2-1))))/(numel(calc.sedshed(val,sur1:(sur2-1)))); % km2, Time weighted calc.sedshed (effective sed contributing DA) between surveys, at ReservoirX
            calc.wSAatdam(val)=(sum(calc.SAatdam(val,sur1:(sur2-1))))/(numel(calc.SAatdam(val,sur1:(sur2-1)))); %km2, Time weighted sediment contributing DA upstream of Reservoir X
            calc.AveTrap(val)=calc.wsedshed(val)/calc.wSAatdam(val);%unitless, TE
            calc.wseddel(val)=(data.cap1(val)-data.cap2(val))/(calc.AveTrap(val)); %m3, volume of sediment delivered to reservoir X, between sur 1 and sur 2 (some passed through)
            calc.wSDR(val)=calc.wseddel(val)/(numel(sur1:sur2-1));%m3/yr, Sediment delivery rate, mean volume of sediment delivered to surveyed reservoir X per year between surveys
            calc.wSDRyield(val)=calc.wSDR(val)/calc.wSAatdam(val); %m3/km2/yr, sediment yield rate, volume of sed per year per km2 between sur 1 and sur 2. EXACT SAME AS reservoir sediment volume/((length sur1:sur2-1)*(calc.wsedshed));. in other works, volume of sediment delivered divided by Us sed DA is same as volume of sediment retained divided by effective sed DA
            
            %capacity sur1 to sur2
            nmin=sur1+1; nmax=sur2-1;
            for n=nmin:nmax %
                calc.cap(val,n)=calc.cap(val,n-1)-(calc.wSDRyield(val)*calc.sedshed(val,n-1)); %m3; capacity is former capacity minus sediment yeild times "effective Sed contributing DA"
            end
            clear n nmin nmax

            %beyond sur2 to end of time or dam removal
            if isnan(calc.removedidx(val))
                nmax=numt;
            else
                nmax=calc.removedidx(val)-1;
            end
            nmin=sur2+1; 
            for n=nmin:nmax
                calc.cap(val,n)=calc.cap(val,n-1)-(calc.wSDRyield(val)*calc.sedshed(val,n-1));
                                
                %don't let capacity go below 0
                if(calc.cap(val,n)<=0)
                    calc.cap(val,n)=0;
                    calc.trap(val,n)=0;
                    calc.sedshed(val,n)=0;
                    calc.SedDAtoDS(val,n)=calc.SAatdam(val,n);
                else
                    calc.trap(val,n)=1-1/(1+calc.kappa(val)*(((calc.cap(val,n))/convert1)/((data.DA(val))/convert2)));
                    calc.sedshed(val,n)=calc.SAatdam(val,n).*calc.trap(val,n);
                    calc.SedDAtoDS(val,n)=calc.SAatdam(val,n)-calc.sedshed(val,n);
                end
            end
            clear n nmin nmax;

            %Between yr comp and sur1. We will now set this to be a straight line for now and back project with final rate 
            if built<sur1
                nmin=built; nmax=sur1-1;
                for n=nmax:-1:nmin
                    calc.cap(val,n)=calc.cap(val,(n+1))+(calc.wSDRyield(val)*calc.sedshed(val,n+1)); %calcs capacity between completion and first survey. this would be the weighted DA at that timestep. (so tot DA-DA of upstream dams)
                    calc.trap(val,n)=1-1/(1+calc.kappa(val)*(((calc.cap(val,n))/convert1)/((data.DA(val))/convert2)));
                    calc.sedshed(val,n)=calc.SAatdam(val,n).*calc.trap(val,n);
                    calc.SedDAtoDS(val,n)=calc.SAatdam(val,n)-calc.sedshed(val,n);
                end
                clear n nmin nmax
            end
            % for sites we want the maxstor data to match the capacity at yrc 
            data.maxstor(val)=calc.cap(val,built);
            data.capp(val)=data.maxstor(val);
            clear sur1 sur2 built;
        end
        if data.term(val)==0 %if it isn't a term dam, move that drainage area downstream
            goesto=data.todam(val);
            thatdam=find(data.SID==goesto);
            calc.SAatdam(thatdam,:)=calc.SAatdam(thatdam,:)-(data.DA(val)-calc.SedDAtoDS(val,:)); %Adjusts SA at dam for next time loop. Rank of this dam must be higher than dam it comes from
            %double check for any drainage area errors, can be caused by bad dam snaps on flow diversions
            if any(calc.SAatdam(thatdam,:)<0)
                error('Drainage area errors')
            end
            clear thatdam goesto
        end
    end
end
    
%% Save data for output
calc3_SY1.AveTrap=AveTrap;
data3SYmodel=data;
calc3_SY1=calc;
calc3_SY1.wSDRyield1=calc3_SY1.wSDRyield;
calc3_SY1=rmfield(calc3_SY1, 'wSDRyield');

save('Yield1.mat','data3SYmodel','calc3_SY1','-v7.3');
disp('Done with First Yield Function')

end
