function [calc1all_wrivers,data1wrivers,data6_delta,calc6_delta,calc5all,fignum] = SDArivers(calc1all_wrivers,data1wrivers,calc5all,data5all,fignum)%

data=data1wrivers;
calc=calc1all_wrivers;
rivers=find(data.isriver==1);
notriver=find(data.isriver==0);
notriverUS=find(data.isriver==0 & data.SID<500000 & data.PermStorag==1);
numall=numel(data.SID);
numt=numel(calc.t);

%years for sea level rate changes
years = [1900 1970 1971 2005 2006 2010 2011 2020 2021 2030 2031 2040 2041 2050];
idxt = arrayfun(@(y) find(calc.t == y), years);

%% Calculate Sediment contributing drainage area at river mouths (keep Canada for DA)

calc.sedDAtoDS=NaN(numall,numt); % the is the DA that moves downstream, start with NaN. will be overwritten
calc.SAatdam=(calc.origDA);% this is the sediment contributing drainage area upstream from a reservoir X (does not include trap efficiency at reservoir X).
calc.sedshed=NaN(numall,numt);

%re-do the ranking loop with rivers included- treat them like 0 trap
%efficiency dams
imin=1;
imax=max(data.Rank);
for i=imin:imax
    ranknum=find(data.Rank==i);
    calc.sedshed(ranknum,:)=calc.SAatdam(ranknum,:).*calc.trap(ranknum,:); %km2, this will calc calc.sedshed except beyond second survey for the rank we are on
    calc.SedDAtoDS(ranknum,:)=calc.SAatdam(ranknum,:)-calc.sedshed(ranknum,:); %km2, this is the volume moving downstream past dam site in any given year, up to second survey.
    for val=ranknum'
        if data.term(val)==0 %if it isn't a term dam, move that drainage area downstream
            goesto=data.todam(val);
            thatdam=find(data.SID==goesto);
            calc.SAatdam(thatdam,:)=calc.SAatdam(thatdam,:)-(data.DA(val)-calc.SedDAtoDS(val,:)); %Adjusts SA at dam for next time loop. Rank of this dam must be higher than dam it comes from
            %double check for any drainage area errors, can be caused by
            %flow diversions
            if any(calc.SAatdam(thatdam,:)<0)
                error('Drainage area errors downstream')
            end
            clear thatdam goesto
        end
    end
end


%% With national sediment model compute sedSed contributing DA and record
calc5all.SAatdamALL=calc.SAatdam(notriver,:);
calc5all.sedshedALL=calc.sedshed(notriver,:);
calc5all.SedDAtoDSALL=calc.SedDAtoDS(notriver,:);

% get sed yield for entire model
now=find(calc.t==2025);
numdam=numel(notriver);
calc5all.YieldAlldt=zeros(numdam,numt);
calc5all.wYield25=zeros(numdam,1);
calc5all.YieldAlldt(:,2:end)=(calc5all.seddt(:,2:end))./(calc5all.SAatdamALL(:,(1:(end-1)))).*calc5all.trap(:,(1:(end-1))); % because the first year trapping and capacity yield predicts the next year sedimentation.

for i=1:numdam
    start=calc5all.builtidx(i);
    if ~isnan(calc5all.removedidx(i))
        stop=(calc5all.removedidx(i))-1;
    else
        stop=now; %either the year before it is removed or the end of the time loop
    end
    calc5all.wsedshed25(i)=(sum(calc5all.sedshedALL(i,start:(stop-1))))/(numel(calc5all.sedshedALL(i,start:(stop-1))));
    calc5all.wYield25(i)=(calc5all.sed(i,stop))/(numel(calc5all.sedshedALL(i,start:(stop-1)))*calc5all.wsedshed25(i));
end

%% Basin trapping above deltas
% this is the percent of the upstream basin that is trapped 
calc.PercentUpstreamBasinTrapping=zeros(numall,numt);
calc.PercentUpstreamBasinTrapping=(calc.origDA-calc.SAatdam)./(calc.origDA)*100;

%% calc delta sediment contributing drainage area above deltas (Keep Canada for DA)

filename= 'DeltaIDs.csv'; %all dams & rivermouths
Deltas=readtable(filename,'Format','auto');
clear filename;
Deltas.Var3=[]; %doesn't show on csv but imports empty
numdelta=numel(Deltas.deltaID);
deltaSAatdam=zeros(numdelta,numt);
deltaorigDA=zeros(numdelta,numt);

%a couple of non-terminal rivers were flagged as deltas in resnet. remove the delta
%IDs from those

a=find(data.deltaID>0 & data.rivertag1~=0);
data.deltaID(a)=0;
clear a;

for i=1:numdelta
    riverin=find(data.deltaID==i); % terminal river with a delta tag
    if isempty(riverin)~=1
        deltaorigDA(i,1)=sum(calc.origDA(riverin,1));
        for j=1:numt
            deltaSAatdam(i,j)=sum(calc.SAatdam(riverin,j));
        end
    end
end
%fill out origDA
deltaorigDA(:,2:numt) = repmat(deltaorigDA(:,1), 1, numt-1);
%Percent of upstream basin trapped by dams above deltas
deltasPercentUpstreamBasinTrapping=(deltaorigDA-deltaSAatdam)./deltaorigDA*100;

%% plot basin trapping above deltas (Keep Canada for DA)

maindeltas=[16,23,25,2,22,12,17,1,18,28,30,19,21,8,4,15];
maindeltas=maindeltas';

newcolors = [0 0 0
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560 
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330
             0.6350 0.0780 0.1840
             1 0 1
             0 1 1 
             0.5 0.5 0.5
             0.25 0.15 0.05  
             0.0000 0.4471 0.7412   
             0.90 0.80 0.95   
             0.8667 0.7961 0.6431   
             0.0000 0.8549 0.6902   
             0.1294 0.3412 0.1961];
         
colororder(newcolors)

str1='Total Basin Trapping by Reservoirs Above Deltas';
figure(fignum)
colororder(newcolors)
plot(calc.t(idxt(1):end),(deltasPercentUpstreamBasinTrapping(maindeltas,idxt(1):end)),'LineWidth',2);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Basin Trapping (%)')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= str1;
            title(caption)
            legend('Mississippi','Colorado','Columbia','Savannah','Rio Grande','Mobile','Trinity','Santee','Brazos','Elwha','Sabine','Colorado (Texas)','Nueces','Apalachicola','Altamaha','Pearl',...
                'Location','southeastoutside')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\DeltaFigs\'];
            saveas(h,fullfile(folder,filename),'pdf')
            
fignum=fignum+1;
clear c; clear h; clear ax; clear filename; clear str1; clear caption;



%% Calculate cumulative sediment above river mouths, and sed rate per year from US dams. 
% rivers, sediment upstream and number of dams Must have permanent storage
calc.sedaboverivers=zeros(numall,numt);
calc.seddtaboverivers=zeros(numall,numt);


for j=rivers'
    riverSID=data.SID(j);
    idx=[];
    for iA=notriverUS' % identify dams that are above this river
        tmpA = data.rivertags{iA};
        if ismember([riverSID],[tmpA])==1
            idx= [idx iA]; %index location for a dam that flows to target river
        end
    end

    if isempty(idx)~=1
        idx=idx';
        for i=1:numt
            calc.sedaboverivers(j,i)=sum(calc.sed(idx,i));
            calc.seddtaboverivers(j,i)=sum(calc.seddt(idx,i));
        end
    end
end

clear j iA i idx tmpA
%% calculate cumulative sediment above deltas, and sed rate per year; excludes canada

sedabovedeltas=zeros(numdelta,numt);
seddtabovedeltas=zeros(numdelta,numt);

for j=1:numdelta
    deltaID=Deltas.deltaID(j);
    idx=[];
    for iA=notriverUS' % identify dams that are above this delta
        tmpA = data.deltatag(iA);
        if ismember([deltaID],[tmpA])==1
            idx= [idx iA]; %index location for a dam that flows to target river
        end
    end

    if isempty(idx)~=1
        idx=idx';
        for i=1:numt
            sedabovedeltas(j,i)=sum(calc.sed(idx,i));
            seddtabovedeltas(j,i)=sum(calc.seddt(idx,i));
        end
    end
end

%% Equivalent depth above deltas (excludes Canada b/c those numbers not in sed and seddt)
%this is the depth of sediment that could be placed on a delta's entire
%area
deltaareas_m2=zeros(numdelta,numt);
deltaareas_m2(:,1)=Deltas.Area*1e6;
deltaareas_m2(:,2:numt) = repmat(deltaareas_m2(:,1), 1, numt-1);

deltatrap=0.3; %30% delta trapping efficiency

EquivDepth_Total=(sedabovedeltas*deltatrap)./deltaareas_m2; %meters
EquivDepth_dt=(seddtabovedeltas*deltatrap)./deltaareas_m2; %meters

%% plot delta figs
str1='Reservoir sedimentation rates above deltas';
figure(fignum)
colororder(newcolors)
semilogy(calc.t(idxt(1):end),seddtabovedeltas(maindeltas,idxt(1):end),'LineWidth',2);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Sedimentation Rate, m^3yr^-^1')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= str1;
            title(caption)
             legend('Mississippi','Colorado','Columbia','Savannah','Rio Grande','Mobile','Trinity','Santee','Brazos','Elwha','Sabine','Colorado (Texas)','Nueces','Apalachicola','Altamaha','Pearl',...
                 'Location','southoutside')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\DeltaFigs\'];
            saveas(h,fullfile(folder,filename),'pdf')
            
fignum=fignum+1;
clear h ax filename str1 caption

str1='Total reservoir sediment volume above deltas';
figure(fignum)
colororder(newcolors)
semilogy(calc.t(idxt(1):end),sedabovedeltas(maindeltas,idxt(1):end)/1e9,'LineWidth',2);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Sediment Volume, km^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= str1;
            title(caption)
            legend('Mississippi','Colorado','Columbia','Savannah','Rio Grande','Mobile','Trinity','Santee','Brazos','Elwha','Sabine','Colorado (Texas)','Nueces','Apalachicola','Altamaha','Pearl',...
                'Location','southoutside')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\DeltaFigs\'];
            saveas(h,fullfile(folder,filename),'pdf')
            
fignum=fignum+1;
clear h ax filename str1 caption

%% plot equivalent depths

%sea level rate Fox-Kemper et al....IPCC 6th assessment, 

slr=zeros(1,numt);
slr(1,idxt(1):idxt(2))=1.35; % mm/yr (palmer et al., 2021, cited in Fox-Kemper)
slr(1,idxt(3):idxt(4))=2.33; % mm/yr (tide gages, https://psmsl.org/ , cited in Fox-Kemper)
slr(1,idxt(5):idxt(6))=3.69; % mm/yr (Frederiske et al 2020, cited in Fox-Kemper)
slr(1,idxt(7):idxt(8))=4.9; % mm/yr (SSP2-4.5, 50th quantile,Fox-Kemper)
slr(1,idxt(9):idxt(10))=4.4; % mm/yr (SSP2-4.5, 50th quantile,Fox-Kemper)
slr(1,idxt(11):idxt(12))=4.9; % mm/yr (SSP2-4.5, 50th quantile,Fox-Kemper)
slr(1,idxt(13):idxt(14))=6.2; % mm/yr (SSP2-4.5, 50th quantile,Fox-Kemper)

x = calc.t(idxt(1):end).';

%% plot not main deltas

notmain=[3,5,6,7,9,10,11,13,14,20,24,26,31,32,33,34,35];
notmain=notmain';

%exluded  = 29, Dungeness (no dams); 27, Quillayute (no dams); 36,
%Sixes (no dams)

newcolors2 = [0 0 0
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560 
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330
             0.6350 0.0780 0.1840
             1 0 1
             0 1 1 
             0.5 0.5 0.5
             0.25 0.15 0.05  
             0.0000 0.4471 0.7412   
             0.90 0.80 0.95   
             0.8667 0.7961 0.6431   
             0.0000 0.8549 0.6902   
             0.1294 0.3412 0.1961
             0.6667 1.0000 0.7647];
         

            
%%
str1='Total equivalent delta aggradation depth deposited in upstream reservoirs';
figure(fignum)
colororder(newcolors)
semilogy(calc.t(idxt(1):end),EquivDepth_Total(maindeltas,idxt(1):end),'LineWidth',2);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Potential Aggradation Depth, m')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= str1;
            title(caption)
            legend('Mississippi','Colorado','Columbia','Savannah','Rio Grande','Mobile','Trinity','Santee','Brazos','Elwha','Sabine','Colorado (Texas)','Nueces','Apalachicola','Altamaha','Pearl','Location','southoutside')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\DeltaFigs\'];
            saveas(h,fullfile(folder,filename),'pdf')
            
fignum=fignum+1;
clear h ax filename str1 caption

%%
str1='Total equivalent delta aggradation depth deposited in upstream reservoirs, supplemental';
figure(fignum)
colororder(newcolors2)
semilogy(calc.t(idxt(1):end),EquivDepth_Total(notmain,idxt(1):end),'LineWidth',2);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Potential Aggradation Depth, m')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= str1;
            title(caption)
            legend('Ogee Chee','Satilla','St Marys','Ochlo','Choctawhatchee',...
                'Blackwater','Dead','W Pascagoula','Wolf','Lavaca','Eel',...
                'Humptulips','Suwannee','Myakka','Peace','Withlacoochee','Keys','Location','southoutside')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\DeltaFigs\'];
            saveas(h,fullfile(folder,filename),'pdf')
            
fignum=fignum+1;
clear h ax filename str1 caption


%% save things
data6_delta=Deltas;
calc6_delta.maindeltas=maindeltas;
calc6_delta.notmaindeltas=notmain;
calc6_delta.SAatdam=deltaSAatdam;
calc6_delta.PercentUpstreamBasinTrapping=deltasPercentUpstreamBasinTrapping;
calc6_delta.origDA=deltaorigDA;
calc6_delta.EquivDepth_dt=EquivDepth_dt;
calc6_delta.EquivDepth_Total=EquivDepth_Total;
calc6_delta.sealevelrise=slr;
calc6_delta.sedabovedeltas=sedabovedeltas;
calc6_delta.seddtabovedeltas=seddtabovedeltas;

%% Added this to limit output in calc1
% List of fields to remove
fieldsToRemove = {'trap1','trap2','kappa','builtidx','removedidx','sur1idx','sur2idx','predictidx','trapp','trapMLR','capMLR','capMLRlo','capMLRup','sedMLR','sedMLRup','sedMLRlo',...
    'seddtMLR','trapSY','capSY','sedSY','seddtSY','trapMLRup','wSDRyieldSY','SAatdam_MLRinput','sedshed_MLRinput','SAatdamSY','sedshedSY','trapMLRlo','seddtMLRup','seddtMLRlo','capdesignSY','capsurveyed'};

% Remove the fields
calc = rmfield(calc, fieldsToRemove);

data1wrivers=data;
calc1all_wrivers=calc;

%%
close all

save('riverSDA.mat','calc1all_wrivers','data1wrivers','data6_delta','calc6_delta','calc5all','-v7.3');
disp('*')
disp('Done with SDArivers function')

end