%% compare MLR & sed yield models, find dams that aren't sites but exist in both models
function [fignum,calc5all,data5all] = CompareSY_MLRAug25(calc5all,data5all,fignum)

%% input data 

%SY model only want to compare output from dams that were modeled using SY rates from permanent storage dams
data=data5all;
calc=calc5all;

numdam=numel(data.SID);
numt=numel(calc.t);

data.compare(:,1)=0;
data.compare(data.inSY == 1 & data.inMLR == 1, 1) = 1;

%% locate dams and clean data. put some same timeseries. Get design timeseries and surveyed timeseries

%locate dams with predictions that were made using sediment yield rates at
%dams upstream from sites that do not have permanent storage. Don't compare
%here.

nonsite=find(data.issite==0 & data.inSY==1);
for j=1:numel(nonsite)
    val=nonsite(j);
    tosite=data.sitetag1(val); %this is the immediate downstream site controlling the SY rate.
    toidx=find(data.SID==tosite);
    if data.PermStorag(toidx)==0
        data.compare(val,1)=0; %don't compare at these sites
    end
end

%Canada.  Don't compare at canada Don't count as "in MLR" or "inSY"
canada=find(data.SID>=500000);
data.compare(canada)=0;
data.inMLR(canada)=0;
data.inSY(canada)=0;

% create a capdesign file for MLR and for SY and for all. Start with the
% cap file b/c will fill in with Nans for dams not in MLR or SY model and
% then zeros for before dam built. Will overwrite. And will be all zeros if
% it was removed.
calc.capdesign=calc.cap;
calc.capdesignMLR=calc.capMLR;
calc.capdesignSY=calc.capSY;

for i=1:numdam
    %set the time for this to run until
    start=calc.builtidx(i);
    %get stop
    if ~isnan(calc.removedidx(i))
        stop=(calc.removedidx(i))-1;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    calc.capdesign(i,start:stop)=calc.cap(i,start);
    calc.capdesignMLR(i,start:stop)=calc.capMLR(i,start);
    calc.capdesignSY(i,start:stop)=calc.capSY(i,start);
    clear start stop
end
clear i start stop

calc.capsurveyed=calc.capdesign;
site=find(data.issite==1);
for i=site'
    %set the time for this to run until
    start1=calc.sur1idx(i);
    start2=calc.sur2idx(i);
    %get stop
    if ~isnan(calc.removedidx(i))
        stop=(calc.removedidx(i))-1;
    else
        stop=numt; %either the year before it is removed or the end of the time loop
    end
    calc.capsurveyed(i,start1:(start2-1))=calc.cap(i,start1);
    calc.capsurveyed(i,start2:stop)=calc.cap(i,start2);
end

    
%% Make figure
tplot=find(calc.t==1900);
now=find(calc.t==2025);
compare=find(data.compare==1);
c=numel(compare);

figure(fignum)
hold on
str1='Compare Model Results';
plot(calc.t(tplot:numt),(sum(calc.capdesignMLR(compare,tplot:numt))/1e9),'k')
plot(calc.t(tplot:numt),(sum(calc.capdesignSY(compare,tplot:numt))/1e9),'b')
plot(calc.t(tplot:numt),(sum(calc.capMLR(compare,tplot:numt))/1e9),'k:')
plot(calc.t(tplot:numt),(sum(calc.capSY(compare,tplot:numt))/1e9),'b--')
            ax = gca;
            ax.FontSize = 16;
            xlabel('Year')
            xlim([1900 2050])
            ax.YColor='k';
            ylabel('Total Capacity, km^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= append(str1,', n = ',num2str(c));
            title(caption)
            legend('DesignMLR','DesignSY','SSUR Model','SYSRModel','Location','southeast')
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\ModelCompare\'];

            saveas(h,fullfile(folder,filename),'pdf')
            hold off
            
fignum=fignum+1;
clear c h ax filename str1 caption str2 strcombo

%% plot 2025 capacity comparison log log
d=numel(compare);

xline=(0.1:10000000:1000000000000);
yline=(0.1:10000000:1000000000000);

figure(fignum)
str1='Compare Capacity SSUR Model to SYSR model 2025';
logplotcapSY=calc.capSY(compare,now);
replace=find(logplotcapSY<1);
logplotcapSY(replace)=0.1;
clear replace
logplotcapMLR=calc.capMLR(compare,now);
replace=find(logplotcapMLR<1);
logplotcapMLR(replace)=0.1;
a=find(logplotcapSY==0.1 & logplotcapMLR>0.1);
b=find(logplotcapMLR==0.1 & logplotcapSY>0.1);

MLRtotalcapwhereSYhasnocap2025_km3=sum(logplotcapMLR(a))/1e9;
SYtotalcapwhereMLRhasnocap2025_km3=sum(logplotcapSY(b))/1e9;
disp(['MLR total capacity where SYSR model has no capacity in 2025: ',num2str(MLRtotalcapwhereSYhasnocap2025_km3),' km3']);
disp(['SYSR total capacity where SSUR Model has no capacity in 2025: ',num2str(SYtotalcapwhereMLRhasnocap2025_km3),' km3']);

loglog(logplotcapSY(:),logplotcapMLR(:),'.')
hold on
loglog(xline,yline,'k','LineWidth',2);
grid on
  ax = gca;
            ax.FontSize = 16;
            xlabel('Capacity SYSR Model, m^3')
            ax.YColor='k';
            ylabel('Capacity SSUR Model, m^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= append(str1,', (n = ',num2str(d),')');
            title(caption)
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\ModelCompare\'];
            
            saveas(h,fullfile(folder,filename),'pdf')
            hold off
  fignum=fignum+1;
  clear h ax filename str1 caption logplotcapSY logplotcapMLR replace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xline=(0.1:1000000:2000000000);
yline=(0.1:1000000:2000000000);
  
figure(fignum)
str1='Compare sedimentation from SSUR Model to SYSR Model 2025';


MLRsed2025_km3=sum(calc.sedMLR(compare,now))/1e9;
SYsed2025_km3=sum(calc.sedSY(compare,now))/1e9;
disp(['MLR, total reservoir sedimentation 2025 at comparison sites: ',num2str(MLRsed2025_km3),' km3']);
disp(['SY, total reservoir sedimentation 2025 at comparison sites: ',num2str(SYsed2025_km3),' km3']);

loglog(calc.sedSY(compare,now),calc.sedMLR(compare,now),'.')
hold on
loglog(xline,yline,'k','LineWidth',2);
grid on
  ax = gca;
            ax.FontSize = 16;
            xlabel('Sediment Volume SYSR Model, m^3')
            ax.YColor='k';
            ylabel('Sediment Volume SSUR Model, m^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= append(str1,', (n = ',num2str(d),')');
            title(caption)
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\ModelCompare\'];
            
            saveas(h,fullfile(folder,filename),'pdf')
            hold off
  fignum=fignum+1;
  clear  h ax filename str1 caption replace
          
%% 2050
xline=(0.1:10000000:1000000000000);
yline=(0.1:10000000:1000000000000);

figure(fignum)
str1='Compare Capacity SSUR Model to SYSR Model 2050';
logplotcapSY=calc.capSY(compare,numt);
replace=find(logplotcapSY<1);
logplotcapSY(replace)=0.1;
clear replace
logplotcapMLR=calc.capMLR(compare,numt);
replace=find(logplotcapMLR<1);
logplotcapMLR(replace)=0.1;
a=find(logplotcapSY==0.1 & logplotcapMLR>0.1);
b=find(logplotcapMLR==0.1 & logplotcapSY>0.1);

MLRtotalcapwhereSYhasnocap2050_km3=sum(logplotcapMLR(a))/1e9;
SYtotalcapwhereMLRhasnocap2050_km3=sum(logplotcapSY(b))/1e9;
disp(['MLR total capacity where SYSR Model has no capacity in 2050  at comparison sites: ',num2str(MLRtotalcapwhereSYhasnocap2050_km3),' km3']);
disp(['SY total capacity where SSUR Model has no capacity in 2050 at comparison sites: ',num2str(SYtotalcapwhereMLRhasnocap2050_km3), ' km3']);

loglog(logplotcapSY(:),logplotcapMLR(:),'.')
hold on
loglog(xline,yline,'k','LineWidth',2);
grid on
  ax = gca;
            ax.FontSize = 16;
            xlabel('Capacity SYSR Model, m^3')
            ax.YColor='k';
            ylabel('Capacity SSUR Model, m^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= append(str1,', (n = ',num2str(d),')');
            title(caption)
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\ModelCompare\'];
            
            saveas(h,fullfile(folder,filename),'pdf')
            hold off
  fignum=fignum+1;
  clear h ax filename str1 caption logplotcapSY logplotcapMLR replace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xline=(0.1:1000000:2000000000);
yline=(0.1:1000000:2000000000);
  
figure(fignum)
str1='Compare sedimentation from SSUR Model to SYSR Model 2050';

MLRsed2050_km3=sum(calc.sedMLR(compare,numt))/1e9;
SYsed2050_km3=sum(calc.sedSY(compare,numt))/1e9;
disp(['MLR, total reservoir sedimentation 2050 at comparison sites: ',num2str(MLRsed2050_km3),' km3']);
disp(['SY, total reservoir sedimentation 2050 at comparison sites: ',num2str(SYsed2050_km3), ' km3']);

loglog(calc.sedSY(compare,numt),calc.sedMLR(compare,numt),'.')
hold on
loglog(xline,yline,'k','LineWidth',2);
grid on
  ax = gca;
            ax.FontSize = 16;
            xlabel('Sediment Volume SYSR Model, m^3')
            ax.YColor='k';
            ylabel('Sediment Volume SSUR Model, m^3')
            
            h=gcf; set(h,'color','w')
            fmt='%s.pdf';
            filename=sprintf(fmt,str1);
            caption= append(str1,', (n = ',num2str(d),')');
            title(caption)
            box on
            set(h,'Position',[50 50 1200 800]);
            set(h,'PaperOrientation','landscape');
            h.PaperPositionMode = 'manual';
            orient(h,'landscape')
                      
            folder=[pwd '\Figures\ModelCompare\'];
            
            saveas(h,fullfile(folder,filename),'pdf')
            hold off
  fignum=fignum+1;
  clear  h ax filename str1 caption replace

  %% save
data5all=data;
calc5all=calc;

close all

save('CompareModels.mat','data5all','calc5all','-v7.3');

disp('*')
disp('Done with Model Comparison')
end


