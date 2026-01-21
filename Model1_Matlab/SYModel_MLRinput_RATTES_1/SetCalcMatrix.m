function [calc2_MLRinput,calc3_SY1,calc4_SY2] = SetCalcMatrix(data2norivers,convert1,convert2)
data=data2norivers;

%% Create variables and empty datasets

%create the 1700ish to 2025 dataset use same timeset as MLR model
%tmin=(min(data.yrc))-1;
tmin=1699;
tmax=2050;
t=(tmin:1:tmax); clear tmin;
numt=numel(t);
numdam=height(data);

capcalc=NaN(numdam,numt); %this is a variable that would need to be stored in the structure ; I want Nans so I can tell if there was an error
sedshed=NaN(numdam,numt); %this is the km2 area of the watershed that has sediment getting trapped in reservoir (so contribut DA * trap efficiency)
SAatdam=NaN(numdam,numt); %this is the km2 area of the watershed that has sediment getting trapped in reservoir (so contribut DA * trap efficiency)
calctrap=zeros(numdam,numt);
calctrap1=NaN(numdam,1);
calctrap2=NaN(numdam,1);
wsedshed=NaN(numdam,1);% Time weighted (between surveys) SAatDam above Reservoir X * calctrap at dam X for Reservoir X ("effective sediment contributing DA")
wSAatdam=NaN(numdam,1); %Time weighted (between surveys) sediment-contributing drainage area above Reservoir X
AveTrap=NaN(numdam,1); %time weighted trap efficiency (between surveys)- used to constrain sediment delivery
wseddel=NaN(numdam,1); %m3, total volume of sediment delivered to reservoir X between survey 1 and survey 2
wSDR=NaN(numdam,1);%m3/yr, Sediment delivery rate, mean volume of sediment delivered to reservoir X per year between surveys
wSDRyield=NaN(numdam,1); %m3/(km3*t). sediment yeild. volume sed per km2 per yr
origDA = repmat(data.DA, 1, numt);
builtidx=NaN(numdam,1);
predictidx=NaN(numdam,1);
removedidx=NaN(numdam,1);
sur1idx=NaN(numdam,1);
sur2idx=NaN(numdam,1);
sedDAtoDS=NaN(numdam,numt); % the is the DA that moves downstream, start with NaN. will be overwritten


%% Trap efficiency

%notsite = data.issite == 0;
site = find(data.issite == 1);
%numsite = numel(site);
kappa = 0.1 * ones(numdam, 1); % Design assumption (silt: 0.1); coarse(sand) = 1, medium(silt)= 0.1, fine(clay)= 0.046;

%Trap efficiency at sites, linear interp for now
calctrap1(site,1)=1-1./(1+kappa(site).*(((data.cap1(site))/convert1)./((data.DA(site))/convert2)));%trap at 1st survey year
calctrap2(site,1)=1-1./(1+kappa(site).*(((data.cap2(site))/convert1)./((data.DA(site))/convert2))); %trap at 2nd survey year

% Find Indices for when dam was built and removed, fill in zeros. make linear trap at sites
for j=1:numdam
    %find pre& post dam indices and fill in 0s before dams built
    predictidx(j)=find(t==data.yrp(j));
    builtidx(j)=find(t==data.yrc(j));
    predam=(1:(builtidx(j)-1));
    calctrap(j,predam)=0;
    capcalc(j,predam)=0;%predam, capacity and trap efficiency = 0
    junk=find(t==data.yrr(j));
    
    %if dam was removed infill 0s after removal and created removed idx
    if ~isempty(junk)
        removedidx(j)=junk;
        removed=(removedidx(j):numt);
        
        % make cap and trap 0 after dam removal
        calctrap(j,removed)=0;
        capcalc(j,removed)=0; %post dam removal
    end
    %create indices for first and second survey years
    if data.issite(j)==1
        sur1idx(j)=find(t==data.yr1(j));
        sur2idx(j)=find(t==data.yr2(j));
    end
end
   

%% save variables
calc2_MLRinput.SID=data.SID;
calc2_MLRinput.trap=calctrap;
calc2_MLRinput.cap=capcalc;
calc2_MLRinput.SAatdam=SAatdam;
calc2_MLRinput.SedDAtoDS=sedDAtoDS;
calc2_MLRinput.sedshed=sedshed;
calc2_MLRinput.t=t;
calc2_MLRinput.AveTrap=AveTrap;
calc2_MLRinput.kappa=kappa;
calc2_MLRinput.trap1=(calctrap1);
calc2_MLRinput.trap2=(calctrap2);
calc2_MLRinput.wSDR=wSDR;
calc2_MLRinput.wSDRyield=wSDRyield;
calc2_MLRinput.builtidx=builtidx;
calc2_MLRinput.removedidx=removedidx;
calc2_MLRinput.sur1idx=sur1idx;
calc2_MLRinput.sur2idx=sur2idx;
calc2_MLRinput.origDA=origDA;
calc2_MLRinput.wsedshed=wsedshed;
calc2_MLRinput.wSAatdam=wSAatdam;
calc2_MLRinput.wseddel=wseddel;
calc2_MLRinput.predictidx=predictidx;

%% remove dams to nowwhere and dams without tags to match data 3
isolated=find(data.head == 1 & data.term == 1 & data.issite == 0); % Remove them from the dataset
notag=find(data.sitetag1 == 0 & data.issite == 0);
remove=vertcat(isolated,notag);

% Remove specified rows from calc2_MLRinput fields
fields = fieldnames(calc2_MLRinput);
calc3_SY1 = struct();

for i = 1:length(fields)
    field_data = calc2_MLRinput.(fields{i});
    if size(field_data, 1) == size(data.SID, 1) % Ensure row count matches
        calc3_SY1.(fields{i}) = field_data;
        calc3_SY1.(fields{i})(remove, :) = []; % Remove rows
    else
        calc3_SY1.(fields{i}) = field_data; % Keep unchanged for other dimensions
    end
end

calc4_SY2=calc3_SY1; %sets up the input for second yield also

%% save matrices
% save data3_SYmodel and mat

save('setcalc.mat','calc2_MLRinput','calc3_SY1','calc4_SY2','-v7.3');
disp('*');
disp('Done with SetCalcMatrix Function');


end




