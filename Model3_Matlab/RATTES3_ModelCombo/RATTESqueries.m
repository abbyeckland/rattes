function [] = RATTESqueries(data5all,calc5all,data6_delta,calc6_delta)
data=data5all;

% Queries for Eckland et al.
data=data5all;
calc=calc5all;

calcdelta=calc6_delta;
datadelta=data6_delta;

smvol=100000000;% break between small / mid-size reservoirs and large reservoirs
canSID=500000; % all SIDs in canada are >500000

now=find(calc.t==2025);

%% Define masks
mask_rattes = (data5all.inMLR == 1 | (data5all.issite == 1 & data5all.PermStorag == 1)); %rattes is us only w/ perm storage (inMLR does not include canada)
mask_notremoved= (data.yrr > 2050); %already confirmed no future dam removal dates past 2025
mask_sm= (data.capcALL<smvol);
mask_lg= (data.capcALL>=smvol);
mask_site=(data.issite==1); %sites are only in us, not canada
mask_perm=(data.PermStorag==1);
mask_noperm=(data.PermStorag==0);
mask_usa=(data.SID<canSID);
mask_mlr=(data.inMLR==1);
mask_far=(data.pathx>=500);
mask_term=(data.term==1);
mask_head=(data.head==1);
mask_yrc=(data.yrcsub==0);
mask_grand=(data.isgrand==1);
mask_SY=(data.inSY==1 & data.PermStorag==1);
temp_owner=strcmp(data.OwnerTypes,'Private');
mask_private=(temp_owner==1);

%regional masks
mask_gulf=(data.countryOut==5);
mask_usatl=(data.countryOut==41);
mask_canatl=(data.countryOut==11);
mask_int=(data.countryOut==0);
mask_gl=(data.countryOut==2);
mask_mexpac=(data.countryOut==32);
mask_uspac=(data.countryOut==42);

%% Queries about total reservoir numbers in RATTES (not specific to 2025)
%Total Reservoirs with sedimentation predictions and their design capacity
rattes=find(mask_rattes & mask_usa);
siteperm=find(mask_perm & mask_site);
sitenoperm=find(mask_noperm & mask_site);
nopermtotal=find(mask_noperm & mask_usa); %no permanent storage excluding canada
inmlr=find(mask_mlr);
headsite=find(mask_perm & mask_site & mask_head);
storc=sum(data.capcALL(rattes))/1e9;
storcsite=sum(data.capcALL(siteperm))/1e9;
mediansitestor=median(data.capcALL(siteperm))/1e6;
medianstorall=median(data.capcALL(rattes))/1e6;
headresnet=find(mask_usa & mask_head);

%small/med reservoirs
smsitesidxall=(find(mask_site & mask_sm & mask_perm));
smsitegrand=(find(mask_site & mask_sm & mask_perm & mask_grand));

disp('*');
disp('*');
disp('*');
disp('****Queries about reservoir numbers in RESNET****');
disp(['Number of headwater, reservoirs in RESNET: ', num2str(numel(headresnet))]);
disp(['Number of reservoirs in RESNET: ', num2str(numel(find(mask_usa)))]);

disp('*');
disp('*');
disp('****Queries about reservoir numbers in RATTES****');
disp(['Number of US Reservoirs with predictions RATTES: ', num2str(numel(rattes))]);
disp(['Number of surveyed US Reservoirs with permanent storage in RATTES: ', num2str(numel(siteperm))]);
disp(['Unsurveyed reservoirs with permanent storage in RATTES: ', num2str(numel(inmlr))]);
disp(['Number of surveyed US Reservoirs without permanent storage in RATTES: ', num2str(numel(sitenoperm))]);
disp(['Total US Reservoirs without permanent storage in RATTES: ', num2str(numel(nopermtotal))]);
disp(['Number of headwater, surveyed reservoirs in RATTES: ', num2str(numel(headsite))]);
disp(['Total Design Capacity in RATTES: ', num2str(storc),' km3']);
disp(['Percent of all reservoirs surveyed : ',num2str(numel(siteperm)/numel(rattes)*100),'%']);
disp(['Percent of designed capacity held in surveyed sites : ',num2str(storcsite/storc*100),'%']);
disp(['Median Design Capacity in Surveyed Reservoirs: ', num2str(mediansitestor),' M m3']);
disp(['Median Design Capacity in RATTES: ', num2str(medianstorall),' M m3']);
disp(['Number of small/ medium surveyed reservoirs: ', num2str(numel(smsitesidxall))]);
disp(['Number of small/ medium surveyed reservoirs not included in GRAND: ', num2str(numel(smsitesidxall)-numel(smsitegrand))]);

disp('*');
disp('****Queries about SY model in RATTES****');
disp(['Median Sediment Yield Rate for Surveyed Sites, SY model result: ', num2str(median(calc.wSDRyieldSY(siteperm))),'m3 per km2 per yr']);


%% 2025 Results Queries of counts and general data  
%Find reservoirs that have not been removed by 2025 for consistent
%comparisons through time. 

%define reservoirs in 2025
rattes25=find(mask_rattes & mask_notremoved); %reservoirs in study that still exist in 2050 (check and same as reservoirs that still exist in 2025)
rattes25yrc=find(mask_rattes & mask_notremoved & mask_yrc);
site25=find(mask_site & mask_perm & mask_notremoved);
private25idx=find(mask_private & mask_rattes & mask_notremoved);

smidx25=(find(mask_rattes & mask_notremoved & mask_sm));
lgidx25=(find(mask_rattes & mask_notremoved & mask_lg));

%query data
sed25=sum(calc.sed(rattes25,now))/1e9; %sum of sediment volume in 2025, km3
sed_smmed25=sum(calc.sed(smidx25,now))/1e9; %sum of sediment volume in 2025, km3
sed_lg25=sum(calc.sed(lgidx25,now))/1e9;
storc25=sum(data.capcALL(rattes25))/1e9; % design storage for reservoirs that still exist, km3
caploss25=sed25/storc25*100; % percent capacity loss
sedmass=((sed25*1e9)*960)/(1e12); % dendy and chamption / minear and kondolf reservoir sediment mass
times_globload=sedmass/8.5; %8.5 gigatonnes is syvitski et al. 2022 annual fluvial sediment load
sed25rattesm3=(calc.sed(rattes25,now));
over1km3=find(calc.sed(rattes25,now)>1e9);
over1kmidx=rattes25(over1km3);

%percent loss
loss25rattes=sed25rattesm3./(data.capcALL(rattes25))*100;
tmp=find(loss25rattes>=25);
tmp2=find(loss25rattes>=50);
over25percentidx=rattes25(tmp);
over50percentidx=rattes25(tmp2);

%age
medyrc_rates25=median(data.yrc(rattes25yrc));

%sed rate
seddt25=sum(calc.seddt(rattes25,now))/1e9;
seddtlg25=sum(calc.seddt(lgidx25,now))/1e9;
seddtsm25=sum(calc.seddt(smidx25,now))/1e9;

% designed water storage large and small
storclg25=sum(data.capcALL(lgidx25))/1e9; % design storage for large reservoirs that still exist, km3
storcsm25=sum(data.capcALL(smidx25))/1e9; % design storage for small / med reservoirs that still exist, km3

disp('*');
disp('*');
disp('****QUERIES ABOUT RATTES 2025 RESULTS****');
disp('****Counts and general stats, 2025****');
disp(['Number of Remaining Reservoirs, 2025: ', num2str(numel(rattes25))]);
disp(['Number of Reservoirs with over 25% loss, 2025: ', num2str(numel(over25percentidx))]);
disp(['Number of Privately-Owned Reservoirs with over 25% loss, 2025: ', num2str(numel(intersect(private25idx,over25percentidx)))]);
disp(['NIDs of reservoirs with over 1km3 of sediment, 2025: ', strjoin(data.NID(over1kmidx), ', ')]);
disp(['NID names of reservoirs with over 1km3 of sediment, 2025: ', strjoin(data.NIDname(over1kmidx), ', ')]);
disp(['Number of headwater dams in 2025: ' num2str(numel(intersect(headresnet,rattes25)))])
disp(['Percentage of dams that are headwater dams in 2025: ' num2str(numel(intersect(headresnet,rattes25))/numel(rattes25)*100), ' %'])

disp('*');
disp(['Number of Large (>=100Mm3) reservoirs in 2025: ',num2str(numel(lgidx25))]);
disp(['Number of Small/Med (<100Mm3) reservoirs in 2025: ',num2str(numel(smidx25))]);
disp(['Median Year Built of Remaining Reservoirs in RATTES, that report yrc: ', num2str(medyrc_rates25)]);
disp(['Cumulative designed water storage of large reservoirs, 2025: ',num2str(storclg25),'km3']);
disp(['Cumulative designed water storage of small / medium reservoirs, 2025: ',num2str(storcsm25),'km3']);
disp(['Percent of total, designed water storage of small / medium reservoirs, 2025: ',num2str(storcsm25/storc25*100),'%']);

disp('*');
disp('****Sediment Volumes and Mass National, 2025****');
disp(['Total Sediment Volume at Remaining Reservoirs, 2025: ', num2str(sed25),' km^3']);
disp(['Percent of Sediment in Small and Mid-Sized Reservoirs, 2025: ', num2str(sed_smmed25/sed25*100),' %']);
disp(['Total Sediment Volume at Small and Mid-Sized Reservoirs, 2025: ', num2str(sed_smmed25),' km^3']);
disp(['Total Sediment Mass at Remaining Reservoirs, 2025: ', num2str(sedmass),' Gigatonnes']);
disp(['Multiple of global annual sediment load: ', num2str(times_globload)]);

disp('*');
disp('****Sedimentation Rates, 2025****');
disp(['Sedimentation Rate at Remaining Reservoirs in RATTES, 2025: ',num2str(seddt25), 'km3/yr']);
disp(['Sedimentation Rate at Large (>=100Mm3), Remaining Reservoirs in RATTES, 2025: ',num2str(seddtlg25), 'km3/yr']);
disp(['Sedimentation Rate at Small / Medium (<100Mm3), Remaining Reservoirs in RATTES, 2025: ',num2str(seddtsm25), 'km3/yr']);
disp(['Number of median sized reservoirs filled in each year: ',num2str(sum(calc.seddt(rattes25,now))/median(data.capcALL(rattes25)))]);

disp('*');
disp('****Capacity Loss, 2025****');
disp(['Total Capacity Loss at Remaining Reservoirs, 2025: ', num2str(caploss25),' %']);
disp(['Median Design Capacity in Remaining Reservoirs: ',num2str(median(data.capcALL(rattes25))), 'm3']);
disp(['Median % Capacity Loss at Remaining RATTES Reservoirs, 2025: ', num2str(median(calc.sed(rattes25,now)./data.capcALL(rattes25))*100),' %']);

disp('*');
disp('****Site data, 2025****');
disp(['Number of Remaining Reservoirs, 2025: ', num2str(numel(site25))]);
disp(['Median % Capacity Loss at Remaining Surveyed Reservoirs, 2025: ', num2str(median(calc.sed(site25,now)./data.capcALL(site25))*100),' %']);
disp(['Cumulative sedimentation rate at surveyed reservoirs, 2025: ', num2str(sum(calc.seddt(site25,now)/1e9)),' km3']);


%% 2025 sediment by distance

% fraction of sediment over 500 km from coast
faridx=find(mask_rattes & mask_far & mask_notremoved);
farsed=sum(calc.sed(faridx,now))/1e9;
perfarsed=farsed/sed25*100;

disp('*');
disp('****Queries about sediment distance from coast, 2025 results****');
disp(['Total Sediment Volume at Reservoirs more than 500 km from coast, 2025: ', num2str(farsed),' km^3']);
disp(['Percent of Sediment Volume at Reservoirs more than 500 km from coast, 2025: ', num2str(perfarsed),' %']);

%% 2025 terminal dams 

%terminal dams
termidx=find(mask_term & mask_rattes & mask_notremoved);
termsed25=sum(calc.sed(termidx,now))/1e9;
pertermsed25=termsed25/sed25*100;

disp('*');
disp('****Queries about terminal dams conus, 2025 results****');
disp(['Percent of Sediment Stored in Terminal Reservoirs in 2025: ', num2str(pertermsed25),'%']);
disp(['Volume of Sediment Stored in Terminal Reservoirs in 2025: ', num2str(termsed25), 'km3']);

%% Regional Analysis 2025 results
%masks: gulf, usatl, canatl, int, gl, mexpac, uspac
gulf25=find(mask_gulf & mask_rattes & mask_notremoved);
usatl25=find(mask_usatl & mask_rattes & mask_notremoved);
canatl25=find(mask_canatl & mask_rattes & mask_notremoved);
int25=find(mask_int & mask_rattes & mask_notremoved);
gl25=find(mask_gl & mask_rattes & mask_notremoved);
mexpac25=find(mask_mexpac & mask_rattes & mask_notremoved);
uspac25=find(mask_uspac & mask_rattes & mask_notremoved);
noteastgulf25=vertcat(canatl25,int25,gl25,mexpac25,uspac25);

gulfsedind25=calc.sed(gulf25,now);
gulfsedind25=sort(gulfsedind25,'descend');

smalleast=intersect(smidx25,usatl25);
smallgulf=intersect(smidx25,gulf25);
smallgl=intersect(smidx25,gl25);
smallcanatl=intersect(smidx25,canatl25);
smallmexpac=intersect(smidx25,mexpac25);
smallint=intersect(smidx25,int25);
smalluspac=intersect(smidx25,uspac25);

%dam density by region
%areas from GIS map for outlets 
area_uspac=974355.572376;
area_int=531503.267382;
area_gulf=4427874.14057;
area_mexpac=656556.286409;
area_canatl=210639.180114;
area_gl=454169.743943;
area_usatl=778793.107418;

disp('*');
disp('****Regional Queries, 2025 results****');
disp(['Dam Density for US to Pacific Outlet : ',num2str(numel(uspac25)/area_uspac*1000),'dams per 1000-km2']);
disp(['Dam Density for US Internally drained : ',num2str(numel(int25)/area_int*1000),'dams per 1000-km2']);
disp(['Dam Density for Gulf of Mexico Outlet : ',num2str(numel(gulf25)/area_gulf*1000),'dams per 1000-km2']);
disp(['Dam Density for Mexico to Pacific Outlet : ',num2str(numel(mexpac25)/area_mexpac*1000),'dams per 1000-km2']);
disp(['Dam Density for Canada to Atlantic Outlet : ',num2str(numel(canatl25)/area_canatl*1000),'dams per 1000-km2']);
disp(['Dam Density for Great Lakes Outlet : ',num2str(numel(gl25)/area_gl*1000),'dams per 1000-km2']);
disp(['Dam Density for US to Atlantic Outlet : ',num2str(numel(usatl25)/area_usatl*1000),'dams per 1000-km2']);

disp('*');
disp('****Total number of dams in RATTES by region, 2025 onward****');
disp(['Total number of dams in RATTES by region, 2025 onward, US to Pacific: ',numel(uspac25)]);
disp(['Total number of dams in RATTES by region, 2025 onward, Internally drained: ',numel(int25)]);
disp(['Total number of dams in RATTES by region, 2025 onward, Gulf: ',numel(mexpac25)]);
disp(['Total number of dams in RATTES by region, 2025 onward, Mexico to Pacific: ',numel(canatl25)]);
disp(['Total number of dams in RATTES by region, 2025 onward, Great Lakes: ',numel(gl25)]);
disp(['Total number of dams in RATTES by region, 2025 onward, US to Atlantic: ',numel(usatl25)]);

%gulf small med
disp('*');
disp('****Gulf****');
disp('*');
disp(['Percent of Sediment Stored Above Gulf in 2025: ', num2str((sum(gulfsedind25)/1e9)/sed25*100),'%']);
disp(['Percent of Gulf Sediment Stored in top 5 sediment storing reservoirs only : ', num2str(sum(gulfsedind25(1:5))/(sum(gulfsedind25))*100),'%']);
disp(['Number of small / med Gulf of Mexico Reservoirs in 2025 : ',num2str(numel(smallgulf))]);
disp(['Number of small / med Gulf of Mexico Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallgulf,over25percentidx)))]);
disp(['Percent of small / med Gulf of Mexico Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallgulf,over25percentidx))/numel(smallgulf)*100),'%'])
disp(['Number of all Gulf of Mexico Reservoirs with 25% cap loss : ',num2str(numel(intersect(gulf25,over25percentidx)))]);
disp(['Number of small / med Gulf of Mexico Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallgulf,over50percentidx)))]);
disp(['Percent of small / med Gulf of Mexico Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallgulf,over50percentidx))/numel(smallgulf)*100),'%'])
disp(['Number of all Gulf of Mexico Reservoirs with 50% cap loss : ',num2str(numel(intersect(gulf25,over50percentidx)))]);
disp(['Median Design Capacity, Gulf remaining in 2025: ',num2str(median(data.capcALL(gulf25)))]);
disp(['Median Design Capacity for small / med, Gulf remaining in 2025: ',num2str(median(data.capcALL(smallgulf)))]);
disp(['Median age of dam construction used by model, Gulf: ',num2str(median(data.yrc(gulf25)))]);
disp(['Percent of all Gulf of Mexico Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(gulf25,over50percentidx))/numel(gulf25)*100),'%'])
disp(['Percentage of Gulf Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(gulf25,over50percentidx),headresnet))/numel(intersect(gulf25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for Gulf Headwater Reservoirs with 50+ % capacity loss in 2025 : ', num2str((median(data.yrc(intersect(intersect(gulf25,over50percentidx),headresnet)))))])

%ustal sm-med
disp('*');
disp('****US to Atlantic****');
disp(['Number of small / med US to Atlantic Reservoirs in 2025 : ',num2str(numel(smalleast))]);
disp(['Number of small / med US to Atlantic Reservoirs with 25% cap loss : ',num2str(numel(intersect(smalleast,over25percentidx)))]);
disp(['Percent of small / med US to Atlantic Reservoirs with 25% cap loss : ',num2str(numel(intersect(smalleast,over25percentidx))/numel(smalleast)*100),'%'])
disp(['Number of small / med US to Atlantic Reservoirs with 50% cap loss : ',num2str(numel(intersect(smalleast,over50percentidx)))]);
disp(['Percent of small / med US to Atlantic Reservoirs with 50% cap loss : ',num2str(numel(intersect(smalleast,over50percentidx))/numel(smalleast)*100),'%'])%gulf sm-med
disp(['Median Design Capacity, US to Atlantic remaining in 2025: ',num2str(median(data.capcALL(usatl25)))]);
disp(['Median Design Capacity for small / med, US to Atlantic remaining in 2025: ',num2str(median(data.capcALL(smalleast)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, US to Atlantic: ',num2str(median(data.yrc(usatl25)))]);
disp(['Percent of all US to Atlantic Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(usatl25,over50percentidx))/numel(usatl25)*100),'%']);
disp(['Percentage of US to Atlantic Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(usatl25,over50percentidx),headresnet))/numel(intersect(usatl25,over50percentidx))*100),'%']);
disp(['Median Year of Dam Completion for US to Atlantic Headwater Reservoirs with 50+ % capacity loss in 2025 : ', num2str((median(data.yrc(intersect(intersect(usatl25,over50percentidx),headresnet)))))])

%great lakes sm-med
disp('*');
disp('****Great Lakes****');
disp(['Number of small / med Great Lakes Reservoirs in 2025 : ',num2str(numel(smallgl))]);
disp(['Number of small / med Great Lakes Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallgl,over25percentidx)))]);
disp(['Percent of small / med Great Lakes Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallgl,over25percentidx))/numel(smallgl)*100),'%'])
disp(['Number of small / med Great Lakes Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallgl,over50percentidx)))]);
disp(['Percent of small / med Great Lakes Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallgl,over50percentidx))/numel(smallgl)*100),'%'])
disp(['Median Design Capacity, Great Lakes remaining in 2025: ',num2str(median(data.capcALL(gl25)))]);
disp(['Median Design Capacity for small / med, Great Lakes remaining in 2025: ',num2str(median(data.capcALL(smallgl)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, Great Lakes: ',num2str(median(data.yrc(gl25)))]);
disp(['Percent of all Great Lakes Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(gl25,over50percentidx))/numel(gl25)*100),'%'])
disp(['Percentage of Great Lakes Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(gl25,over50percentidx),headresnet))/numel(intersect(gl25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for Great Lakes Headwater Reservoirs with 50+ % capacity loss in 2025: ', num2str((median(data.yrc(intersect(intersect(gl25,over50percentidx),headresnet)))))])

%canada atlantic sm-med
disp('*');
disp('****Canada to Atlantic****');
disp(['Number of small / med Canada to Atlantic Reservoirs in 2025 : ',num2str(numel(smallcanatl))]);
disp(['Number of small / med Canada to Atlantic Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallcanatl,over25percentidx)))]);
disp(['Percent of small / med Canada to Atlantic Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallcanatl,over25percentidx))/numel(smallcanatl)*100),'%'])
disp(['Number of small / med Canada to Atlantic Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallcanatl,over50percentidx)))]);
disp(['Percent of small / med Canada to Atlantic Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallcanatl,over50percentidx))/numel(smallcanatl)*100),'%'])
disp(['Median Design Capacity, Canada atlantic remaining in 2025: ',num2str(median(data.capcALL(canatl25)))]);
disp(['Median Design Capacity for small / med, Canada atlantic remaining in 2025: ',num2str(median(data.capcALL(smallcanatl)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, Canada atlantic: ',num2str(median(data.yrc(canatl25)))]);
disp(['Percent of all Canada to Atlantic Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(canatl25,over50percentidx))/numel(canatl25)*100),'%'])
disp(['Percentage of Canada to Atlantic Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(canatl25,over50percentidx),headresnet))/numel(intersect(canatl25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for Canada to Atlantic Headwater Reservoirs with 50+ % capacity loss in 2025: ', num2str((median(data.yrc(intersect(intersect(canatl25,over50percentidx),headresnet)))))])

%mexico pacific sm-med
disp('*');
disp('****Mexico to Pacific****');
disp(['Number of small / med Mexico to Pacific Reservoirs in 2025 : ',num2str(numel(smallmexpac))]);
disp(['Number of small / med Mexico to Pacific Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallmexpac,over25percentidx)))]);
disp(['Percent of small / med Mexico to Pacific Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallmexpac,over25percentidx))/numel(smallmexpac)*100),'%'])
disp(['Number of small / med Mexico to Pacific Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallmexpac,over50percentidx)))]);
disp(['Percent of small / med Mexico to Pacific Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallmexpac,over50percentidx))/numel(smallmexpac)*100),'%'])
disp(['Median Design Capacity, Mexico Pacific remaining in 2025: ',num2str(median(data.capcALL(mexpac25)))]);
disp(['Median Design Capacity for small / med, Mexico Pacific remaining in 2025: ',num2str(median(data.capcALL(smallmexpac)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, Mexico Pacific: ',num2str(median(data.yrc(mexpac25)))]);
disp(['Percent of all Mexico to Pacific Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(mexpac25,over50percentidx))/numel(mexpac25)*100),'%'])
disp(['Percentage of Mexico to Pacific Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(mexpac25,over50percentidx),headresnet))/numel(intersect(mexpac25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for Mexico to Pacific Headwater Reservoirs with 50+ % capacity loss in 2025: ', num2str((median(data.yrc(intersect(intersect(mexpac25,over50percentidx),headresnet)))))])

%US pacific sm-med
disp('*');
disp(['Number of small / med US to Pacific Reservoirs in 2025 : ',num2str(numel(smalluspac))]);
disp(['Number of small / med US to Pacific Reservoirs with 25% cap loss : ',num2str(numel(intersect(smalluspac,over25percentidx)))]);
disp(['Percent of small / med US to Pacific Reservoirs with 25% cap loss : ',num2str(numel(intersect(smalluspac,over25percentidx))/numel(smalluspac)*100),'%'])
disp(['Number of small / med US to Pacific Reservoirs with 50% cap loss : ',num2str(numel(intersect(smalluspac,over50percentidx)))]);
disp(['Percent of small / med US to Pacific Reservoirs with 50% cap loss : ',num2str(numel(intersect(smalluspac,over50percentidx))/numel(smalluspac)*100),'%'])
disp(['Median Design Capacity, US to Pacific remaining in 2025: ',num2str(median(data.capcALL(uspac25)))]);
disp(['Median Design Capacity for small / med, US to Pacific remaining in 2025: ',num2str(median(data.capcALL(smalluspac)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, US to Pacific: ',num2str(median(data.yrc(uspac25)))]);
disp(['Percent of all US to Pacific Reservoirs with 50% cap loss in 2025: ',num2str(numel(intersect(uspac25,over50percentidx))/numel(uspac25)*100),'%'])
disp(['Percentage of US to Pacific Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(uspac25,over50percentidx),headresnet))/numel(intersect(uspac25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for US to Pacific Headwater Reservoirs with 50+ % capacity loss in 2025: ', num2str((median(data.yrc(intersect(intersect(uspac25,over50percentidx),headresnet)))))])


%internal sm-med
disp('*');
disp('****Internally Drained****');
disp(['Number of small / med Internally Drained Reservoirs in 2025 : ',num2str(numel(smallint))]);
disp(['Number of small / med Internally Drained Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallint,over25percentidx)))]);
disp(['Percent of small / med Internally Drained Reservoirs with 25% cap loss : ',num2str(numel(intersect(smallint,over25percentidx))/numel(smallint)*100),'%'])
disp(['Number of small / med Internally Drained Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallint,over50percentidx)))]);
disp(['Percent of small / med Internally Drained Reservoirs with 50% cap loss : ',num2str(numel(intersect(smallint,over50percentidx))/numel(smallint)*100),'%'])
disp(['Median Design Capacity, Internally Drained remaining in 2025: ',num2str(median(data.capcALL(int25)))]);
disp(['Median Design Capacity for small / med, Internally Drained remaining in 2025: ',num2str(median(data.capcALL(smallint)))]);
disp(['Median age of dam construction used by model including yrc subsituted data, Internally Drained: ',num2str(median(data.yrc(int25)))]);
disp(['Percent of all Internally Drained Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(int25,over50percentidx))/numel(int25)*100),'%'])
disp(['Percentage of Internally Drained Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(int25,over50percentidx),headresnet))/numel(intersect(int25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for Internally Drained Headwater Reservoirs with 50+ % capacity loss in 2025 : ', num2str((median(data.yrc(intersect(intersect(int25,over50percentidx),headresnet)))))])

%CONUS sm-med
disp('*');
disp('****CONUS****');
disp(['Number of small / med RATTES Reservoirs in 2025 : ',num2str(numel(smidx25))]);
disp(['Number of small / med RATTES Reservoirs with 25% cap loss : ',num2str(numel(intersect(smidx25,over25percentidx)))]);
disp(['Percent of small / med RATTES Reservoirs with 25% cap loss : ',num2str(numel(intersect(smidx25,over25percentidx))/numel(smidx25)*100),'%'])
disp(['Number of small / med RATTES Reservoirs with 50% cap loss : ',num2str(numel(intersect(smidx25,over50percentidx)))]);
disp(['Percent of small / med RATTES Reservoirs with 50% cap loss : ',num2str(numel(intersect(smidx25,over50percentidx))/numel(smidx25)*100),'%'])%conus
disp(['Median Design Capacity, RATTES remaining in 2025: ',num2str(median(data.capcALL(rattes25)))]);
disp(['Median Design Capacity of small / med, RATTES remaining in 2025: ',num2str(median(data.capcALL(smidx25)))]);
disp(['Median age of dam construction used by model, RATTES- including yrc subsituted data: ',num2str(median(data.yrc(rattes25)))]);
disp(['Percent of all RATTES Reservoirs with 50% cap loss in 2025 : ',num2str(numel(intersect(rattes25,over50percentidx))/numel(rattes25)*100),'%'])
disp(['Percentage of RATTES Reservoirs with 50% cap loss that are also headwater in 2025 : ',num2str(numel(intersect(intersect(rattes25,over50percentidx),headresnet))/numel(intersect(rattes25,over50percentidx))*100),'%'])
disp(['Median Year of Dam Completion for RATTES Headwater Reservoirs with 50+ % capacity loss in 2025 : ', num2str((median(data.yrc(intersect(intersect(rattes25,over50percentidx),headresnet)))))])


%Sedrates by region
disp('*');
disp(['Sedimentation Rate in Gulf of Mexico 2025 : ',num2str(sum(calc.seddt(gulf25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in US to Atlantic 2025 : ',num2str(sum(calc.seddt(usatl25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in Great Lakes 2025 : ',num2str(sum(calc.seddt(gl25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in Canada to Atlantic 2025 : ',num2str(sum(calc.seddt(canatl25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in Mexico to Pacific 2025 : ',num2str(sum(calc.seddt(mexpac25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in US to Pacific 2025 : ',num2str(sum(calc.seddt(uspac25,now))/1e6), ' M m3/yr']);
disp(['Sedimentation Rate in Internally Drained 2025 : ',num2str(sum(calc.seddt(int25,now))/1e6), ' M m3/yr']);

%Yield rates by region
disp('*');
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, Gulf of Mexico 2025 : ',num2str(sum(calc.seddt(gulf25,now))/sum(calc.wsedshed25(gulf25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, US Atlantic 2025 : ',num2str(sum(calc.seddt(usatl25,now))/sum(calc.wsedshed25(usatl25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, Great Lakes 2025 : ',num2str(sum(calc.seddt(gl25,now))/sum(calc.wsedshed25(gl25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, Canada to Atlantic 2025 : ',num2str(sum(calc.seddt(canatl25,now))/sum(calc.wsedshed25(canatl25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, Mexico to Pacific 2025 : ',num2str(sum(calc.seddt(mexpac25,now))/sum(calc.wsedshed25(mexpac25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, US to Pacific 2025 : ',num2str(sum(calc.seddt(uspac25,now))/sum(calc.wsedshed25(uspac25))), ' m3/km2/yr']);
disp(['Regional Sediment Yield per unit of effective contributing drainage area, RATTES, Internally Drained 2025 : ',num2str(sum(calc.seddt(int25,now))/sum(calc.wsedshed25(int25))), ' m3/km2/yr']);
disp(['National Sediment Yield per unit of effective contributing drainage area, RATTES 2025 : ',num2str(sum(calc.seddt(rattes25,now))/sum(calc.wsedshed25(rattes25))), ' m3/km2/yr']);

disp('*');
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, Gulf of Mexico 2025 : ',num2str(median(calc.wYield25(gulf25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, US Atlantic 2025 : ',num2str(median(calc.wYield25(usatl25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, Great Lakes 2025 : ',num2str(median(calc.wYield25(gl25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, Canada to Atlantic 2025 : ',num2str(median(calc.wYield25(canatl25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, Mexico to Pacific 2025 : ',num2str(median(calc.wYield25(mexpac25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, US to Pacific 2025 : ',num2str(median(calc.wYield25(uspac25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, Internally Drained 2025 : ',num2str(median(calc.wYield25(int25))), ' m3/km2/yr']);
disp(['Median Sediment Yield per unit of effective contributing drainage area, RATTES, National 2025 : ',num2str(median(calc.wYield25(rattes25))), ' m3/km2/yr']);


%% Delta 2025 data & rivers

%Delta IDs
mississippi=16;
colorado=23;
rio=22;
columbia=25;
elwha=28;

%2025 sealevel rise
slr25=calcdelta.sealevelrise(now); 

%deltas with more than 1 km3 upstream
deltw1km3=find((calcdelta.sedabovedeltas(:,now)/1e9)>=1);

%active lobe area mississippi
actareamiss=12730*(1000*1000); % km2 converted to m2

%deltatrapping
deltatrap=0.3;

%Adjust Miss river equiv heights for just the active lobes
adjustmiss_H=calcdelta.sedabovedeltas(mississippi,now)/actareamiss*deltatrap;%meters
adjustmiss_hdt=(calcdelta.seddtabovedeltas(mississippi,now)/actareamiss)*1000*deltatrap; %mm/yr

%Elwha before dam removals
t2011=find(calc.t==2011);

disp('*');
disp('****Delta Data Queries, 2025 results****');
disp(['Total Sediment Volume above Mississippi River Delta, 2025: ', num2str(calcdelta.sedabovedeltas(mississippi,now)/1e9),' km^3']);
disp(['Total Sediment Volume above Colorado River Delta, 2025: ', num2str(calcdelta.sedabovedeltas(colorado,now)/1e9),' km^3']);
disp(['Total Sediment Volume above Columbia River Delta (in US only), 2025: ', num2str(calcdelta.sedabovedeltas(columbia,now)/1e9),' km^3']);
disp(['Total Sediment Volume above Rio Grande River Delta (in US only), 2025: ', num2str(calcdelta.sedabovedeltas(rio,now)/1e9),' km^3']);
disp(['Number of deltas where the potential equivalent aggradation depth outpaces sea level rise, 2025: ', num2str(numel(find((calcdelta.EquivDepth_dt(:,now)*1000)>=slr25)))]);
disp(['Numbers of deltas where upstream reservoirs trap more than 1km3, 2025: ', num2str(numel(deltw1km3))]);
disp(['Numbers of deltas where upstream reservoirs trap more than 1km3 that are located along Gulf Coast, 2025: ', num2str(numel(find(datadelta.countryOut(deltw1km3)==5)))]);
disp(['Names of deltas where upstream reservoirs trap more than 1km3, 2025: ', strjoin(datadelta.DeltaName(deltw1km3))]);
disp(['Cumulative number of dams above deltas storing more than 1km3, 2025: ', num2str(sum(calcdelta.DamsAboveDeltas(deltw1km3,now)))]);
disp(['Number of deltas on free-flowing river systems (no dams upstreams): ', num2str(numel(find(calcdelta.DamsAboveDeltas(:,now)==0)))]);
disp(['Equivalent height (H) on active lobes of Mississippi River Delta: ', num2str(adjustmiss_H),' meters']);
disp(['Annual equivalent height (h) on active lobes of Mississippi River Delta: ', num2str(adjustmiss_hdt),' mm/yr']);
disp(['Sedimentation rate in 2011 above Elwha delta prior to dam removals: ', num2str(calcdelta.seddtabovedeltas(elwha,t2011)/1e6),' M m3/yr']);

%Rivers
t2023=find(calc.t==2023);
klaremidx=find(data.SID==288941 | data.SID==37025 | data.SID==37027 | data.SID==274759); %sids of removed klamath dams

disp('*');
disp('****River Data Queries, 2025 results****');
disp(['Sedimentation rate in 2023 above removed Klamath mouth delta prior to removal: ', num2str(sum(calc.seddt(klaremidx,t2023)/1e6)),' M m3/yr']);

%% 2025 value of sand and gravel

% Value of sand and gravel
dolperton=(12e9)/(890e6); %12 billion dollars for 890 million tons of sand and gravel
value15=(sedmass/1e9)*0.15*dolperton*1e9; %15% of mass is sand and gravel- put sed mass in tons (from gigatons), keep 15%, calc dollars per ton, report billons of dollars
disp(['Value of Sand and Gravel, if 15% concentration:is sand and gravel $',num2str(value15)]);


%% 2050 results capacity loss

sed50=sum(calc.sed(rattes25,end))/1e9; %sum of sediment volume in 2025, km3
storc50=sum(data.capcALL(rattes25))/1e9; % design storage for reservoirs that still exist, km3
caploss50=sed50/storc50*100; % percent capacity loss
%sed rate
seddt50=sum(calc.seddt(rattes25,end))/1e9;
disp('*');
disp('*');
disp('****QUERIES ABOUT RATTES 2050 RESULTS****');
disp(['Total Sediment Volume at Remaining Reservoirs, 2050: ', num2str(sed50),' km^3']);
disp(['Total Capacity Loss at Remaining Reservoirs, 2050: ', num2str(caploss50),' %']);
disp(['Sedimentation Rate at Remaining Reservoirs in RATTES, 2050: ',num2str(seddt50), 'km3/yr']);
disp(['Median % Capacity Loss at Remaining RATTES Reservoirs, 2050: ', num2str(median(calc.sed(rattes25,end)./data.capcALL(rattes25))*100),' %']);

%% Create Supplemental Delta table
numdelta = height(datadelta);  % 
main=calcdelta.maindeltas;

deltasupp.deltaID_Resnet = datadelta.deltaID;
deltasupp.delta_name    = datadelta.DeltaName;

deltasupp.InMainText = zeros(numdelta,1);
deltasupp.InMainText(main) = 1;

deltasupp.Delta_trapping = 0.3 * ones(numdelta,1);
deltasupp.PercentBasinTrapped = calcdelta.PercentUpstreamBasinTrapping(:,now);
deltasupp.NumDams = calcdelta.DamsAboveDeltas(:,now);
deltasupp.SedVol_Mm3 = calcdelta.sedabovedeltas(:,now) / 1e6;
deltasupp.deltaArea = datadelta.Area;
deltasupp.totH = calcdelta.EquivDepth_Total(:,now);
deltasupp.sedrate_Mm3peryr = calcdelta.seddtabovedeltas(:,now) / 1e6;
deltasupp.annual_h_mmperyr = calcdelta.EquivDepth_dt(:,now) * 1000;
deltasupp.uptermdams = calcdelta.maxTermDamsAboveDeltas;
deltasupp.percentsedinterm = calcdelta.PercentSedInTermDams(:,now);

for k = 1:5
    fname = sprintf('TopSed%d', k);
    deltasupp.(fname) = zeros(numdelta,1);
end

topsed=NaN(numdelta,5);
sedvols=NaN(numdelta,5);

%Top 5 sed storing reservoirs
for i=datadelta.deltaID'
    mask_tmp=(data.deltatag==i);
    idx=find(mask_tmp & mask_rattes & mask_notremoved);
    junk(:,1)=data.SID(idx);
    junk(:,2)=calc.sed(idx,now);
    if isempty(junk)~=1
    junk=sortrows(junk,2,'descend');
    n = min(5, size(junk,1));   % how many actually exist
    topsed(i,1:n) = junk(1:n,1)';
    sedvols(i,1:n) = (junk(1:n,2)')/1e6;
    end
    clear junk idx
end
for k = 1:5
    deltasupp.(sprintf('TopSed%d',k)) = topsed(:,k);
end

%Overwrite Elwha entries for 2011 data, prior to dam removal
deltasupp.PercentBasinTrapped(elwha) = calcdelta.PercentUpstreamBasinTrapping(elwha,t2011);
deltasupp.NumDams(elwha) = calcdelta.DamsAboveDeltas(elwha,t2011);
deltasupp.SedVol_Mm3(elwha) = calcdelta.sedabovedeltas(elwha,t2011) / 1e6;
deltasupp.totH(elwha) = calcdelta.EquivDepth_Total(elwha,t2011);
deltasupp.sedrate_Mm3peryr(elwha) = calcdelta.seddtabovedeltas(elwha,t2011) / 1e6;
deltasupp.annual_h_mmperyr(elwha) = calcdelta.EquivDepth_dt(elwha,t2011) * 1000;
deltasupp.percentsedinterm(elwha) = calcdelta.PercentSedInTermDams(elwha,t2011);

for i=elwha
    mask_tmp=(data.deltatag==i);
    idx=find(mask_tmp & mask_rattes);
    junk(:,1)=data.SID(idx);
    junk(:,2)=calc.sed(idx,t2011);
    if isempty(junk)~=1
    junk=sortrows(junk,2,'descend');
    n = min(5, size(junk,1));   % how many actually exist
    topsed(i,1:n) = junk(1:n,1)';
    sedvols(i,1:n) = (junk(1:n,2)')/1e6;
    end
    clear junk idx
end
for k = 1:5
    deltasupp.(sprintf('TopSed%d',k)) = topsed(:,k);
end

deltasupp=struct2table(deltasupp);

%% create Zenodo files
usa=find(mask_usa);
numdam=numel(usa);
numt=numel(calc.t);
%timeseries output, for MLR input
Sediment_m3=NaN((numdam+1),(numt+1));
Sediment_m3(1,2:end)=calc.t;
Sediment_m3(2:end,1)=data.SID(usa);

SedimentHi_m3=Sediment_m3; SedimentLo_m3=Sediment_m3;
Capacity_m3=Sediment_m3;CapacityHi_m3=Sediment_m3;CapacityLo_m3=Sediment_m3;

Sediment_m3(2:end,2:end)=calc.sed(usa,:);
SedimentHi_m3(2:end,2:end)=calc.sedup(usa,:);
SedimentLo_m3(2:end,2:end)=calc.sedlo(usa,:);

Capacity_m3(2:end,2:end)=calc.cap(usa,:);
CapacityHi_m3(2:end,2:end)=calc.capup(usa,:);
CapacityLo_m3(2:end,2:end)=calc.caplo(usa,:);

writematrix(Sediment_m3,'Sediment_m3_010626.csv','Delimiter',',');
writematrix(SedimentHi_m3,'SedimentHi_m3_010626.csv','Delimiter',',');
writematrix(SedimentLo_m3,'SedimentLo_m3_010626.csv','Delimiter',',');

writematrix(Capacity_m3,'Capacity_m3_010626.csv','Delimiter',',');
writematrix(CapacityHi_m3,'CapacityHi_m3_010626.csv','Delimiter',',');
writematrix(CapacityLo_m3,'CapacityLo_m3_010626.csv','Delimiter',',');

end