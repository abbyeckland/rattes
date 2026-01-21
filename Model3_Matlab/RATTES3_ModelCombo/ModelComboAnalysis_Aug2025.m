%Combine MLR and SY models, back project MLR results to year of dam
%completion if yrp>yrc, analysis
%figure and table creation
%figures will save into pre-set 'Figures' subfolders

clear
tic

load SYModelResults_MLRinput.mat
loadresults2="no"; %load final results if don't need to rerun

if loadresults2=="yes"
    load FinalRATTEScombo.mat
    clear loadresults2;
    return
end
clear loadresults2;
%% setup
inputdateMLR='010626'; %MLR model run date

modelcombo="no";
compare="no";
riversetup="no";
riverDA="no";
termdeltadams="no";
queries="no";
paperfigs="no";

fignum=1; %figure number for plotting

%% conversions
%conversions,
convert1=1233.482; %converts m3 to ac-ft is convert1*AF=m3, or from m3 is m3/convert1=AF
convert2= 2.59; %converts km2 and mi2 is convert2*mi2= km2 or from km2 is km2/convert2=mi2

%% Calc MLR backwards and combine with site data to match original input dam length
%need to add 95% confidence interval going backwards

if modelcombo == "yes"
    fmt='mlr_trap_decimal_%s.csv';
    filenameMLRtrap=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_up_95_%s.csv';
    filenameMLRcapup=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_low_95_%s.csv';
    filenameMLRcaplow=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_%s.csv';
    filenameMLRcap=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_up_95_%s.csv';
    filenameMLRsvup=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_low_95_%s.csv';
    filenameMLRsvlow=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_%s.csv';
    filenameMLRsv=sprintf(fmt,inputdateMLR); clear fmt;
   [data5all,calc5all] = CombineMLRSY(filenameMLRtrap,filenameMLRcap,filenameMLRcapup,filenameMLRcaplow,filenameMLRsv,filenameMLRsvup,filenameMLRsvlow,data2norivers,calc2_MLRinput,data4SYmodel,calc4_SY2,convert1,convert2);
   clear filenameMLRtrap filenameMLRcap filenameMLRcapup filenameMLRcaplow filenameMLRsv filenameMLRsvup filenameMLRsvlow
else
    load ComboModel.mat
end
clear modelcombo;

%% compare model results
%at dams with perm storage in the US only. 
if compare=="yes"
    [fignum,calc5all,data5all] = CompareSY_MLRAug25(calc5all,data5all,fignum);
else
    load CompareModels.mat;
end
clear compare;

%% Calc sed contributing DA at rivers and all dams
%includes all dams, including lock dams without permanent storage.

if riversetup == "yes"
    [data1wrivers,calc1all_wrivers] = RiverSedDAinput(data1wrivers,calc5all,data5all);%
else
    load riversetup.mat
end
clear riversetup;

%% Calc sed contributing DA at rivers, and upstream sediment
%Sed contributing DA includes all dams, including lock dams without
%permanent storage in US for sediment and yield

if riverDA == "yes"
    [calc1all_wrivers,data1wrivers,data6_delta,calc6_delta,calc5all,fignum] = SDArivers(calc1all_wrivers,data1wrivers,calc5all,data5all,fignum);
    clear doneriver;
else
    load riverSDA.mat
end
clear riverDA

%% Terminal Dams & number of dams above rivers
% this is with permanent storage us dams
if termdeltadams=="yes"
    [calc1all_wrivers,calc6_delta] = DeltaTerm(calc1all_wrivers,calc6_delta,calc5all,data5all,data1wrivers,data6_delta);
else
    load deltaterm.mat
end
clear termdeltadams


%% Queries for paper
if queries=="yes"
    RATTESqueries(data5all,calc5all,data6_delta,calc6_delta);
end
clear queries

%% Eckland et al Fig 2
if paperfigs=="yes"
    [fignum] = Eckland_PaperFigure2(data5all,calc5all,fignum);
    [fignum] = Eckland_PaperFigure3(data5all,calc5all,fignum,calc6_delta,data6_delta);
end
clear paperfigs

%%
clear inputdateMLR
fignum=1;
save('FinalRATTEScombo.mat','-v7.3');

toc