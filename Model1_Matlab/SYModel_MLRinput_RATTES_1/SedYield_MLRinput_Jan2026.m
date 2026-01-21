%% Sediment Yield model (following Minear and Kondolf) and file prep for Multiple Linear Regression model component of RATTES
clear
tic

loadpart1results="no";

if loadpart1results=="yes"
    load SYModelResults_MLRinput.mat
    return
end
 clear loadpart1results
%% setup

inputdate='052225';
inputdateMLR='010626';

inputs= "no";
RankTag= "no";
setcalc= "no";
FirstYield="no";
FinalYield="no";
MLRinput="no";

%% conversions
convert1=1233.482; %converts m3 to ac-ft is convert1*AF=m3, or from m3 is m3/convert1=AF
convert2= 2.59; %converts km2 and mi2 is convert2*mi2= km2 or from km2 is km2/convert2=mi2

%% Add inputs and data cleaning, fill in missing yrc within this function.
if inputs == "yes"
    fmt="ResNetInput_SitesCanada_%s.csv";
    filename=sprintf(fmt,inputdate); clear fmt;
    [InputAll] = AddInputsUpdates(filename);
    clear filename; 
else
    load Inputs.mat;
end
clear inputs
%% Rank dams from headwaters to terminal dams
if RankTag == "yes"
    [data1wrivers,data2norivers,data3SYmodel] = RankTagApr25(InputAll);
else
    load RankTag.mat;
end
clear RankTag

%% Set up initial calc matrices, fill with pre & post removal zeros. 

if setcalc == "yes"
    [calc2_MLRinput,calc3_SY1,calc4_SY2] = SetCalcMatrix(data2norivers,convert1,convert2); %calc2 goes with data2; calc3 goes with data3
else
    load setcalc.mat
end
clear setcalc

%% Sediment yield function

if FirstYield == "yes"
    [data3SYmodel,calc3_SY1] = FirstYieldFunctionApr25(data3SYmodel,calc3_SY1,convert1,convert2);
else
    load Yield1.mat;
end
clear FirstYield

%% Final Yield Function

if FinalYield == "yes"
    [data4SYmodel,calc4_SY2] = FinalYieldFunctionApr25(data3SYmodel,calc4_SY2,calc3_SY1,convert1,convert2);
else
    load YieldFinal.mat
end
clear FinalYield 


%% MLR input- create input files for Multiple Linear Regression Model
if MLRinput == "yes"
    fmt='MLR_Input%s.csv';
    MLRfilename=sprintf(fmt,inputdate); clear fmt;
    fmt='MLR_SAatDamInitial%s.csv';
    MLRfilename2=sprintf(fmt,inputdate); clear fmt;
    [calc2_MLRinput] = MLRinputApr25(convert1,convert2,data2norivers,MLRfilename,MLRfilename2,calc2_MLRinput);
    clear MLRfilename; clear MLRfilename2;
else
    load MLRdata.mat;
end
clear MLRinput;

%% Save SY Model and MLR inputs for later model combination 

clear inputdate inputdateMLR InputAll
save('SYModelResults_MLRinput.mat','-v7.3');

toc


