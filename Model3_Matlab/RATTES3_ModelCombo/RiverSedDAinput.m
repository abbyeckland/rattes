function [data1wrivers,calc1all_wrivers] = RiverSedDAinput(data1wrivers,calc5all,data5all)%

%% set up data

data=data1wrivers;
a=find(data.SID<0);
numriv=numel(a);
damstart=numriv+1;
numall=numel(data.SID);

%fill in with the updated reservoir data for end columns and for capc
%(which has back projections)
data.inMLR(:)=0;
data.inMLR(damstart:numall)=data5all.inMLR;
data.inSY(:)=0;
data.inSY(damstart:numall)=data5all.inSY;
data.capcMLR(:)=0;
data.capcMLR(damstart:numall)=data5all.capcMLR;
data.capcALL(:)=0;
data.capcALL(damstart:numall)=data5all.capcALL;
data.capcSY(:)=0;
data.capcSY(damstart:numall)=data5all.capcSY;
data.compare(:)=0;
data.compare(damstart:numall)=data5all.compare;
data.maxstor=[];

%% set up calc
%create calc matrix by prepending rows to the top for data that won't
%change 
calc=calc5all;
numt=numel(calc.t);

% Number of rows to prepend
nPrepend = numriv;

% Get all field names
fields = fieldnames(calc);

% Loop over each field
for i = 1:numel(fields)
    field = fields{i};

    % Skip the 't' field
    if strcmp(field, 't')
        continue
    end
    
    % Get current data
    calcdata = calc.(field);
    
    % Determine number of columns
    nCols = size(calcdata,2);
    
    % Prepend zeros (or NaNs if preferred)
    prependData = zeros(nPrepend, nCols); % or NaN(nPrepend,nCols)
    
    % Concatenate
    calc.(field) = [prependData; calcdata];
end

%Infill with river data where can't be 0
builtriver=find(calc.t==data.yrc(1)); %all rivers have same "year built" 
calc.SID(1:numriv)=data.SID(1:numriv);
calc.builtidx(1:numriv)=builtriver;
calc.removedidx(1:numriv)=NaN;
calc.sur1idx(1:numriv)=NaN;
calc.predictidx(1:numriv)=builtriver;

calc.origDA(1:numriv,1)=data.DA(1:numriv);
calc.origDA(:,2:numt) = repmat(calc.origDA(:,1), 1, numt-1);

%check make sure SID matches
a=find(data.SID ~= calc.SID);
if isempty(a)~=1
    Error="SID mismatch"
    return
end

%% finalize 

calc1all_wrivers=calc;
data1wrivers=data;

save('riversetup.mat','calc1all_wrivers','data1wrivers','-v7.3');
disp('*')
disp('Done with setting up River struct')
end


