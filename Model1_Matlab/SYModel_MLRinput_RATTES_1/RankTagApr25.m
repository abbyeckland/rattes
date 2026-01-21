%2024 dam ordering for the entire nation. This will probably become the
%data cleaning function

%tag and rank dams with and without rivers

function [data1wrivers,data2norivers,data3SYmodel] = RankTagApr25(InputAll)

data=InputAll;
nRows=height(data);
rivers=find(data.isriver==1);
%% Create additional river tags for rivers downstream from rivers

tags(:,1)=data.rivertag1;
j=1;
tag=find(tags(:,j)~=0);


%Make a rivertags cell
data.rivertags=num2cell(data.rivertag1);

while ~isempty(tag)
    tags(:,(j+1))=zeros(nRows,1);
    for taglocation = tag'
        SIDoftaggedriver = tags(taglocation, j);
        taggedriverlocation=find(data.SID==SIDoftaggedriver); %row of the tagged river
        tags(taglocation,(j+1))=tags(taggedriverlocation,1);%found original location of the river- its tag will be in column 1, if one exists.
    end
    j=j+1;
    temp = tags(:, j);
    tag=find(tags(:,j)~=0); %rows that have river tags
    for idx=tag'
        sts_mat=cell2mat(data.rivertags(idx));
        sts_mat(end + 1) = temp(idx);
        data.rivertags(idx) = {sts_mat};
    end
    clear idx; clear sts_mat; clear new; clear temp;
end

vars={'j','tags','taglocation','SIDoftaggedriver','taggedriverlocation','tag'};
clear(vars{:}); clear vars;

%% create a additional site tags for sites downstream from other sites.

tags(:,1)=data.sitetag1;
j=1;
tag=find(tags(:,j)~=0);


%Make a sitetags cell
data.sitetags=num2cell(data.sitetag1);

while ~isempty(tag)
    tags(:,(j+1))=zeros(nRows,1);
    for taglocation = tag'
        SIDoftaggedsite = tags(taglocation, j);
        taggedsitelocation=find(data.SID==SIDoftaggedsite); %row of the tagged site
        tags(taglocation,(j+1))=tags(taggedsitelocation,1);%found original location of the site- its tag will be in column 1, if one exists.
    end
    j=j+1;
    temp = tags(:, j);
    tag=find(tags(:,j)>0); %rows that have site tags
    for idx=tag'
        sts_mat=cell2mat(data.sitetags(idx));
        sts_mat(end + 1) = temp(idx);
        data.sitetags(idx) = {sts_mat};
    end
    clear idx; clear sts_mat; clear new; clear temp;
end

vars={'j','tags','taglocation','SIDoftaggedsite','taggedsitelocation','tag'};
clear(vars{:}); clear vars;

%% Make a ranking file of dam order WITH rivers
%Make a ranking file where headwater dams = rank 1, and so on down the line
% Ranking will be incorrect for dams (b/c includes rivers) but SA will be correct
  
Rank = NaN(nRows, 1);
Rank(data.head == 1) = 1;
DAerrorNumber = 0;

i=1;
ranknum=find(Rank==i);

while ~isempty(ranknum)
    jmax=numel(ranknum); % number of dams with a rank of i
    for j=1:jmax
        thisdam=ranknum(j); % The index we are on within array j
        if data.SID(thisdam)<0
            toriv=find(data.todam==data.SID(thisdam));
            totDAtoriv=sum(data.DA(toriv));
            if totDAtoriv>data.DA(thisdam)
                data.DA(thisdam)=totDAtoriv;
            end
        end
        %process nonterminal dams
        if  data.term(thisdam)==0 % if dam j is not a terminal dam, then
            thisdamDA=data.DA(thisdam); % identify DA of j
            goesto=data.todam(thisdam); % identify where dam j goes to
            thatdam=find(data.SID==goesto); % find the SID of the dam that dam j goes to
            thatdamDA=data.DA(thatdam); % identify the DA that dam j goes to

            %do a check to make sure downstream DA doesn't exceed upstream
            if (thatdamDA<thisdamDA && data.SID(thatdam)>0) %identify if upstream drainage area is greater than downstream drainage area, and not a river, b/c that will be fixed
                DAerrorNumber=DAerrorNumber+1;
            else
                Rank(thatdam)=Rank(thisdam)+1;
            end
        end
    end
    i=i+1;
    ranknum=find(Rank==i);
end

if DAerrorNumber~=0
    error("Manually investigate DA errors in GIS and find solution")
end
data.Rank=Rank;

vars={'Rank','check','DAerrorNumber','goesto','i','imax','j','jmax','n','ranknum','thatdam','thatdamDA','thisdam','thisdamDA'};
clear(vars{:}); clear vars;

%% save data for mat 
 data1wrivers=data;

%% Now remove rivers from the "To Dam" category, update terminal and headwater flags

%remove rivers/deltas from the "To Dam" category. May need to loop through
%only once if downstream dams have their "ToDam" replaced before the
%upstream dam in encountered in the loop

while true
    % Find dams that go to a river
    numdamstoriver = find(data.todam < 0);
    if isempty(numdamstoriver)
        break; % Exit the loop if no such dams exist
    end
    
    % Update "to dam" values for these dams
    for damlocationtofix = numdamstoriver'
        theriver = data.todam(damlocationtofix); % SID of the river
        riverlocation = find(data.SID == theriver); % Location of the river
        data.todam(damlocationtofix) = data.todam(riverlocation); % Update to dam
    end
end

% identify "headwater" rivers (rivers that were flagged like headwater dams). If they go to a dam, then the dam then needs
% to be potentially reclassified as a headwater dam, if not downstream from another dam on antoher flowline. 
headrivers=find(data.isriver==1 & data.head==1);
fromHR=data.todam(headrivers);
realdams=find(fromHR>0);

if ~isempty(realdams)
    realdamlocation = find(ismember(data.SID, fromHR(realdams)));
    error("Check these dams: see if they need to be reclassified as headwater dams"); %if situation exists, create a solution
end

% now that to dam has been replaced on real dams, Remove "to dam" for rivers and update terminal dam flags
data.todam(rivers) = NaN;

% Update terminal dam flags
data.term(data.todam > 0) = 0; % Non-terminal dams
data.term(data.todam == 0) = 1; % Terminal dams
data.term(rivers) = NaN; % Rivers are not dams

clear realdams
%%  Make a ranking file of dam order WITHOUT RIVERS
%Make a ranking file ; Let the higher number win

Rank = NaN(nRows, 1);
Rank(data.head == 1) = 1;
Rank(rivers)=0;
DAerrorNumber = 0;

i=1;
ranknum=find(Rank==i);

while ~isempty(ranknum)
    jmax=numel(ranknum); % number of dams with a rank of i
    for j=1:jmax
        thisdam=ranknum(j); % The index we are on within array j
        
        %process nonterminal dams
        if  data.term(thisdam)==0 % if dam j is not a terminal dam, then
            thisdamDA=data.DA(thisdam); % identify DA of j
            goesto=data.todam(thisdam); % identify where dam j goes to
            thatdam=find(data.SID==goesto); % find the SID of the dam that dam j goes to
            thatdamDA=data.DA(thatdam); % identify the DA that dam j goes to

            %do a check to make sure downstream DA doesn't exceed upstream
            if (thatdamDA<thisdamDA && data.SID(thatdam)>0) %identify if upstream drainage area is greater than downstream drainage area, and not a river, b/c that will be fixed
                DAerrorNumber=DAerrorNumber+1;
            else
                Rank(thatdam)=Rank(thisdam)+1;
            end
        end
    end
    i=i+1;
    ranknum=find(Rank==i);
end

if DAerrorNumber~=0
    error("Manually investigate DA errors in GIS and find solution")
end
data.Rank=Rank;

vars={'Rank','check','DAerrorNumber','goesto','i','imax','j','jmax','n','ranknum','thatdam','thatdamDA','thisdam','thisdamDA'};
clear(vars{:}); clear vars;

%% Get Rid of Rivers from data

% Remove rivers from data
rivers = data(data.isriver == 1, :); % Extract rows where isriver is 1
data(data.isriver == 1, :) = []; % Remove those rows from data

%% save data for mat without rivers

data2norivers=data;
%save('data2.mat','data2norivers');

%% Print things 
disp(['Headwater Dams (excluding rivers): ', num2str(sum(data.head == 1 & data.isriver == 0))]);
disp(['Terminal Dams (excluding rivers): ', num2str(sum(data.term == 1 & data.isriver == 0))]);

%% Get Rid of dams that go nowhere, unless they are sites. get rid of dams without site tags for SY model
% Remove dams that go nowhere, unless they are sites
IsolatedDams = data(data.head == 1 & data.term == 1 & data.issite == 0, :); % Extract isolated dams - don't need in SY model
data(data.head == 1 & data.term == 1 & data.issite == 0, :) = []; % Remove them from the dataset


%% Setup: Also remove non-site dams without a site tag- because won't go to a site, not part of sed yield model. For MLR model only

% Identify and remove non-site rows with no site tags
NoTag = data(data.sitetag1 == 0 & data.issite == 0, :);
data(data.sitetag1 == 0 & data.issite == 0, :) = [];

%% save data3_SYmodel for mat, Save mat

data3SYmodel=data;
save('RankTag.mat','data1wrivers','data2norivers','data3SYmodel');

%% Print things
% Print the number of isolated dams
disp('*');
disp('Done with AddInputsUpdats Function');
disp(['Isolated Dams (no todam or fromdam): ', num2str(height(IsolatedDams))]);
disp(['Dams Outside Surveyed Basins: ', num2str(height(NoTag))]);
disp('*');
disp('Done with RankTag Function');


end




