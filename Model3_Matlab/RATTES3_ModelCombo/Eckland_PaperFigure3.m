function [fignum] = Eckland_PaperFigure3(data5all,calc5all,fignum,calc6_delta,data6_delta)

%% Define 
calc6=calc6_delta;
data6=data6_delta;
data=data5all;
calc=calc5all;

numt=numel(calc.t);
t1900=find(calc.t==1900);

slr=calc6.sealevelrise;
EquivDepth_dt=calc6.EquivDepth_dt;
main=calc6.maindeltas;
x = calc.t(t1900:numt).';
y1=EquivDepth_dt(main,t1900:numt);
y2 = slr(t1900:numt);

%% create figure
figure(fignum)
set(gcf,'InvertHardcopy','off', ...
        'Color',[1 1 1], ...
        'OuterPosition',[762.3333 389.6667 716.6667 692.6667]);

ax2 = axes('Position', ...
    [0.51160444 0.26747720 0.38133406 0.65752280]);

hold(ax2,'on')
set(ax2, ...
    'Color',[0.9 0.9 0.9], ...
    'FontName','Arial', ...
    'FontSize',9, ...
    'Layer','top', ...
    'XColor','k', ...
    'YColor','k', ...
    'XTick',[1900 1925 1950 1975 2000 2025 2050], ...
    'YScale','log', ...
    'YMinorTick','on');

semilogy(ax2, x, (y1.'*1e3), 'LineWidth',1.5)
colororder([0 0 0;0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;1 0 1;0 1 1;0.5 0.5 0.5;0.25 0.15 0.05;0 0.4471 0.7412;0.9 0.8 0.95;0.8667 0.7961 0.6431;0 0.8549 0.6902;0.1294 0.3412 0.1961]);

semilogy(ax2, x, y2.', 'r--','LineWidth',2)

xlabel(ax2,'Year')
ylabel(ax2,'\it h \rm (mm yr^{-1}), equivalent aggradation rate','Interpreter','tex')
xlim(ax2,[1900 2050])
ylim(ax2,[1e-2 600])
box(ax2,'on')
grid(ax2,'on')

lgd = legend(ax2, ...
    {'Mississippi','Colorado','Columbia','Savannah','Rio Grande','Mobile', ...
     'Trinity','Santee','Brazos','Elwha','Sabine','Colorado (Texas)', ...
     'Nueces','Apalachicola','Altamaha','Pearl','Sea level rise'});

set(lgd, ...
    'Position',[0.06797 0.10491 0.89863 0.10472], ...
    'NumColumns',5, ...
    'IconColumnWidth',35, ...
    'FontSize',8, ...
    'EdgeColor','none', ...
    'Color',[1 1 1]);

%% Figure 3A set up
%PathX

% Unique delta IDs
deltaid = data6.deltaID;
numdelta = numel(deltaid);

% Distance vector rounded to nearest 1 km
maxval = ceil(max(data.pathx));   % nearest 1 km
dist = 0:1:maxval;

% Preallocate arrays
pathxsed25 = NaN(numdelta, numel(dist));
topSedX = NaN(numdelta, 5);    % X positions of top raw sediment
topSedY = NaN(numdelta, 5);    % Cumulative sediment at top raw sediment positions

% Time indices
now = find(calc.t == 2025);
time = now + zeros(numdelta,1);

% Loop through each delta ID
for j = 1:numdelta
    currentID = deltaid(j);

    % Relevant indices
    idx = find(data.SID < 500000 & data.deltatag == currentID & data.PermStorag==1);

    if ~isempty(idx)
        % Sediment and X data
        sed.x = data.pathx(idx);       % distance from coast
        sed.y = calc.sed(idx, time(j)); % sediment volume at given time
        [sed.x, sortIdx] = sort(sed.x);
        sed.y = sed.y(sortIdx);

        % Cumulative sediment along the path
        cumSed = cumsum(sed.y);

        % Maximum distance for this delta, rounded up to nearest km
        maxSedX = ceil(max(sed.x));

        % Only loop over distances that exist in dist
        validK = dist <= maxSedX;

        % Compute cumulative sediment up to each distance
        for k = find(validK)
            pathxsed25(j, k) = sum(sed.y(sed.x <= dist(k)));
        end

        % Optional: fill remaining entries with NaN to avoid plotting beyond data
        pathxsed25(j, ~validK) = NaN;

        % Find top 5 raw sediment values
        [sortedSedY, sortIdx] = sort(sed.y, 'descend');
        nTop = min(5, numel(sortedSedY));

        % Store cumulative sediment at positions of top raw sediment
        topSedY(j, 1:nTop) = cumSed(sortIdx(1:nTop));
        topSedX(j, 1:nTop) = round(sed.x(sortIdx(1:nTop)));  % round to nearest 1 km
    end
end

%convert cumulative sediment values to m3 million m3
topSedY=topSedY/1e6; 
pathxsed25=pathxsed25/1e6;
%%

a=find(main==28); %delta ID Elwha == 28
main(a)=[]; %eliminate Elwha
newcolors=colororder;
newcolors(a,:)=[];

ax1 = axes('Position',...
    [0.13 0.269503546099291 0.294823410696266 0.655496453900711]);
hold(ax1,'on');

% Force log-log axes explicitly
ax1 = gca;
% Force log-log axes
set(ax1, 'XScale', 'log', 'YScale', 'log');
hold on
% Max values for axes
maxy=max(pathxsed25(main,end));
maxploty=ceil(maxy/10000)*10000;
maxplotx = 6000; %km

minploty=50;
minplotx=100;

hLines = gobjects(numel(main),1);
for i = 1:numel(main)
    val=main(i);
    hassed=find(pathxsed25(val,:)>0);
    hLines(i) = loglog(dist(hassed), pathxsed25(val,hassed), 'LineWidth', 1.5, 'Color',newcolors(i,:));
end

% pearl (15), nueces (21), have top 5 sed points with less than 1M m3 sediment.  don't plot those
topSedX(15,2:5)=NaN;
topSedY(15,2:5)=NaN;
topSedX(21,4:5)=NaN;
topSedY(21,4:5)=NaN;

% Plot squares using the same line colors
for i = 1:numel(main)
    val=main(i);
    hassedpt=find(topSedY(val,:)>0);
    loglog(topSedX(val,hassedpt), topSedY(val,hassedpt), ...
        's', 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', newcolors(i,:), ...
        'MarkerSize', 6, 'LineWidth', 1.25)
end

% Force log-log axes explicitly
ax1 = gca;
% Force log-log axes
set(ax1, 'XScale', 'log', 'YScale', 'log');

% Set custom tick values
ax1.XTick = [minplotx 500 1000 maxplotx];
ax1.YTick = [minploty 100 500 1000 5000 10000 maxploty];

% Optional: make tick labels look nice
ax1.XTickLabel = arrayfun(@num2str, ax1.XTick, 'UniformOutput', false);
ax1.YTickLabel = arrayfun(@num2str, ax1.YTick, 'UniformOutput', false);

% Axes properties
ax1.FontSize = 16;
ax1.Color = [0.9 0.9 0.9];
ax1.Box = 'on';
ax1.Layer = 'top';
ax1.XColor = 'k';
ax1.YColor = 'k';

xlabel('Distance from coast (km)')
ylabel('Cumulative Sediment Volume, Million m^3')

% Set limits safely above zero
xlim([minplotx maxplotx])
ylim([minploty maxploty])

box(ax1,'on');
grid(ax1,'on');
hold(ax1,'off');
% Set the remaining axes properties
set(ax1,'Color',[0.9 0.9 0.9],'FontName','Arial','FontSize',9,'GridAlpha',...
    0.5,'GridColor',[0.7 0.7 0.7],'Layer','top','XColor',[0 0 0],'XMinorTick',...
    'on','XScale','log','XTick',[100 300 500 1000 3000 6000],'XTickLabel',...
    {'100','300','500','1000','3000','6000'},'YColor',[0 0 0],'YMinorTick','on',...
    'YScale','log','YTick',[50 100 300 1000 3000 5000 10000 30000 40000],...
    'YTickLabel',...
    {'50','100','300','1000','3000','5000','10000','30000','40000'});
% Create textbox
annotation(gcf,'textbox',...
    [0.365464692482916 0.271529888551166 0.0605034168564921 0.0587639311043565],...
    'String',{'A'},...
    'FontSize',24,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.830726651480639 0.272036474164135 0.0605034168564921 0.0587639311043565],...
    'String','B',...
    'FontSize',24,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');

savefig(gcf,'Figure3Eckland.fig');
fignum=fignum+1;
end
