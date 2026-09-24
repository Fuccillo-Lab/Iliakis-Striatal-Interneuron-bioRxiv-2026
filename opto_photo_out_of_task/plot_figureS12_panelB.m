function results = plot_figureS12_panelB( ...
    streamKeyFile,h5File,outputFolder,varargin)
%PLOT_FIGURES12_PANELB Reproduce the optogenetic-inhibition photometry panel.
%
% results = plot_figureS12_panelB(streamKeyFile,h5File,outputFolder)
% results = plot_figureS12_panelB(...,Name,Value)
%
% The function selects same-cell eNpHR3.0/GCaMP recordings from the processed
% out-of-task opto-photometry dataset, aligns the stored Z-scored dF/F signal
% to stimulation onset, collapses epochs and bilateral streams within each
% animal, and plots animal-level mean +/- SEM for SST and TH interneurons.
%
% Required inputs
% ---------------
% streamKeyFile : stream_key_build_clean_02.csv
% h5File        : streams_build_clean_02.h5
% outputFolder  : destination for figures and provenance tables
%
% Name-value options
% ------------------
% PowerMw       : stimulation power to plot (default 10)
% SignalName    : H5 signal dataset (default dff_debleached_artcorr_z)
% EpochWindow   : seconds around stimulation onset (default [-2.5 5])
% BaselineWindow: baseline-centering window (default [-2 0])
% StimWindow    : displayed light interval (default [0 1])
% EpochDt       : interpolation interval in seconds (default 0.025)
% Aggregation   : "mean" or "median" within stream and animal (default median)

narginchk(3,inf);

streamKeyFile = requireFile(streamKeyFile,"Stream-key CSV");
h5File = requireFile(h5File,"Processed-stream H5");
outputFolder = string(outputFolder);
assert(isscalar(outputFolder) && strlength(outputFolder)>0, ...
    'outputFolder must be a nonempty path.');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

cfg = parseInputs(varargin{:});
epochTime = (cfg.EpochWindow(1):cfg.EpochDt:cfg.EpochWindow(2))';

T = readtable(streamKeyFile,'TextType','string');
requiredVariables = ["stream_id","animalID","siteSide","power_mW", ...
    "stimProtocol","usable_site","opsin","opto_cell","photo_cell", ...
    "import_status"];
requireVariables(T,requiredVariables,"Stream-key CSV");

T.animalID = numericColumn(T.animalID);
T.power_mW = numericColumn(T.power_mW);

optoCell = normalizeCellName(T.opto_cell);
photoCell = normalizeCellName(T.photo_cell);
opsin = lower(strtrim(string(T.opsin)));
protocol = lower(strtrim(string(T.stimProtocol)));
importStatus = lower(strtrim(string(T.import_status)));
usable = logicalColumn(T.usable_site);

T.cohort = repmat("",height(T),1);
sameCellHalo = contains(opsin,"halo") & optoCell==photoCell;
T.cohort(sameCellHalo & optoCell=="sst") = "SST";
T.cohort(sameCellHalo & optoCell=="th") = "TH";

selected = importStatus=="ok" & usable & ...
    ismember(protocol,["contin","continuous"]) & ...
    abs(T.power_mW-cfg.PowerMw)<1e-9 & strlength(T.cohort)>0;
T = T(selected,:);

assert(~isempty(T), ...
    'No usable same-cell Halo streams matched power %g mW.',cfg.PowerMw);
T = sortrows(T,{'cohort','animalID','siteSide'});

fprintf('Figure S12B: selected %d streams at %g mW.\n', ...
    height(T),cfg.PowerMw);

%% Extract and epoch every selected stream

nStreams = height(T);
traceCells = cell(nStreams,1);
nTTL = zeros(nStreams,1);
nValidEpochs = zeros(nStreams,1);
included = false(nStreams,1);
exclusionReason = repmat("",nStreams,1);

for i = 1:nStreams
    streamID = string(T.stream_id(i));
    groupPath = "/streams/"+streamID;

    try
        timestamps = double(h5read( ...
            char(h5File),char(groupPath+"/timestamps")));
        ttlTimestamps = double(h5read( ...
            char(h5File),char(groupPath+"/ttl_timestamps")));
        signal = double(h5read( ...
            char(h5File),char(groupPath+"/"+cfg.SignalName)));

        ttlTimestamps = ttlTimestamps(isfinite(ttlTimestamps));
        nTTL(i) = numel(ttlTimestamps);
        [streamTrace,nValidEpochs(i)] = epochAndAggregate( ...
            timestamps,signal,ttlTimestamps,epochTime,cfg);

        if nValidEpochs(i)==0
            exclusionReason(i) = "no_complete_epochs";
        else
            traceCells{i} = streamTrace;
            included(i) = true;
        end
    catch ME
        exclusionReason(i) = "read_or_epoch_failed: "+string(ME.message);
        warning('Skipping stream %s: %s',streamID,ME.message);
    end
end

streamSummary = T(:,{'stream_id','animalID','cohort','siteSide','power_mW'});
streamSummary.signal = repmat(cfg.SignalName,nStreams,1);
streamSummary.nTTL = nTTL;
streamSummary.nValidEpochs = nValidEpochs;
streamSummary.included = included;
streamSummary.exclusionReason = exclusionReason;

assert(any(included),'No selected stream contained a complete valid epoch.');

%% Collapse streams within animal

animalCohort = strings(0,1);
animalID = zeros(0,1);
animalNStreams = zeros(0,1);
animalNEpochs = zeros(0,1);
animalTraces = zeros(numel(epochTime),0);

cohortOrder = ["SST","TH"];
for c = 1:numel(cohortOrder)
    cohort = cohortOrder(c);
    animalList = unique(T.animalID(included & T.cohort==cohort),'stable');

    for a = 1:numel(animalList)
        useStream = find(included & T.cohort==cohort & ...
            T.animalID==animalList(a));
        traces = zeros(numel(epochTime),numel(useStream));
        for j = 1:numel(useStream)
            traces(:,j) = traceCells{useStream(j)};
        end

        animalCohort(end+1,1) = cohort; %#ok<AGROW>
        animalID(end+1,1) = animalList(a); %#ok<AGROW>
        animalNStreams(end+1,1) = numel(useStream); %#ok<AGROW>
        animalNEpochs(end+1,1) = sum(nValidEpochs(useStream)); %#ok<AGROW>
        animalTraces(:,end+1) = aggregateMatrix( ...
            traces,cfg.Aggregation,2); %#ok<AGROW>
    end
end

nAnimalsTotal = numel(animalID);
assert(nAnimalsTotal>0,'No animal traces could be constructed.');

animalTraceTable = table( ...
    repelem(animalCohort,numel(epochTime)), ...
    repelem(animalID,numel(epochTime)), ...
    repelem(animalNStreams,numel(epochTime)), ...
    repelem(animalNEpochs,numel(epochTime)), ...
    repmat(epochTime,nAnimalsTotal,1), ...
    animalTraces(:), ...
    'VariableNames',{'cohort','animalID','nStreams','nEpochs', ...
    'timeFromStimSec','zScoredDff'});

%% Compute animal-weighted group summaries

summaryCohort = strings(0,1);
summaryTime = zeros(0,1);
summaryMean = zeros(0,1);
summarySem = zeros(0,1);
summaryN = zeros(0,1);

for c = 1:numel(cohortOrder)
    cohort = cohortOrder(c);
    useAnimal = animalCohort==cohort;
    assert(any(useAnimal),'No included animals remain for the %s cohort.',cohort);

    traces = animalTraces(:,useAnimal);
    nAnimals = size(traces,2);
    meanTrace = mean(traces,2,'omitnan');
    semTrace = std(traces,0,2,'omitnan')./sqrt(nAnimals);

    summaryCohort = [summaryCohort; ...
        repmat(cohort,numel(epochTime),1)]; %#ok<AGROW>
    summaryTime = [summaryTime;epochTime]; %#ok<AGROW>
    summaryMean = [summaryMean;meanTrace]; %#ok<AGROW>
    summarySem = [summarySem;semTrace]; %#ok<AGROW>
    summaryN = [summaryN;repmat(nAnimals,numel(epochTime),1)]; %#ok<AGROW>

    fprintf('Figure S12B %s: N=%d animals (%s).\n', ...
        cohort,nAnimals,strjoin(string(animalID(useAnimal)),', '));
end

groupSummary = table( ...
    summaryCohort,summaryTime,summaryMean,summarySem,summaryN, ...
    'VariableNames',{'cohort','timeFromStimSec','meanZScoredDff', ...
    'semZScoredDff','nAnimals'});

%% Plot Figure S12B

fig = figure('Color','w','Units','inches','Position',[1 1 2.45 4.10]);

sstColor = [37 61 143]./255;
thColor = [152 56 148]./255;
titleColor = [0.95 0.20 0.16];
lightColor = [1.00 0.68 0.70];

axSst = subplot(2,1,1,'Parent',fig);
plotCohort(axSst,"SST",sstColor,titleColor,lightColor, ...
    epochTime,animalCohort,animalTraces,cfg,false);

axTh = subplot(2,1,2,'Parent',fig);
plotCohort(axTh,"TH",thColor,titleColor,lightColor, ...
    epochTime,animalCohort,animalTraces,cfg,true);

set(axSst,'Position',[0.22 0.56 0.72 0.35]);
set(axTh,'Position',[0.22 0.13 0.72 0.35]);

annotation(fig,'textbox',[0.005 0.31 0.055 0.40], ...
    'String','\DeltaF/F (Z-scored)', ...
    'Interpreter','tex','Rotation',90, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Arial','FontSize',9,'EdgeColor','none');

%% Export figure and provenance tables

outputStem = fullfile(outputFolder,"figureS12_panelB_opto_photometry");
exportFigurePair(fig,outputStem);
writetable(animalTraceTable,outputStem+"_animal_traces.csv");
writetable(groupSummary,outputStem+"_group_summary.csv");
writetable(streamSummary,outputStem+"_stream_summary.csv");

results = struct();
results.figure = fig;
results.animalTraces = animalTraceTable;
results.groupSummary = groupSummary;
results.streamSummary = streamSummary;
results.parameters = cfg;
results.inputFiles = struct( ...
    'streamKey',streamKeyFile,'processedStreams',h5File);
results.outputStem = outputStem;

end


function cfg = parseInputs(varargin)
p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'PowerMw',5,@isFiniteScalar);
addParameter(p,'SignalName',"dff_debleached_artcorr_z",@isTextScalar);
addParameter(p,'EpochWindow',[-2.5 5],@isTwoElementWindow);
addParameter(p,'BaselineWindow',[-2 0],@isTwoElementWindow);
addParameter(p,'StimWindow',[0 1],@isTwoElementWindow);
addParameter(p,'EpochDt',0.025,@(x)isFiniteScalar(x) && x>0);
addParameter(p,'Aggregation',"median", ...
    @(x)any(strcmpi(string(x),["mean","median"])));
parse(p,varargin{:});
cfg = p.Results;
cfg.PowerMw = double(cfg.PowerMw);
cfg.SignalName = string(cfg.SignalName);
cfg.EpochWindow = double(cfg.EpochWindow(:)');
cfg.BaselineWindow = double(cfg.BaselineWindow(:)');
cfg.StimWindow = double(cfg.StimWindow(:)');
cfg.EpochDt = double(cfg.EpochDt);
cfg.Aggregation = lower(string(cfg.Aggregation));

assert(cfg.BaselineWindow(1)>=cfg.EpochWindow(1) && ...
    cfg.BaselineWindow(2)<=cfg.EpochWindow(2), ...
    'BaselineWindow must fall within EpochWindow.');
assert(cfg.StimWindow(1)>=cfg.EpochWindow(1) && ...
    cfg.StimWindow(2)<=cfg.EpochWindow(2), ...
    'StimWindow must fall within EpochWindow.');
end


function [trace,nValid] = epochAndAggregate( ...
    timestamps,signal,ttlTimestamps,epochTime,cfg)

timestamps = timestamps(:);
signal = signal(:);
assert(numel(timestamps)==numel(signal), ...
    'Signal and timestamp lengths differ.');

validSample = isfinite(timestamps) & isfinite(signal);
timestamps = timestamps(validSample);
signal = signal(validSample);
[timestamps,order] = sort(timestamps);
signal = signal(order);
[timestamps,uniqueIndex] = unique(timestamps,'stable');
signal = signal(uniqueIndex);

ttlTimestamps = ttlTimestamps(:);
ttlTimestamps = ttlTimestamps(isfinite(ttlTimestamps));
epochs = zeros(numel(epochTime),0);

for k = 1:numel(ttlTimestamps)
    queryTime = ttlTimestamps(k)+epochTime;
    if queryTime(1)<timestamps(1) || queryTime(end)>timestamps(end)
        continue;
    end
    epoch = interp1(timestamps,signal,queryTime,'linear',NaN);
    if all(isfinite(epoch))
        epochs(:,end+1) = epoch; %#ok<AGROW>
    end
end

nValid = size(epochs,2);
if nValid==0
    trace = [];
    return;
end

baselineIndex = epochTime>=cfg.BaselineWindow(1) & ...
    epochTime<cfg.BaselineWindow(2);
assert(any(baselineIndex),'BaselineWindow contains no epoch samples.');
baseline = aggregateMatrix(epochs(baselineIndex,:),cfg.Aggregation,1);
epochs = bsxfun(@minus,epochs,baseline);
trace = aggregateMatrix(epochs,cfg.Aggregation,2);
end


function plotCohort(ax,cohort,lineColor,titleColor,lightColor, ...
    epochTime,animalCohort,animalTraces,cfg,showXAxis)

useAnimal = animalCohort==cohort;
traces = animalTraces(:,useAnimal);
nAnimals = size(traces,2);
meanTrace = mean(traces,2,'omitnan');
semTrace = std(traces,0,2,'omitnan')./sqrt(nAnimals);

hold(ax,'on');
plot(ax,[cfg.EpochWindow(1) cfg.EpochWindow(2)],[0 0],'-', ...
    'Color',[0.82 0.82 0.82],'LineWidth',0.6);
fill(ax,[epochTime;flipud(epochTime)], ...
    [meanTrace-semTrace;flipud(meanTrace+semTrace)], ...
    lineColor,'FaceAlpha',0.18,'EdgeColor','none');
plot(ax,epochTime,meanTrace,'Color',lineColor,'LineWidth',1.35);

lowerEnvelope = meanTrace-semTrace;
upperEnvelope = meanTrace+semTrace;
dataLower = min(lowerEnvelope(isfinite(lowerEnvelope)));
dataUpper = max(upperEnvelope(isfinite(upperEnvelope)));
dataRange = dataUpper-dataLower;
if ~isfinite(dataRange) || dataRange<=0
    dataRange = 1;
end
yLower = dataLower-0.08*dataRange;
yUpper = dataUpper+0.24*dataRange;
ylim(ax,[yLower yUpper]);
xlim(ax,cfg.EpochWindow);
xticks(ax,[-2 0 2 4]);

barBottom = dataUpper+0.08*dataRange;
barTop = dataUpper+0.13*dataRange;
patch(ax, ...
    [cfg.StimWindow(1) cfg.StimWindow(2) cfg.StimWindow(2) cfg.StimWindow(1)], ...
    [barBottom barBottom barTop barTop], ...
    lightColor,'EdgeColor','none');
text(ax,mean(cfg.StimWindow),barTop+0.015*dataRange,'Light ON', ...
    'Color',titleColor,'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontName','Arial','FontSize',8);

if cohort=="SST"
    titleText = 'SST-eNpHR3.0';
    traceLabel = 'SST-GCaMP';
else
    titleText = 'TH-eNpHR3.0';
    traceLabel = 'TH-GCaMP';
end
title(ax,titleText,'Color',titleColor,'FontName','Arial', ...
    'FontSize',9,'FontWeight','normal','HorizontalAlignment','left');

labelText = sprintf('%s\nN = %d animals',traceLabel,nAnimals);
text(ax,cfg.EpochWindow(2)-0.05*diff(cfg.EpochWindow), ...
    yLower+0.08*(yUpper-yLower),labelText, ...
    'Color',lineColor,'Interpreter','none', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontName','Arial','FontSize',8);

set(ax,'Box','off','TickDir','out','FontName','Arial','FontSize',8, ...
    'LineWidth',0.75,'Layer','top');
if showXAxis
    xlabel(ax,'Time from stim (s)','FontName','Arial','FontSize',9);
else
    set(ax,'XTickLabel',[]);
end
end


function out = aggregateMatrix(X,aggregation,dimension)
switch aggregation
    case "mean"
        out = mean(X,dimension,'omitnan');
    case "median"
        out = median(X,dimension,'omitnan');
    otherwise
        error('Unknown aggregation: %s',aggregation);
end
end


function values = numericColumn(values)
if isnumeric(values) || islogical(values)
    values = double(values);
else
    values = str2double(string(values));
end
values = values(:);
end


function values = logicalColumn(values)
if islogical(values)
    values = values(:);
elseif isnumeric(values)
    values = isfinite(values) & values~=0;
else
    values = ismember(lower(strtrim(string(values))), ...
        ["true","1","y","yes"]);
end
values = values(:);
end


function names = normalizeCellName(values)
names = lower(strtrim(string(values)));
names(ismember(names,["ltsi","sst"])) = "sst";
names(ismember(names,["thin","th"])) = "th";
end


function requireVariables(T,variables,label)
missing = setdiff(variables,string(T.Properties.VariableNames));
assert(isempty(missing),'%s is missing required variable(s): %s', ...
    label,strjoin(missing,', '));
end


function file = requireFile(file,label)
file = string(file);
assert(isscalar(file) && strlength(file)>0, ...
    '%s path must be a nonempty text scalar.',label);
assert(isfile(file),'%s not found: %s',label,file);
end


function exportFigurePair(fig,outputStem)
pdfFile = outputStem+".pdf";
pngFile = outputStem+".png";

if exist('exportgraphics','file')==2
    exportgraphics(fig,pdfFile,'ContentType','vector');
    exportgraphics(fig,pngFile,'Resolution',600);
else
    oldPaperUnits = get(fig,'PaperUnits');
    oldPaperSize = get(fig,'PaperSize');
    oldPaperPosition = get(fig,'PaperPosition');
    oldPaperPositionMode = get(fig,'PaperPositionMode');
    oldUnits = get(fig,'Units');

    set(fig,'Units','inches');
    figurePosition = get(fig,'Position');
    figureSize = figurePosition(3:4);
    set(fig,'PaperUnits','inches', ...
        'PaperSize',figureSize, ...
        'PaperPosition',[0 0 figureSize], ...
        'PaperPositionMode','manual');
    print(fig,char(pdfFile),'-dpdf','-painters');
    print(fig,char(pngFile),'-dpng','-r600');

    set(fig,'Units',oldUnits, ...
        'PaperUnits',oldPaperUnits, ...
        'PaperSize',oldPaperSize, ...
        'PaperPosition',oldPaperPosition, ...
        'PaperPositionMode',oldPaperPositionMode);
end

fprintf('Saved Figure S12B PDF: %s\n',pdfFile);
fprintf('Saved Figure S12B PNG: %s\n',pngFile);
end


function tf = isTextScalar(x)
tf = (ischar(x) && (isrow(x) || isempty(x))) || ...
    (isstring(x) && isscalar(x));
end


function tf = isFiniteScalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x);
end


function tf = isTwoElementWindow(x)
tf = isnumeric(x) && numel(x)==2 && all(isfinite(x)) && x(2)>x(1);
end
