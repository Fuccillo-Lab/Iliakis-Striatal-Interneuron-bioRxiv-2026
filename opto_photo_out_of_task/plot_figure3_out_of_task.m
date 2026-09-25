function results = plot_figure3_out_of_task( ...
    streamKeyFile,h5File,outputFolder,varargin)
%PLOT_FIGURE3_OUT_OF_TASK Reproduce Figure 3E-G out-of-task photometry.
%
% results = plot_figure3_out_of_task(streamKeyFile,h5File,outputFolder)
% results = plot_figure3_out_of_task(...,Name,Value)
%
% The function compares 5-mW TH-ChrimsonR stimulation with opsin-negative
% controls while recording SST, D1, or A2A populations. Stimulation epochs
% and bilateral streams are collapsed within animal before group plotting
% and cluster-permutation statistics.
%
% Required inputs
% ---------------
% streamKeyFile : stream_key_build_clean_02.csv
% h5File        : streams_build_clean_02.h5
% outputFolder  : destination for figures and provenance tables
%
% Name-value options
% ------------------
% PowerMw          : stimulation power (default 5)
% SignalName       : H5 signal (default dff_debleached_artcorr)
% EpochWindow      : extracted seconds around onset (default [-2.5 5])
% PlotWindow       : displayed seconds around onset (default [-1 3])
% BaselineWindow   : baseline-centering window (default [-2 0])
% StimWindow       : displayed stimulation interval (default [0 1])
% TestWindow       : cluster-test interval (default [0 1.5])
% EpochDt          : interpolation interval in seconds (default 0.025)
% Aggregation      : "mean" or "median" (default mean)
% RunClusterTest   : run excitation-control permutation tests (default true)
% NPermutations    : number of label permutations (default 5000)
% RandomSeed       : deterministic permutation seed (default 1)
% ClusterFormingAlpha : pointwise cluster threshold (default 0.05)
% ClusterAlpha     : corrected cluster threshold (default 0.05)

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
status = lower(strtrim(string(T.import_status)));
usable = logicalColumn(T.usable_site);

isExcitation = contains(opsin,"chrimson") & optoCell=="th" & ...
    protocol=="pulsed";
isControl = controlMask(T,opsin) & ...
    ismember(protocol,["pulsed","unspecified","unspec",""]);

baseSelection = status=="ok" & usable & ...
    abs(T.power_mW-cfg.PowerMw)<1e-9 & ...
    ismember(photoCell,["sst","d1","a2a"]);

T.cellType = upper(photoCell);
T.group = repmat("",height(T),1);
T.group(isExcitation) = "excitation";
T.group(isControl) = "control";
T = T(baseSelection & strlength(T.group)>0,:);

assert(~isempty(T), ...
    'No Figure 3 streams matched the requested %g-mW condition.', ...
    cfg.PowerMw);
T = sortrows(T,{'cellType','group','animalID','siteSide'});

cellOrder = ["SST","D1","A2A"];
groupOrder = ["control","excitation"];

fprintf('Figure 3E-G: selected %d streams at %g mW.\n', ...
    height(T),cfg.PowerMw);
for c = 1:numel(cellOrder)
    for g = 1:numel(groupOrder)
        use = T.cellType==cellOrder(c) & T.group==groupOrder(g);
        fprintf('  %s %s: %d streams, %d animals.\n', ...
            cellOrder(c),groupOrder(g),sum(use),numel(unique(T.animalID(use))));
    end
end

%% Extract and aggregate each selected stream

nStreams = height(T);
traceCells = cell(nStreams,1);
responseAmplitude = nan(nStreams,1);
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
            responseIndex = epochTime>=cfg.StimWindow(1) & ...
                epochTime<=cfg.StimWindow(2);
            responseAmplitude(i) = aggregateMatrix( ...
                streamTrace(responseIndex),cfg.Aggregation,1);
            included(i) = true;
        end
    catch ME
        exclusionReason(i) = "read_or_epoch_failed: "+string(ME.message);
        warning('Skipping stream %s: %s',streamID,ME.message);
    end
end

streamSummary = T(:,{'stream_id','animalID','cellType','group', ...
    'siteSide','power_mW','opto_cell','photo_cell','opsin'});
streamSummary.signal = repmat(cfg.SignalName,nStreams,1);
streamSummary.nTTL = nTTL;
streamSummary.nValidEpochs = nValidEpochs;
streamSummary.responseAmplitude = responseAmplitude;
streamSummary.included = included;
streamSummary.exclusionReason = exclusionReason;

assert(any(included),'No selected stream contained a complete valid epoch.');

%% Collapse bilateral streams within animal

animalCellType = strings(0,1);
animalGroup = strings(0,1);
animalID = zeros(0,1);
animalNStreams = zeros(0,1);
animalNEpochs = zeros(0,1);
animalResponse = zeros(0,1);
animalTraces = zeros(numel(epochTime),0);

for c = 1:numel(cellOrder)
    for g = 1:numel(groupOrder)
        cellType = cellOrder(c);
        group = groupOrder(g);
        animalList = unique(T.animalID(included & ...
            T.cellType==cellType & T.group==group),'stable');

        for a = 1:numel(animalList)
            useStream = find(included & T.cellType==cellType & ...
                T.group==group & T.animalID==animalList(a));
            traces = zeros(numel(epochTime),numel(useStream));
            for j = 1:numel(useStream)
                traces(:,j) = traceCells{useStream(j)};
            end

            animalCellType(end+1,1) = cellType; %#ok<AGROW>
            animalGroup(end+1,1) = group; %#ok<AGROW>
            animalID(end+1,1) = animalList(a); %#ok<AGROW>
            animalNStreams(end+1,1) = numel(useStream); %#ok<AGROW>
            animalNEpochs(end+1,1) = sum(nValidEpochs(useStream)); %#ok<AGROW>
            animalResponse(end+1,1) = aggregateMatrix( ...
                responseAmplitude(useStream),cfg.Aggregation,1); %#ok<AGROW>
            animalTraces(:,end+1) = aggregateMatrix( ...
                traces,cfg.Aggregation,2); %#ok<AGROW>
        end
    end
end

nAnimalsTotal = numel(animalID);
assert(nAnimalsTotal>0,'No animal traces could be constructed.');

animalTraceTable = table( ...
    repelem(animalCellType,numel(epochTime)), ...
    repelem(animalGroup,numel(epochTime)), ...
    repelem(animalID,numel(epochTime)), ...
    repelem(animalNStreams,numel(epochTime)), ...
    repelem(animalNEpochs,numel(epochTime)), ...
    repmat(epochTime,nAnimalsTotal,1), ...
    animalTraces(:), ...
    'VariableNames',{'cellType','group','animalID','nStreams', ...
    'nEpochs','timeFromStimSec','dff'});

animalResponseTable = table(animalCellType,animalGroup,animalID, ...
    animalNStreams,animalNEpochs,animalResponse, ...
    'VariableNames',{'cellType','group','animalID','nStreams', ...
    'nEpochs','responseAmplitude'});

%% Group summaries and cluster-permutation tests

summaryCellType = strings(0,1);
summaryGroup = strings(0,1);
summaryTime = zeros(0,1);
summaryMean = zeros(0,1);
summarySem = zeros(0,1);
summaryN = zeros(0,1);

for c = 1:numel(cellOrder)
    for g = 1:numel(groupOrder)
        useAnimal = animalCellType==cellOrder(c) & ...
            animalGroup==groupOrder(g);
        assert(any(useAnimal), ...
            'No included animals remain for %s %s.', ...
            cellOrder(c),groupOrder(g));

        traces = animalTraces(:,useAnimal);
        nAnimals = size(traces,2);
        meanTrace = mean(traces,2,'omitnan');
        semTrace = std(traces,0,2,'omitnan')./sqrt(nAnimals);

        summaryCellType = [summaryCellType; ...
            repmat(cellOrder(c),numel(epochTime),1)]; %#ok<AGROW>
        summaryGroup = [summaryGroup; ...
            repmat(groupOrder(g),numel(epochTime),1)]; %#ok<AGROW>
        summaryTime = [summaryTime;epochTime]; %#ok<AGROW>
        summaryMean = [summaryMean;meanTrace]; %#ok<AGROW>
        summarySem = [summarySem;semTrace]; %#ok<AGROW>
        summaryN = [summaryN;repmat(nAnimals,numel(epochTime),1)]; %#ok<AGROW>

        fprintf('Figure 3 %s %s: N=%d animals (%s).\n', ...
            cellOrder(c),groupOrder(g),nAnimals, ...
            strjoin(string(animalID(useAnimal)),', '));
    end
end

groupSummary = table(summaryCellType,summaryGroup,summaryTime, ...
    summaryMean,summarySem,summaryN, ...
    'VariableNames',{'cellType','group','timeFromStimSec', ...
    'meanDff','semDff','nAnimals'});

clusterTable = emptyClusterTable();
clusterDetails = cell(numel(cellOrder),1);
significantClusters = cell(numel(cellOrder),1);

if cfg.RunClusterTest
    previousRng = rng;
    rngCleanup = onCleanup(@()rng(previousRng)); %#ok<NASGU>
    rng(cfg.RandomSeed,'twister');

    for c = 1:numel(cellOrder)
        useExcitation = animalCellType==cellOrder(c) & ...
            animalGroup=="excitation";
        useControl = animalCellType==cellOrder(c) & ...
            animalGroup=="control";

        [clusterDetails{c},cellClusterTable,significantClusters{c}] = ...
            clusterPermutationTest( ...
            animalTraces(:,useExcitation)', ...
            animalTraces(:,useControl)',epochTime,cfg,cellOrder(c));
        clusterTable = [clusterTable;cellClusterTable]; %#ok<AGROW>
    end
else
    for c = 1:numel(cellOrder)
        clusterDetails{c} = struct();
        significantClusters{c} = emptySignificantClusters();
    end
end

%% Plot Figure 3E-G

fig = figure('Color','w','Units','inches','Position',[1 1 7.20 2.70]);
panelColors = [37 61 143; 113 190 74; 236 28 36]./255;
yLimits = [-1.5 1.0; -0.4 0.8; -0.4 0.8];

axesHandles = gobjects(1,numel(cellOrder));
for c = 1:numel(cellOrder)
    axesHandles(c) = subplot(1,3,c,'Parent',fig);
    plotCellPanel(axesHandles(c),cellOrder(c),panelColors(c,:), ...
        yLimits(c,:),epochTime,animalCellType,animalGroup,animalTraces, ...
        significantClusters{c},cfg,c==1);
    text(axesHandles(c),-0.18,1.08,char(double('D')+c), ...
        'Units','normalized','FontName','Arial','FontSize',11, ...
        'FontWeight','bold','Color','k','Clipping','off');
end

set(axesHandles(1),'Position',[0.085 0.22 0.255 0.63]);
set(axesHandles(2),'Position',[0.390 0.22 0.255 0.63]);
set(axesHandles(3),'Position',[0.695 0.22 0.255 0.63]);

annotation(fig,'textbox',[0.28 0.90 0.44 0.07], ...
    'String','TH Connectivity out of task', ...
    'Interpreter','none','HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontName','Arial','FontSize',11, ...
    'FontWeight','bold','EdgeColor','none');

%% Export figure and provenance tables

outputStem = fullfile(outputFolder,"figure3_out_of_task_photometry");
exportFigurePair(fig,outputStem);
writetable(animalTraceTable,outputStem+"_animal_traces.csv");
writetable(animalResponseTable,outputStem+"_animal_responses.csv");
writetable(groupSummary,outputStem+"_group_summary.csv");
writetable(streamSummary,outputStem+"_stream_summary.csv");
if cfg.RunClusterTest
    writetable(clusterTable,outputStem+"_cluster_statistics.csv");
end

results = struct();
results.figure = fig;
results.animalTraces = animalTraceTable;
results.animalResponses = animalResponseTable;
results.groupSummary = groupSummary;
results.streamSummary = streamSummary;
results.clusterStatistics = clusterTable;
results.clusterDetails = clusterDetails;
results.parameters = cfg;
results.inputFiles = struct( ...
    'streamKey',streamKeyFile,'processedStreams',h5File);
results.outputStem = outputStem;

end


function cfg = parseInputs(varargin)
p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'PowerMw',5,@isFiniteScalar);
addParameter(p,'SignalName',"dff_debleached_artcorr",@isTextScalar);
addParameter(p,'EpochWindow',[-2.5 5],@isTwoElementWindow);
addParameter(p,'PlotWindow',[-1 3],@isTwoElementWindow);
addParameter(p,'BaselineWindow',[-2 0],@isTwoElementWindow);
addParameter(p,'StimWindow',[0 1],@isTwoElementWindow);
addParameter(p,'TestWindow',[0 1.5],@isTwoElementWindow);
addParameter(p,'EpochDt',0.025,@(x)isFiniteScalar(x) && x>0);
addParameter(p,'Aggregation',"mean", ...
    @(x)any(strcmpi(string(x),["mean","median"])));
addParameter(p,'RunClusterTest',true,@isLogicalScalar);
addParameter(p,'NPermutations',5000,@isPositiveInteger);
addParameter(p,'RandomSeed',1,@isNonnegativeInteger);
addParameter(p,'ClusterFormingAlpha',0.05,@isProbability);
addParameter(p,'ClusterAlpha',0.05,@isProbability);
parse(p,varargin{:});

cfg = p.Results;
cfg.PowerMw = double(cfg.PowerMw);
cfg.SignalName = string(cfg.SignalName);
cfg.EpochWindow = double(cfg.EpochWindow(:)');
cfg.PlotWindow = double(cfg.PlotWindow(:)');
cfg.BaselineWindow = double(cfg.BaselineWindow(:)');
cfg.StimWindow = double(cfg.StimWindow(:)');
cfg.TestWindow = double(cfg.TestWindow(:)');
cfg.EpochDt = double(cfg.EpochDt);
cfg.Aggregation = lower(string(cfg.Aggregation));
cfg.RunClusterTest = logical(cfg.RunClusterTest);
cfg.NPermutations = double(cfg.NPermutations);
cfg.RandomSeed = double(cfg.RandomSeed);

assert(windowWithin(cfg.BaselineWindow,cfg.EpochWindow), ...
    'BaselineWindow must fall within EpochWindow.');
assert(windowWithin(cfg.StimWindow,cfg.EpochWindow), ...
    'StimWindow must fall within EpochWindow.');
assert(windowWithin(cfg.TestWindow,cfg.EpochWindow), ...
    'TestWindow must fall within EpochWindow.');
assert(windowWithin(cfg.PlotWindow,cfg.EpochWindow), ...
    'PlotWindow must fall within EpochWindow.');
end


function mask = controlMask(T,opsin)
mask = false(height(T),1);
candidateVariables = ["opsin","genotype","virus_opsin","stimProtocol"];

for i = 1:numel(candidateVariables)
    variable = candidateVariables(i);
    if ~ismember(variable,string(T.Properties.VariableNames))
        continue;
    end
    values = lower(strtrim(string(T.(char(variable)))));
    mask = mask | contains(values,"ctrl") | ...
        contains(values,"control") | contains(values,"fluor") | ...
        contains(values,"gfp") | contains(values,"no opsin") | ...
        values=="none" | values=="" | values=="nan";
end

mask(contains(opsin,"chrimson") | contains(opsin,"halo")) = false;
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


function [stats,clusterTable,significantClusters] = ...
    clusterPermutationTest(XA,XB,time,cfg,cellType)

testIndex = time>=cfg.TestWindow(1) & time<=cfg.TestWindow(2);
testTime = time(testIndex);
XA = XA(:,testIndex);
XB = XB(:,testIndex);

XA = XA(all(isfinite(XA),2),:);
XB = XB(all(isfinite(XB),2),:);
nA = size(XA,1);
nB = size(XB,1);
assert(nA>=2 && nB>=2, ...
    'Cluster test for %s requires at least two animals per group.',cellType);

[tObserved,pObserved] = pooledTwoSampleT(XA,XB);
degreesFreedom = nA+nB-2;
tThreshold = studentTInverse( ...
    1-cfg.ClusterFormingAlpha/2,degreesFreedom);

positiveObserved = findClusters1D(tObserved>tThreshold);
negativeObserved = findClusters1D(tObserved<-tThreshold);
positiveMass = clusterMass(positiveObserved,tObserved,"positive");
negativeMass = clusterMass(negativeObserved,tObserved,"negative");

combined = [XA;XB];
nTotal = size(combined,1);
maxPositiveNull = zeros(cfg.NPermutations,1);
maxNegativeNull = zeros(cfg.NPermutations,1);

for permutation = 1:cfg.NPermutations
    order = randperm(nTotal);
    permutedA = combined(order(1:nA),:);
    permutedB = combined(order(nA+1:end),:);
    tPermuted = pooledTwoSampleT(permutedA,permutedB);

    positivePermuted = findClusters1D(tPermuted>tThreshold);
    negativePermuted = findClusters1D(tPermuted<-tThreshold);
    positivePermutedMass = clusterMass( ...
        positivePermuted,tPermuted,"positive");
    negativePermutedMass = clusterMass( ...
        negativePermuted,tPermuted,"negative");

    if ~isempty(positivePermutedMass)
        maxPositiveNull(permutation) = max(positivePermutedMass);
    end
    if ~isempty(negativePermutedMass)
        maxNegativeNull(permutation) = max(negativePermutedMass);
    end
end

positiveP = correctedClusterP(positiveMass,maxPositiveNull);
negativeP = correctedClusterP(negativeMass,maxNegativeNull);

clusterTable = emptyClusterTable();
significantClusters = emptySignificantClusters();

for i = 1:numel(positiveObserved)
    indices = positiveObserved{i};
    isSignificant = positiveP(i)<cfg.ClusterAlpha;
    clusterTable = [clusterTable;table(cellType,"positive", ...
        testTime(indices(1)),testTime(indices(end)),positiveMass(i), ...
        positiveP(i),isSignificant, ...
        'VariableNames',clusterTable.Properties.VariableNames)]; %#ok<AGROW>
    if isSignificant
        significantClusters(end+1) = makeCluster( ...
            testTime(indices(1)),testTime(indices(end)), ...
            positiveMass(i),positiveP(i),"positive"); %#ok<AGROW>
    end
end

for i = 1:numel(negativeObserved)
    indices = negativeObserved{i};
    isSignificant = negativeP(i)<cfg.ClusterAlpha;
    clusterTable = [clusterTable;table(cellType,"negative", ...
        testTime(indices(1)),testTime(indices(end)),negativeMass(i), ...
        negativeP(i),isSignificant, ...
        'VariableNames',clusterTable.Properties.VariableNames)]; %#ok<AGROW>
    if isSignificant
        significantClusters(end+1) = makeCluster( ...
            testTime(indices(1)),testTime(indices(end)), ...
            negativeMass(i),negativeP(i),"negative"); %#ok<AGROW>
    end
end

if ~isempty(significantClusters)
    [~,order] = sort([significantClusters.tStart]);
    significantClusters = significantClusters(order);
end

stats = struct();
stats.cellType = cellType;
stats.time = testTime;
stats.tObserved = tObserved;
stats.pObserved = pObserved;
stats.nExcitation = nA;
stats.nControl = nB;
stats.degreesFreedom = degreesFreedom;
stats.tThreshold = tThreshold;
stats.nPermutations = cfg.NPermutations;

fprintf(['Figure 3 %s cluster test: nExc=%d, nCtrl=%d, ' ...
    '%d significant cluster(s).\n'], ...
    cellType,nA,nB,numel(significantClusters));
end


function [tStatistic,pValue] = pooledTwoSampleT(XA,XB)
nA = size(XA,1);
nB = size(XB,1);
degreesFreedom = nA+nB-2;
meanA = mean(XA,1);
meanB = mean(XB,1);
varianceA = var(XA,0,1);
varianceB = var(XB,0,1);
pooledVariance = ((nA-1).*varianceA+(nB-1).*varianceB)./degreesFreedom;
standardError = sqrt(pooledVariance.*(1/nA+1/nB));
tStatistic = (meanA-meanB)./standardError;
tStatistic(~isfinite(tStatistic)) = 0;

if nargout>1
    pValue = betainc(degreesFreedom./ ...
        (degreesFreedom+tStatistic.^2),degreesFreedom/2,0.5);
end
end


function value = studentTInverse(probability,degreesFreedom)
assert(probability>0.5 && probability<1, ...
    'studentTInverse currently expects a probability between 0.5 and 1.');
betaQuantile = betaincinv(2*(1-probability),degreesFreedom/2,0.5);
value = sqrt(degreesFreedom.*(1./betaQuantile-1));
end


function clusters = findClusters1D(mask)
mask = logical(mask(:))';
edges = diff([false mask false]);
starts = find(edges==1);
stops = find(edges==-1)-1;
clusters = cell(numel(starts),1);
for i = 1:numel(starts)
    clusters{i} = starts(i):stops(i);
end
end


function masses = clusterMass(clusters,tValues,direction)
masses = nan(numel(clusters),1);
for i = 1:numel(clusters)
    indices = clusters{i};
    if direction=="positive"
        masses(i) = sum(tValues(indices));
    else
        masses(i) = sum(-tValues(indices));
    end
end
end


function probabilities = correctedClusterP(masses,nullDistribution)
probabilities = nan(numel(masses),1);
for i = 1:numel(masses)
    probabilities(i) = ...
        (1+sum(nullDistribution>=masses(i)))./(numel(nullDistribution)+1);
end
end


function tableOut = emptyClusterTable()
tableOut = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),false(0,1), ...
    'VariableNames',{'cellType','direction','startSec','endSec', ...
    'clusterMass','pCorrected','significant'});
end


function clusters = emptySignificantClusters()
clusters = struct('tStart',{},'tEnd',{},'mass',{},'p',{},'direction',{});
end


function cluster = makeCluster(tStart,tEnd,mass,pValue,direction)
cluster = struct('tStart',tStart,'tEnd',tEnd,'mass',mass, ...
    'p',pValue,'direction',direction);
end


function plotCellPanel(ax,cellType,excitationColor,yLimits, ...
    epochTime,animalCellType,animalGroup,animalTraces, ...
    significantClusters,cfg,showYAxis)

controlColor = [0.58 0.58 0.58];
controlFill = [0.78 0.78 0.78];

useControl = animalCellType==cellType & animalGroup=="control";
useExcitation = animalCellType==cellType & animalGroup=="excitation";
controlTraces = animalTraces(:,useControl);
excitationTraces = animalTraces(:,useExcitation);

controlMean = mean(controlTraces,2,'omitnan');
controlSem = std(controlTraces,0,2,'omitnan')./sqrt(size(controlTraces,2));
excitationMean = mean(excitationTraces,2,'omitnan');
excitationSem = std(excitationTraces,0,2,'omitnan')./ ...
    sqrt(size(excitationTraces,2));

hold(ax,'on');
plot(ax,cfg.PlotWindow,[0 0],'-','Color',[0.82 0.82 0.82], ...
    'LineWidth',0.6);
fillBand(ax,epochTime,controlMean,controlSem,controlFill,0.25);
plot(ax,epochTime,controlMean,'Color',controlColor,'LineWidth',1.15);
fillBand(ax,epochTime,excitationMean,excitationSem,excitationColor,0.18);
plot(ax,epochTime,excitationMean,'Color',excitationColor,'LineWidth',1.35);

xlim(ax,cfg.PlotWindow);
ylim(ax,yLimits);
xticks(ax,-1:1:3);

yRange = diff(yLimits);
barBottom = yLimits(2)-0.075*yRange;
barTop = yLimits(2)-0.035*yRange;
drawPulseTrain(ax,cfg.StimWindow,barBottom,barTop,[1.00 0.62 0.28]);
text(ax,mean(cfg.StimWindow),barTop+0.015*yRange,'TH exc', ...
    'Color',[0.95 0.35 0.10],'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontName','Arial','FontSize',7);

addClusterBars(ax,significantClusters,yLimits);

title(ax,cellType+" Recording",'Color',excitationColor, ...
    'FontName','Arial','FontSize',10,'FontWeight','bold');
xlabel(ax,'Time from stim onset (s)','FontName','Arial','FontSize',8);
if showYAxis
    ylabel(ax,'\DeltaF/F','Interpreter','tex','FontName','Arial','FontSize',9);
end

labelX = cfg.PlotWindow(2)-0.05*diff(cfg.PlotWindow);
text(ax,labelX,yLimits(1)+0.18*yRange, ...
    sprintf('TH Opsin N = %d',size(excitationTraces,2)), ...
    'Color',excitationColor,'HorizontalAlignment','right', ...
    'FontName','Arial','FontSize',7.5,'FontWeight','bold');
text(ax,labelX,yLimits(1)+0.08*yRange, ...
    sprintf('Ctrl N = %d',size(controlTraces,2)), ...
    'Color',controlColor,'HorizontalAlignment','right', ...
    'FontName','Arial','FontSize',7.5);

set(ax,'Box','off','TickDir','out','FontName','Arial','FontSize',8, ...
    'LineWidth',0.75,'Layer','top');
end


function drawPulseTrain(ax,stimWindow,yBottom,yTop,color)
% The stimulation protocol used 20-Hz pulses with 20-ms pulse width.
pulsePeriod = 1/20;
pulseWidth = 0.020;
pulseStarts = stimWindow(1):pulsePeriod:(stimWindow(2)-pulseWidth);
for i = 1:numel(pulseStarts)
    pulseStart = pulseStarts(i);
    pulseEnd = min(pulseStart+pulseWidth,stimWindow(2));
    patch(ax,[pulseStart pulseEnd pulseEnd pulseStart], ...
        [yBottom yBottom yTop yTop],color,'EdgeColor','none');
end
end


function fillBand(ax,time,meanTrace,semTrace,color,alpha)
fill(ax,[time;flipud(time)], ...
    [meanTrace-semTrace;flipud(meanTrace+semTrace)], ...
    color,'FaceAlpha',alpha,'EdgeColor','none');
end


function addClusterBars(ax,clusters,yLimits)
if isempty(clusters)
    return;
end
yRange = diff(yLimits);
for i = 1:numel(clusters)
    if clusters(i).direction=="positive"
        y = yLimits(2)-0.14*yRange;
    else
        y = yLimits(1)+0.04*yRange;
    end
    plot(ax,[clusters(i).tStart clusters(i).tEnd],[y y], ...
        'k-','LineWidth',2,'HandleVisibility','off');
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
    set(fig,'PaperUnits','inches','PaperSize',figureSize, ...
        'PaperPosition',[0 0 figureSize],'PaperPositionMode','manual');
    print(fig,char(pdfFile),'-dpdf','-painters');
    print(fig,char(pngFile),'-dpng','-r600');

    set(fig,'Units',oldUnits,'PaperUnits',oldPaperUnits, ...
        'PaperSize',oldPaperSize,'PaperPosition',oldPaperPosition, ...
        'PaperPositionMode',oldPaperPositionMode);
end

fprintf('Saved Figure 3 PDF: %s\n',pdfFile);
fprintf('Saved Figure 3 PNG: %s\n',pngFile);
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


function tf = isLogicalScalar(x)
tf = (islogical(x) && isscalar(x)) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && ismember(x,[0 1]));
end


function tf = isPositiveInteger(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x>=1 && x==fix(x);
end


function tf = isNonnegativeInteger(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x==fix(x);
end


function tf = isProbability(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<1;
end


function tf = windowWithin(innerWindow,outerWindow)
tf = innerWindow(1)>=outerWindow(1) && innerWindow(2)<=outerWindow(2);
end
