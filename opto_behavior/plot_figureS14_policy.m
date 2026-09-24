function resultsS14 = plot_figureS14_policy(trialFile, metadataFile, posteriorFile, animalFile, outputFolder)
%PLOT_FIGURES14_POLICY Reproduce the behavioral analyses in Figure S14.
%
% resultsS14 = plot_figureS14_policy(trialFile, metadataFile, ...
%     posteriorFile, animalFile, outputFolder)
%
% INPUTS
%   trialFile      modelTfinal_for_rlhmm.csv
%                  Required: animalID, iOrig, jOrig, Y, currChoice
%   metadataFile   rlhmm_plotting_metadata.csv
%                  Required: animalID, iOrig, jOrig, Y, assignment, deltaQ
%   posteriorFile  rlhmm_posteriors.csv
%                  Required: animalID, iOrig, jOrig, Y
%   animalFile     animals.csv
%                  Required: animal_id, cell_type, condition
%   outputFolder   Optional folder for PNG, PDF, and CSV exports. Pass ""
%                  or omit it to suppress file export.
%
% PANELS
%   A-B  P(push) as a function of deltaQ, separately for SST and TH mice.
%   C    Schematic comparison of single- and split-beta policies.
%   D    Per-animal positive- and negative-deltaQ logistic slopes.
%
% The posterior file is used only to retain the exact trial cohort from the
% selected final RL-GLM-HMM fit. State probabilities are not used.
%
% Statistics and Machine Learning Toolbox is required (fitglm and fitlme).

if nargin < 5
    outputFolder = "";
end

trialFile = string(trialFile);
metadataFile = string(metadataFile);
posteriorFile = string(posteriorFile);
animalFile = string(animalFile);
outputFolder = string(outputFolder);

assertFileExists(trialFile);
assertFileExists(metadataFile);
assertFileExists(posteriorFile);
assertFileExists(animalFile);

phaseOrder = ["PRE","OPTO1","OPTO2","POST1","POST2"];
cellTypeOrder = ["sst","th"];
conditionOrder = ["ctrl","ib"];

Ttrial = readtable(trialFile, 'TextType','string');
Tmeta = readtable(metadataFile, 'TextType','string');
Tpost = readtable(posteriorFile, 'TextType','string');
Troster = readtable(animalFile, 'TextType','string');

requireVariables(Ttrial, ["animalID","iOrig","jOrig","Y","currChoice"], "trial file");
requireVariables(Tmeta, ["animalID","iOrig","jOrig","Y","assignment","deltaQ"], "metadata file");
requireVariables(Tpost, ["animalID","iOrig","jOrig","Y"], "posterior file");
requireVariables(Troster, ["animal_id","cell_type","condition"], "animal file");

keyVars = {'animalID','iOrig','jOrig','Y'};
assertUniqueKeys(Ttrial, keyVars, "trial file");
assertUniqueKeys(Tmeta, keyVars, "metadata file");
assertUniqueKeys(Tpost, keyVars, "posterior file");

% Keep only columns used here so joins cannot create ambiguous duplicates.
Ttrial = Ttrial(:, [keyVars, {'currChoice'}]);
Tmeta = Tmeta(:, [keyVars, {'assignment','deltaQ'}]);
Tpost = Tpost(:, keyVars);
Troster = Troster(:, {'animal_id','cell_type','condition'});
Troster.Properties.VariableNames{'animal_id'} = 'animalID';
assertUniqueKeys(Troster, {'animalID'}, "animal file");

T = innerjoin(Ttrial, Tmeta, 'Keys',keyVars);
nTrialMetadata = height(T);
T = innerjoin(T, Tpost, 'Keys',keyVars);
nPosteriorMatched = height(T);

fprintf('Figure S14: %d/%d model-ready trials matched posterior keys.\n', ...
    nPosteriorMatched, nTrialMetadata);

T = innerjoin(T, Troster, 'Keys','animalID');
T.animalID = string(T.animalID);
T.cell_type = lower(strtrim(string(T.cell_type)));
T.condition = lower(strtrim(string(T.condition)));
T.assignment = string(T.assignment);
T.phase = mapPhase(T.assignment);
T.currChoice = toDoubleColumn(T.currChoice);
T.deltaQ = toDoubleColumn(T.deltaQ);

keep = ismember(T.cell_type, cellTypeOrder) & ...
       ismember(T.condition, conditionOrder) & ...
       ismember(T.phase, phaseOrder) & ...
       ismember(T.currChoice, [-1 1]) & ...
       isfinite(T.deltaQ);
T = T(keep,:);

fprintf('Figure S14: %d selected trials, %d animals after phase/cohort filters.\n', ...
    height(T), numel(unique(T.animalID)));

[psychAnimal, psychGroup] = psychometricSummary(T, phaseOrder, ...
    cellTypeOrder, conditionOrder);

[slopeAnimal, slopeGroup, interactionStats, coefficientStats, qcSummary] = ...
    splitBetaAnalysis(T, phaseOrder, cellTypeOrder, conditionOrder);

figPsychometric = plotPsychometric(psychGroup, phaseOrder, ...
    cellTypeOrder, conditionOrder);
figSchematic = plotSchematic();
figSlopes = plotSlopes(slopeAnimal, slopeGroup, interactionStats, ...
    coefficientStats, ...
    phaseOrder, cellTypeOrder, conditionOrder);

if strlength(outputFolder) > 0
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    exportFigurePair(figPsychometric, outputFolder, ...
        "figureS14_panelsA_B_psychometric");
    exportFigurePair(figSchematic, outputFolder, ...
        "figureS14_panelC_schematic");
    exportFigurePair(figSlopes, outputFolder, ...
        "figureS14_panelD_split_beta");

    writetable(psychAnimal, fullfile(outputFolder, ...
        "figureS14_psychometric_animal_bin.csv"));
    writetable(psychGroup, fullfile(outputFolder, ...
        "figureS14_psychometric_group_summary.csv"));
    writetable(slopeAnimal, fullfile(outputFolder, ...
        "figureS14_split_beta_animal_phase.csv"));
    writetable(slopeGroup, fullfile(outputFolder, ...
        "figureS14_split_beta_group_summary.csv"));
    writetable(interactionStats, fullfile(outputFolder, ...
        "figureS14_split_beta_interaction_stats.csv"));
    writetable(coefficientStats, fullfile(outputFolder, ...
        "figureS14_split_beta_coefficient_stats.csv"));
    writetable(qcSummary, fullfile(outputFolder, ...
        "figureS14_split_beta_qc_summary.csv"));
end

resultsS14 = struct();
resultsS14.figures = struct( ...
    'psychometric',figPsychometric, ...
    'schematic',figSchematic, ...
    'splitBeta',figSlopes);
resultsS14.psychometricAnimalBin = psychAnimal;
resultsS14.psychometricGroupSummary = psychGroup;
resultsS14.splitBetaAnimalPhase = slopeAnimal;
resultsS14.splitBetaGroupSummary = slopeGroup;
resultsS14.interactionStats = interactionStats;
resultsS14.coefficientStats = coefficientStats;
resultsS14.qcSummary = qcSummary;
resultsS14.nSelectedTrials = height(T);
resultsS14.nSelectedAnimals = numel(unique(T.animalID));
end


function [animalT, groupT] = psychometricSummary(T, phaseOrder, cellTypes, conditions)
edges = -9:2:9;
centers = edges(1:end-1) + diff(edges)./2;
minTrialsPerAnimalBin = 10;

dqBin = discretize(T.deltaQ, edges);
rows = cell(0,9);

for ct = 1:numel(cellTypes)
    for c = 1:numel(conditions)
        animals = unique(T.animalID(T.cell_type == cellTypes(ct) & ...
            T.condition == conditions(c)));
        for a = 1:numel(animals)
            for p = 1:numel(phaseOrder)
                for b = 1:numel(centers)
                    idx = T.cell_type == cellTypes(ct) & ...
                          T.condition == conditions(c) & ...
                          T.animalID == animals(a) & ...
                          T.phase == phaseOrder(p) & dqBin == b;
                    n = sum(idx);
                    if n == 0
                        continue;
                    end
                    pPush = mean(T.currChoice(idx) == 1);
                    rows(end+1,:) = {cellTypes(ct), conditions(c), ...
                        phaseOrder(p), b, centers(b), animals(a), n, ...
                        pPush, n >= minTrialsPerAnimalBin}; %#ok<AGROW>
                end
            end
        end
    end
end

animalT = cell2table(rows, 'VariableNames', ...
    {'cell_type','condition','phase','deltaQ_bin','deltaQ_center', ...
     'animalID','nTrials','pPush','included'});
animalT.deltaQ_bin = unpackScalarColumn(animalT.deltaQ_bin);
animalT.deltaQ_center = unpackScalarColumn(animalT.deltaQ_center);
animalT.nTrials = unpackScalarColumn(animalT.nTrials);
animalT.pPush = unpackScalarColumn(animalT.pPush);
animalT.included = logical(unpackScalarColumn(animalT.included));
animalT.cell_type = string(animalT.cell_type);
animalT.condition = string(animalT.condition);
animalT.phase = string(animalT.phase);
animalT.animalID = string(animalT.animalID);

groupRows = cell(0,8);
for ct = 1:numel(cellTypes)
    for c = 1:numel(conditions)
        for p = 1:numel(phaseOrder)
            for b = 1:numel(centers)
                idx = animalT.included & ...
                      animalT.cell_type == cellTypes(ct) & ...
                      animalT.condition == conditions(c) & ...
                      animalT.phase == phaseOrder(p) & ...
                      animalT.deltaQ_bin == b;
                values = animalT.pPush(idx);
                if isempty(values)
                    mu = NaN;
                    sem = NaN;
                    nAnimals = 0;
                else
                    mu = mean(values, 'omitnan');
                    nAnimals = sum(isfinite(values));
                    sem = std(values, 'omitnan') ./ sqrt(nAnimals);
                end
                groupRows(end+1,:) = {cellTypes(ct), conditions(c), ...
                    phaseOrder(p), b, centers(b), mu, sem, nAnimals}; %#ok<AGROW>
            end
        end
    end
end

groupT = cell2table(groupRows, 'VariableNames', ...
    {'cell_type','condition','phase','deltaQ_bin','deltaQ_center', ...
     'mean_pPush','sem_pPush','nAnimals'});
groupT.deltaQ_bin = unpackScalarColumn(groupT.deltaQ_bin);
groupT.deltaQ_center = unpackScalarColumn(groupT.deltaQ_center);
groupT.mean_pPush = unpackScalarColumn(groupT.mean_pPush);
groupT.sem_pPush = unpackScalarColumn(groupT.sem_pPush);
groupT.nAnimals = unpackScalarColumn(groupT.nAnimals);
groupT.cell_type = string(groupT.cell_type);
groupT.condition = string(groupT.condition);
groupT.phase = string(groupT.phase);
end


function [slopeT, groupT, interactionT, coefficientT, qcT] = ...
        splitBetaAnalysis(T, phaseOrder, cellTypes, conditions)
minTrialsPerFit = 40;
minPushFrac = 0.05;
maxPushFrac = 0.95;
minSideTrials = 15;
minPerCell = 3;
maxAbsBeta = 6;

slopeRows = cell(0,17);
qcRows = cell(0,8);

for ct = 1:numel(cellTypes)
    Tct = T(T.cell_type == cellTypes(ct),:);
    dqScale = std(Tct.deltaQ, 'omitnan');
    Tct.push = double(Tct.currChoice == 1);
    Tct.deltaQ_pos_z = max(Tct.deltaQ,0) ./ dqScale;
    Tct.deltaQ_neg_z = min(Tct.deltaQ,0) ./ dqScale;

    skippedTrials = 0;
    skippedBalance = 0;
    skippedSides = 0;
    skippedCells = 0;
    skippedExtreme = 0;
    fitFailed = 0;
    accepted = 0;

    animals = unique(Tct.animalID);
    for a = 1:numel(animals)
        for p = 1:numel(phaseOrder)
            idx = Tct.animalID == animals(a) & Tct.phase == phaseOrder(p);
            Ta = Tct(idx,:);

            if height(Ta) < minTrialsPerFit
                skippedTrials = skippedTrials + 1;
                continue;
            end

            pPush = mean(Ta.push);
            if pPush <= minPushFrac || pPush >= maxPushFrac
                skippedBalance = skippedBalance + 1;
                continue;
            end

            nPos = sum(Ta.deltaQ > 0);
            nNeg = sum(Ta.deltaQ < 0);
            if nPos < minSideTrials || nNeg < minSideTrials
                skippedSides = skippedSides + 1;
                continue;
            end

            pushPos = sum(Ta.deltaQ > 0 & Ta.push == 1);
            pullPos = sum(Ta.deltaQ > 0 & Ta.push == 0);
            pushNeg = sum(Ta.deltaQ < 0 & Ta.push == 1);
            pullNeg = sum(Ta.deltaQ < 0 & Ta.push == 0);
            if any([pushPos pullPos pushNeg pullNeg] < minPerCell)
                skippedCells = skippedCells + 1;
                continue;
            end

            try
                mdl = fitglm(Ta, ...
                    'push ~ deltaQ_pos_z + deltaQ_neg_z', ...
                    'Distribution','binomial', 'Link','logit');
                betaPos = coefficientValue(mdl.Coefficients, ...
                    'deltaQ_pos_z', 'Estimate');
                betaNeg = coefficientValue(mdl.Coefficients, ...
                    'deltaQ_neg_z', 'Estimate');
                sePos = coefficientValue(mdl.Coefficients, ...
                    'deltaQ_pos_z', 'SE');
                seNeg = coefficientValue(mdl.Coefficients, ...
                    'deltaQ_neg_z', 'SE');

                if abs(betaPos) > maxAbsBeta || abs(betaNeg) > maxAbsBeta
                    skippedExtreme = skippedExtreme + 1;
                    continue;
                end

                cond = unique(Ta.condition);
                if numel(cond) ~= 1
                    error('FigureS14:ConditionMismatch', ...
                        'Animal %s has multiple conditions.', animals(a));
                end

                accepted = accepted + 1;
                slopeRows(end+1,:) = {cellTypes(ct), animals(a), ...
                    phaseOrder(p), cond(1), height(Ta), pPush, dqScale, ...
                    nPos, nNeg, pushPos, pullPos, pushNeg, pullNeg, ...
                    betaPos, betaNeg, sePos, seNeg}; %#ok<AGROW>
            catch ME
                fitFailed = fitFailed + 1;
                warning('FigureS14:FitFailed', ...
                    'Split-beta fit failed for %s, %s: %s', ...
                    animals(a), phaseOrder(p), ME.message);
            end
        end
    end

    qcRows(end+1,:) = {cellTypes(ct), accepted, skippedTrials, ...
        skippedBalance, skippedSides, skippedCells, skippedExtreme, ...
        fitFailed}; %#ok<AGROW>
end

slopeT = cell2table(slopeRows, 'VariableNames', ...
    {'cell_type','animalID','phase','condition','nTrials','pPush', ...
     'deltaQ_scale','nDeltaQ_pos','nDeltaQ_neg','nPush_posDQ', ...
     'nPull_posDQ','nPush_negDQ','nPull_negDQ','beta_pos','beta_neg', ...
     'se_pos','se_neg'});
slopeT.cell_type = string(slopeT.cell_type);
slopeT.animalID = string(slopeT.animalID);
slopeT.phase = string(slopeT.phase);
slopeT.condition = string(slopeT.condition);
numericVars = {'nTrials','pPush','deltaQ_scale','nDeltaQ_pos', ...
    'nDeltaQ_neg','nPush_posDQ','nPull_posDQ','nPush_negDQ', ...
    'nPull_negDQ','beta_pos','beta_neg','se_pos','se_neg'};
for i = 1:numel(numericVars)
    slopeT.(numericVars{i}) = unpackScalarColumn(slopeT.(numericVars{i}));
end

qcT = cell2table(qcRows, 'VariableNames', ...
    {'cell_type','acceptedFits','skippedTrialCount','skippedChoiceBalance', ...
     'skippedDeltaQSides','skippedChoiceCells','skippedExtremeBeta', ...
     'failedFits'});
qcT.cell_type = string(qcT.cell_type);
for i = 2:width(qcT)
    variable = qcT.Properties.VariableNames{i};
    qcT.(variable) = unpackScalarColumn(qcT.(variable));
end

groupRows = cell(0,8);
metrics = ["beta_pos","beta_neg"];
for ct = 1:numel(cellTypes)
    for m = 1:numel(metrics)
        for c = 1:numel(conditions)
            for p = 1:numel(phaseOrder)
                idx = slopeT.cell_type == cellTypes(ct) & ...
                      slopeT.condition == conditions(c) & ...
                      slopeT.phase == phaseOrder(p);
                values = slopeT.(metrics(m))(idx);
                n = sum(isfinite(values));
                if n == 0
                    mu = NaN;
                    sem = NaN;
                else
                    mu = mean(values, 'omitnan');
                    sem = std(values, 'omitnan') ./ sqrt(n);
                end
                groupRows(end+1,:) = {cellTypes(ct), metrics(m), ...
                    conditions(c), phaseOrder(p), mu, sem, n, ...
                    sum(slopeT.nTrials(idx))}; %#ok<AGROW>
            end
        end
    end
end

groupT = cell2table(groupRows, 'VariableNames', ...
    {'cell_type','metric','condition','phase','meanSlope','semSlope', ...
     'nAnimals','nTrials'});
groupT.cell_type = string(groupT.cell_type);
groupT.metric = string(groupT.metric);
groupT.condition = string(groupT.condition);
groupT.phase = string(groupT.phase);
for i = 5:width(groupT)
    variable = groupT.Properties.VariableNames{i};
    groupT.(variable) = unpackScalarColumn(groupT.(variable));
end

interactionRows = cell(0,5);
coefficientRows = cell(0,11);
for ct = 1:numel(cellTypes)
    for m = 1:numel(metrics)
        idx = slopeT.cell_type == cellTypes(ct);
        Tm = slopeT(idx, {'animalID','phase','condition'});
        Tm.response = slopeT.(metrics(m))(idx);
        Tm.animalID = categorical(Tm.animalID);
        Tm.phase = categorical(Tm.phase, phaseOrder, 'Ordinal',true);
        Tm.condition = categorical(Tm.condition, conditions);

        mdl = fitlme(Tm, ...
            'response ~ condition * phase + (1|animalID)');
        a = anova(mdl);
        pInteraction = termPValue(a, "condition:phase");
        interactionRows(end+1,:) = {cellTypes(ct), metrics(m), ...
            pInteraction, height(Tm), numel(categories(Tm.animalID))}; %#ok<AGROW>

        c = mdl.Coefficients;
        names = coefficientNames(mdl, c);
        for k = 1:numel(names)
            coefficientRows(end+1,:) = {cellTypes(ct), metrics(m), ...
                names(k), statValue(c,k,'Estimate'), ...
                statValue(c,k,'SE'), statValue(c,k,'tStat'), ...
                statValue(c,k,'DF'), statValue(c,k,'pValue'), ...
                statValue(c,k,'Lower'), statValue(c,k,'Upper'), ...
                coefficientPhase(names(k), phaseOrder)}; %#ok<AGROW>
        end
    end
end

interactionT = cell2table(interactionRows, 'VariableNames', ...
    {'cell_type','metric','pInteraction','nAnimalPhaseFits','nAnimals'});
interactionT.cell_type = string(interactionT.cell_type);
interactionT.metric = string(interactionT.metric);
interactionT.pInteraction = unpackScalarColumn(interactionT.pInteraction);
interactionT.nAnimalPhaseFits = unpackScalarColumn(interactionT.nAnimalPhaseFits);
interactionT.nAnimals = unpackScalarColumn(interactionT.nAnimals);

coefficientT = cell2table(coefficientRows, 'VariableNames', ...
    {'cell_type','metric','coefficient','estimate','SE','tStat','DF', ...
     'pValue','lowerCI','upperCI','phase'});
coefficientT.cell_type = string(coefficientT.cell_type);
coefficientT.metric = string(coefficientT.metric);
coefficientT.coefficient = string(coefficientT.coefficient);
coefficientT.phase = string(coefficientT.phase);
for i = 4:10
    variable = coefficientT.Properties.VariableNames{i};
    coefficientT.(variable) = unpackScalarColumn(coefficientT.(variable));
end
end


function fig = plotPsychometric(groupT, phaseOrder, cellTypes, conditions)
colors = struct('ctrl',[0 0 0], 'ib',[0.93 0.10 0.12]);
fig = figure('Color','w','Position',[100 100 1300 520]);
tiledlayout(2,5, 'TileSpacing','compact', 'Padding','compact');

for ct = 1:numel(cellTypes)
    for p = 1:numel(phaseOrder)
        ax = nexttile;
        hold(ax,'on');
        for c = 1:numel(conditions)
            idx = groupT.cell_type == cellTypes(ct) & ...
                  groupT.condition == conditions(c) & ...
                  groupT.phase == phaseOrder(p);
            S = sortrows(groupT(idx,:), 'deltaQ_center');
            errorbar(ax, S.deltaQ_center, S.mean_pPush, S.sem_pPush, ...
                '-o', 'Color',colors.(conditions(c)), ...
                'MarkerFaceColor',colors.(conditions(c)), ...
                'MarkerEdgeColor','none', 'MarkerSize',4, ...
                'LineWidth',0.8, 'CapSize',3);
        end
        xlim(ax,[-8 8]);
        ylim(ax,[0 1]);
        xticks(ax,[-8 -4 0 4 8]);
        yticks(ax,[0 0.5 1]);
        title(ax, phaseDisplayName(phaseOrder(p)));
        xlabel(ax,'\DeltaQ');
        if p == 1
            ylabel(ax, {char(upper(cellTypes(ct))); 'P(push)'}, ...
                'Interpreter','none');
        end
        box(ax,'off');
        set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.8);
    end
end
end


function fig = plotSchematic()
x = linspace(-8,8,500);
beta = 0.75;
pSingle = 1 ./ (1 + exp(-beta.*x));

betaPos = 0.20;
betaNeg = 1.20;
xPos = max(x,0);
xNeg = min(x,0);
pSplit = 1 ./ (1 + exp(-(betaPos.*xPos + betaNeg.*xNeg)));

green = [0.35 0.75 0.20];
cyan = [0.00 0.68 0.72];
magenta = [0.65 0.20 0.65];

fig = figure('Color','w','Position',[100 100 320 520]);
tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

ax = nexttile;
plot(ax,x,pSingle,'Color',[0.60 0.60 0.60],'LineWidth',1.0);
hold(ax,'on');
plot(ax,x,pSingle,'--','Color',green,'LineWidth',1.6);
text(ax,1.6,0.28,'\beta','Color',green,'FontWeight','bold');
formatSchematicAxis(ax,'Single \beta');

ax = nexttile;
plot(ax,x,pSplit,'Color',[0.45 0.45 0.45],'LineWidth',1.0);
hold(ax,'on');
neg = x <= 0;
pos = x >= 0;
plot(ax,x(neg),pSplit(neg),'--','Color',magenta,'LineWidth',1.6);
plot(ax,x(pos),pSplit(pos),'--','Color',cyan,'LineWidth',1.6);
text(ax,-5.4,0.25,'\beta^-','Color',magenta,'FontWeight','bold');
text(ax,3.2,0.68,'\beta^+','Color',cyan,'FontWeight','bold');
formatSchematicAxis(ax,'Split \beta');
end


function fig = plotSlopes(animalT, groupT, interactionT, coefficientT, ...
        phaseOrder, cellTypes, conditions)
groupColors = struct('ctrl',[0 0 0], 'ib',[0.93 0.10 0.12]);
lineColors = struct('ctrl',[0.72 0.72 0.72], 'ib',[1.00 0.72 0.72]);
metrics = ["beta_pos","beta_neg"];

fig = figure('Color','w','Position',[100 100 800 650]);
tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

for ct = 1:numel(cellTypes)
    for m = 1:numel(metrics)
        ax = nexttile;
        hold(ax,'on');
        patch(ax,[1.5 3.5 3.5 1.5],[0 0 5 5], ...
            [1 0.85 0.85], 'EdgeColor','none', ...
            'FaceAlpha',0.45, 'HandleVisibility','off');

        for c = 1:numel(conditions)
            animals = unique(animalT.animalID( ...
                animalT.cell_type == cellTypes(ct) & ...
                animalT.condition == conditions(c)));
            for a = 1:numel(animals)
                idx = animalT.cell_type == cellTypes(ct) & ...
                      animalT.condition == conditions(c) & ...
                      animalT.animalID == animals(a);
                A = animalT(idx,:);
                [tf,loc] = ismember(A.phase, phaseOrder);
                x = loc(tf);
                y = A.(metrics(m))(tf);
                [x,ord] = sort(x);
                y = y(ord);
                if numel(x) >= 2
                    plot(ax,x,y,'-','Color',lineColors.(conditions(c)), ...
                        'LineWidth',0.55,'HandleVisibility','off');
                end
            end
        end

        for c = 1:numel(conditions)
            idx = groupT.cell_type == cellTypes(ct) & ...
                  groupT.metric == metrics(m) & ...
                  groupT.condition == conditions(c);
            S = groupT(idx,:);
            [tf,loc] = ismember(S.phase, phaseOrder);
            x = loc(tf);
            [x,ord] = sort(x);
            S = S(tf,:);
            S = S(ord,:);
            errorbar(ax,x,S.meanSlope,S.semSlope,'-o', ...
                'Color',groupColors.(conditions(c)), ...
                'MarkerFaceColor',groupColors.(conditions(c)), ...
                'MarkerEdgeColor','none','MarkerSize',4, ...
                'LineWidth',0.9,'CapSize',3);
        end

        p = phaseContrastPValues(coefficientT, cellTypes(ct), ...
            metrics(m), phaseOrder);
        for k = 1:numel(phaseOrder)
            if isfinite(p(k)) && p(k) < 0.05
                text(ax,k,4.60,'*','HorizontalAlignment','center', ...
                    'FontWeight','bold','FontSize',13);
            end
        end

        xlim(ax,[0.65 5.35]);
        ylim(ax,[0 5]);
        xticks(ax,1:5);
        xticklabels(ax,["Pre","Opto I","Opto II","Post I","Post II"]);
        xtickangle(ax,45);
        ylabel(ax,'Logistic slope');
        if metrics(m) == "beta_pos"
            title(ax,'\beta^+ (Q_{push} > Q_{pull})', ...
                'Color',[0 0.60 0.65]);
        else
            title(ax,'\beta^- (Q_{pull} > Q_{push})', ...
                'Color',[0.60 0.15 0.60]);
        end
        statIdx = interactionT.cell_type == cellTypes(ct) & ...
                  interactionT.metric == metrics(m);
        if sum(statIdx) == 1
            pInteraction = interactionT.pInteraction(statIdx);
            text(ax,0.52,0.92, ...
                sprintf('p_{interaction} = %.3g',pInteraction), ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontSize',8);
        end
        if m == 1
            text(ax,-0.22,0.5,upper(cellTypes(ct)), ...
                'Units','normalized','Rotation',90, ...
                'HorizontalAlignment','center','FontWeight','bold');
        end
        box(ax,'off');
        set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.8);
    end
end
end


function formatSchematicAxis(ax, ttl)
xlim(ax,[-8 8]);
ylim(ax,[0 1]);
xticks(ax,[-8 0 8]);
yticks(ax,[0 1]);
xlabel(ax,'\DeltaQ');
ylabel(ax,'P(push)');
title(ax,ttl);
text(ax,-7.8,-0.16,'pull better','Color',[0.65 0.20 0.65], ...
    'FontAngle','italic','FontSize',8);
text(ax,2.0,-0.16,'push better','Color',[0.00 0.68 0.72], ...
    'FontAngle','italic','FontSize',8);
box(ax,'off');
set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.8);
end


function p = phaseContrastPValues(coefT, cellType, metric, phaseOrder)
p = nan(size(phaseOrder));
S = coefT(coefT.cell_type == cellType & coefT.metric == metric,:);
for k = 1:numel(phaseOrder)
    idx = S.phase == phaseOrder(k);
    if sum(idx) == 1
        p(k) = S.pValue(idx);
    end
end
end


function phase = coefficientPhase(name, phaseOrder)
name = string(name);
phase = "";
hasCondition = contains(lower(name),'condition_ib');
if ~hasCondition
    return;
end
if strcmpi(name,'condition_ib')
    phase = "PRE";
    return;
end
for k = 2:numel(phaseOrder)
    if contains(upper(name), "PHASE_" + phaseOrder(k))
        phase = phaseOrder(k);
        return;
    end
end
end


function value = coefficientValue(coefT, rowName, variable)
vars = statVariableNames(coefT);
varIdx = find(strcmpi(vars,variable),1);
if isempty(varIdx)
    error('FigureS14:MissingCoefficientField', ...
        'Coefficient table is missing %s.', variable);
end
rowNames = string(coefT.Properties.RowNames);
rowIdx = find(strcmp(rowNames,rowName),1);
if isempty(rowIdx)
    error('FigureS14:MissingCoefficient', ...
        'Coefficient %s was not returned by fitglm.', rowName);
end
values = coefT.(vars{varIdx});
value = values(rowIdx);
end


function names = coefficientNames(mdl, coefT)
try
    names = string(mdl.CoefficientNames(:));
catch
    names = string(coefT.Properties.RowNames);
end
if isempty(names)
    names = "coefficient_" + string((1:height(coefT))');
end
end


function value = statValue(stats, row, variable)
vars = statVariableNames(stats);
idx = find(strcmpi(vars,variable),1);
if isempty(idx)
    value = NaN;
else
    values = stats.(vars{idx});
    value = values(row);
end
end


function p = termPValue(stats, targetTerm)
vars = statVariableNames(stats);
termIdx = find(strcmpi(vars,'Term'),1);
pIdx = find(strcmpi(vars,'pValue'),1);
if isempty(termIdx) || isempty(pIdx)
    error('FigureS14:UnexpectedAnovaOutput', ...
        'ANOVA output does not contain Term and pValue.');
end

terms = lower(erase(string(stats.(vars{termIdx})), " "));
target = lower(erase(string(targetTerm), " "));
row = find(terms == target,1);
if isempty(row)
    targetParts = split(target,":");
    row = find(contains(terms,targetParts(1)) & ...
               contains(terms,targetParts(2)) & contains(terms,":"),1);
end
if isempty(row)
    error('FigureS14:MissingInteraction', ...
        'Could not find %s in the ANOVA output.', targetTerm);
end
values = stats.(vars{pIdx});
p = values(row);
end


function vars = statVariableNames(stats)
try
    vars = stats.Properties.VariableNames;
    return;
catch
end
try
    vars = stats.Properties.VarNames;
    return;
catch
end
error('FigureS14:UnknownStatsContainer', ...
    'Could not determine variable names in statistical output.');
end


function phase = mapPhase(assignment)
a = upper(strtrim(string(assignment)));
phase = strings(size(a));
phase(ismember(a,["PRE","PRE1","PRE2"])) = "PRE";
phase(ismember(a,["O1","O2","O3","OPTO1"])) = "OPTO1";
phase(ismember(a,["O4","O5","O6","BINTRA-OPTO"])) = "OPTO2";
phase(ismember(a,["POST1","POST2"])) = "POST1";
phase(ismember(a,["POST3","POST4"])) = "POST2";
end


function label = phaseDisplayName(phase)
switch string(phase)
    case "PRE"
        label = "Pre";
    case "OPTO1"
        label = "Opto I";
    case "OPTO2"
        label = "Opto II";
    case "POST1"
        label = "Post I";
    case "POST2"
        label = "Post II";
end
end


function exportFigurePair(fig, outputFolder, stem)
pngFile = fullfile(outputFolder, stem + ".png");
pdfFile = fullfile(outputFolder, stem + ".pdf");
if exist('exportgraphics','file') == 2
    exportgraphics(fig, pngFile, 'Resolution',300);
    exportgraphics(fig, pdfFile, 'ContentType','vector');
else
    set(fig, 'PaperPositionMode','auto');
    print(fig, pngFile, '-dpng', '-r300');
    print(fig, pdfFile, '-dpdf', '-painters', '-bestfit');
end
end


function assertFileExists(file)
if ~isfile(file)
    error('FigureS14:FileNotFound','File not found: %s',file);
end
end


function requireVariables(T, required, label)
missing = setdiff(required, string(T.Properties.VariableNames));
if ~isempty(missing)
    error('FigureS14:MissingVariables', ...
        '%s is missing required variables: %s', ...
        label, strjoin(missing,', '));
end
end


function assertUniqueKeys(T, keyVars, label)
keyT = T(:,keyVars);
if height(unique(keyT,'rows')) ~= height(keyT)
    error('FigureS14:DuplicateKeys', ...
        '%s contains duplicate trial keys.',label);
end
end


function values = unpackScalarColumn(values)
if iscell(values)
    values = vertcat(values{:});
end
end


function values = toDoubleColumn(values)
if isnumeric(values) || islogical(values)
    values = double(values);
else
    values = str2double(string(values));
end
end
