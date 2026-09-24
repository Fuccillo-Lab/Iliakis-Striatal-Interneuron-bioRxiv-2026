function resultsS13 = plot_figureS13_wsls(trialFile, metadataFile, posteriorFile, animalFile, outputFolder)
%PLOT_FIGURES13_WSLS Reproduce the behavioral analyses in Figure S13.
%
% resultsS13 = plot_figureS13_wsls(trialFile, metadataFile, ...
%     posteriorFile, animalFile, outputFolder)
%
% INPUTS
%   trialFile      modelTfinal_for_rlhmm.csv
%                  Required: animalID, iOrig, jOrig, Y, currChoice,
%                  currReward, prevChoice, prevReward, isOptoSession,
%                  previousOptoTrial
%   metadataFile   rlhmm_plotting_metadata.csv
%                  Required: animalID, iOrig, jOrig, Y, assignment
%   posteriorFile  rlhmm_posteriors.csv
%                  Required: animalID, iOrig, jOrig, Y
%   animalFile     animals.csv
%                  Required: animal_id, cell_type, condition
%   outputFolder   Optional output folder. Pass "" or omit to suppress
%                  file export.
%
% PANELS
%   A-H  PRE group comparisons and within-session Light OFF/ON behavior.
%        The displayed opto-session p value tests the condition main effect
%        in the condition-by-light GLME. Light and interaction terms, plus
%        direct OFF-versus-ON contrasts within each group, are exported.
%   I-L  Lose-switch behavior across PRE, OPTO I, OPTO II, POST I, POST II.
%
% The posterior file is used only to retain the exact historical cohort.
% State probabilities are not used in Figure S13.
%
% Statistics and Machine Learning Toolbox is required (fitglme).

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
actionOrder = ["Push","Pull"];

Ttrial = readtable(trialFile, 'TextType','string');
Tmeta = readtable(metadataFile, 'TextType','string');
Tpost = readtable(posteriorFile, 'TextType','string');
Troster = readtable(animalFile, 'TextType','string');

trialRequired = ["animalID","iOrig","jOrig","Y","currChoice", ...
    "currReward","prevChoice","prevReward","isOptoSession", ...
    "previousOptoTrial"];
requireVariables(Ttrial, trialRequired, "trial file");
requireVariables(Tmeta, ["animalID","iOrig","jOrig","Y","assignment"], ...
    "metadata file");
requireVariables(Tpost, ["animalID","iOrig","jOrig","Y"], ...
    "posterior file");
requireVariables(Troster, ["animal_id","cell_type","condition"], ...
    "animal file");

keyVars = {'animalID','iOrig','jOrig','Y'};
assertUniqueKeys(Ttrial, keyVars, "trial file");
assertUniqueKeys(Tmeta, keyVars, "metadata file");
assertUniqueKeys(Tpost, keyVars, "posterior file");

Ttrial = Ttrial(:, [keyVars, {'currChoice','currReward','prevChoice', ...
    'prevReward','isOptoSession','previousOptoTrial'}]);
Tmeta = Tmeta(:, [keyVars, {'assignment'}]);
Tpost = Tpost(:, keyVars);
Troster = Troster(:, {'animal_id','cell_type','condition'});
Troster.Properties.VariableNames{'animal_id'} = 'animalID';
assertUniqueKeys(Troster, {'animalID'}, "animal file");

T = innerjoin(Ttrial, Tmeta, 'Keys',keyVars);
nTrialMetadata = height(T);
T = innerjoin(T, Tpost, 'Keys',keyVars);
nPosteriorMatched = height(T);

fprintf('Figure S13: %d/%d model-ready trials matched posterior keys.\n', ...
    nPosteriorMatched, nTrialMetadata);

T = innerjoin(T, Troster, 'Keys','animalID');
T.animalID = string(T.animalID);
T.cell_type = lower(strtrim(string(T.cell_type)));
T.condition = lower(strtrim(string(T.condition)));
T.phase = mapPhase(T.assignment);

numericVars = ["currChoice","currReward","prevChoice","prevReward", ...
    "isOptoSession","previousOptoTrial","iOrig","jOrig"];
for k = 1:numel(numericVars)
    T.(numericVars(k)) = toDoubleColumn(T.(numericVars(k)));
end

valid = ismember(T.cell_type,cellTypeOrder) & ...
        ismember(T.condition,conditionOrder) & ...
        ismember(T.phase,phaseOrder) & ...
        ismember(T.currChoice,[-1 1]) & ...
        ismember(T.prevChoice,[-1 1]) & ...
        ismember(T.prevReward,[0 1]);
T = T(valid,:);

T.stay = double(T.currChoice == T.prevChoice);
T.switch = double(T.currChoice ~= T.prevChoice);
T.prevAction = strings(height(T),1);
T.prevAction(T.prevChoice == 1) = "Push";
T.prevAction(T.prevChoice == -1) = "Pull";
T.sessionID = T.animalID + "_" + string(T.iOrig);

fprintf('Figure S13: %d selected trials, %d animals after phase/cohort filters.\n', ...
    height(T), numel(unique(T.animalID)));

[acuteAnimal, acuteModelStats, acuteCoefficientStats, acuteLightEffects] = ...
    acuteAnalysis(T, cellTypeOrder, conditionOrder, actionOrder);

[sustainedAnimal, sustainedGroup, sustainedModelStats, ...
    sustainedCoefficientStats] = sustainedAnalysis(T, phaseOrder, ...
    cellTypeOrder, conditionOrder, actionOrder);

fig = plotFigureS13(acuteAnimal, acuteModelStats, sustainedAnimal, ...
    sustainedGroup, sustainedModelStats, sustainedCoefficientStats, ...
    phaseOrder, cellTypeOrder, conditionOrder, actionOrder);

if strlength(outputFolder) > 0
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    exportFigurePair(fig, outputFolder, "figureS13_wsls");
    writetable(acuteAnimal, fullfile(outputFolder, ...
        "figureS13_acute_animal_summary.csv"));
    writetable(acuteModelStats, fullfile(outputFolder, ...
        "figureS13_acute_model_stats.csv"));
    writetable(acuteCoefficientStats, fullfile(outputFolder, ...
        "figureS13_acute_coefficient_stats.csv"));
    writetable(acuteLightEffects, fullfile(outputFolder, ...
        "figureS13_acute_light_simple_effects.csv"));
    writetable(sustainedAnimal, fullfile(outputFolder, ...
        "figureS13_sustained_animal_phase.csv"));
    writetable(sustainedGroup, fullfile(outputFolder, ...
        "figureS13_sustained_group_summary.csv"));
    writetable(sustainedModelStats, fullfile(outputFolder, ...
        "figureS13_sustained_model_stats.csv"));
    writetable(sustainedCoefficientStats, fullfile(outputFolder, ...
        "figureS13_sustained_coefficient_stats.csv"));
end

resultsS13 = struct();
resultsS13.figure = fig;
resultsS13.acuteAnimalSummary = acuteAnimal;
resultsS13.acuteModelStats = acuteModelStats;
resultsS13.acuteCoefficientStats = acuteCoefficientStats;
resultsS13.acuteLightSimpleEffects = acuteLightEffects;
resultsS13.sustainedAnimalPhase = sustainedAnimal;
resultsS13.sustainedGroupSummary = sustainedGroup;
resultsS13.sustainedModelStats = sustainedModelStats;
resultsS13.sustainedCoefficientStats = sustainedCoefficientStats;
resultsS13.nSelectedTrials = height(T);
resultsS13.nSelectedAnimals = numel(unique(T.animalID));
end


function [animalT, modelStatsT, coefficientT, lightEffectsT] = ...
        acuteAnalysis(T, cellTypes, conditions, actions)
metricNames = ["winStay","loseSwitch"];
animalRows = cell(0,9);
modelStatRows = cell(0,12);
coefficientRows = cell(0,13);
lightRows = cell(0,10);

for ct = 1:numel(cellTypes)
    Tct = T(T.cell_type == cellTypes(ct),:);

    for m = 1:numel(metricNames)
        metric = metricNames(m);
        if metric == "winStay"
            outcomeValue = 1;
            responseName = "stay";
        else
            outcomeValue = 0;
            responseName = "switch";
        end

        for a = 1:numel(actions)
            action = actions(a);
            baseMask = Tct.phase == "PRE" & ...
                Tct.prevReward == outcomeValue & ...
                Tct.prevAction == action;
            optoMask = ismember(Tct.phase,["OPTO1","OPTO2"]) & ...
                Tct.isOptoSession == 1 & ...
                ismember(Tct.previousOptoTrial,[0 1]) & ...
                Tct.prevReward == outcomeValue & ...
                Tct.prevAction == action;

            Tpre = Tct(baseMask,:);
            Topto = Tct(optoMask,:);
            Tpre.response = Tpre.(responseName);
            Topto.response = Topto.(responseName);

            % Animal-level values used for plotting.
            for c = 1:numel(conditions)
                cond = conditions(c);
                animals = unique(Tpre.animalID(Tpre.condition == cond));
                for i = 1:numel(animals)
                    idx = Tpre.condition == cond & ...
                          Tpre.animalID == animals(i);
                    animalRows(end+1,:) = {cellTypes(ct),metric,action, ...
                        cond,"PRE","PRE",animals(i),sum(idx), ...
                        mean(Tpre.response(idx),'omitnan')}; %#ok<AGROW>
                end

                animals = unique(Topto.animalID(Topto.condition == cond));
                for i = 1:numel(animals)
                    for light = 0:1
                        idx = Topto.condition == cond & ...
                              Topto.animalID == animals(i) & ...
                              Topto.previousOptoTrial == light;
                        if sum(idx) == 0
                            continue;
                        end
                        lightName = "Light Off";
                        if light == 1
                            lightName = "Light On";
                        end
                        animalRows(end+1,:) = {cellTypes(ct),metric,action, ...
                            cond,"OPTO",lightName,animals(i),sum(idx), ...
                            mean(Topto.response(idx),'omitnan')}; %#ok<AGROW>
                    end
                end
            end

            % PRE condition comparison.
            Tpre.condition = categorical(Tpre.condition,conditions);
            Tpre.animalID = categorical(Tpre.animalID);
            Tpre.sessionID = categorical(Tpre.sessionID);
            mdlPre = fitglme(Tpre, ...
                'response ~ condition + (1|animalID) + (1|sessionID)', ...
                'Distribution','Binomial','Link','logit');
            modelStatRows = appendAnovaRows(modelStatRows,anova(mdlPre), ...
                cellTypes(ct),"PRE",metric,action,height(Tpre), ...
                numel(categories(Tpre.animalID)));
            coefficientRows = appendCoefficientRows(coefficientRows,mdlPre, ...
                cellTypes(ct),"PRE",metric,action,height(Tpre));

            % OPTO condition-by-light model.
            Topto.condition = categorical(Topto.condition,conditions);
            Topto.prevLight = categorical(Topto.previousOptoTrial, ...
                [0 1],["Light Off","Light On"]);
            Topto.animalID = categorical(Topto.animalID);
            Topto.sessionID = categorical(Topto.sessionID);
            mdlOpto = fitglme(Topto, ...
                ['response ~ condition * prevLight + ' ...
                 '(1+prevLight|animalID) + (1|sessionID)'], ...
                'Distribution','Binomial','Link','logit');
            modelStatRows = appendAnovaRows(modelStatRows,anova(mdlOpto), ...
                cellTypes(ct),"OPTO",metric,action,height(Topto), ...
                numel(categories(Topto.animalID)));
            coefficientRows = appendCoefficientRows(coefficientRows,mdlOpto, ...
                cellTypes(ct),"OPTO",metric,action,height(Topto));
            lightRows = [lightRows; lightSimpleEffectRows(mdlOpto, ...
                cellTypes(ct),metric,action)]; %#ok<AGROW>
        end
    end
end

animalT = cell2table(animalRows,'VariableNames', ...
    {'cell_type','metric','prevAction','condition','period','light', ...
     'animalID','nTrials','value'});
animalT = normalizeMixedTable(animalT,1:7,8:9);

modelStatsT = cell2table(modelStatRows,'VariableNames', ...
    {'cell_type','period','metric','prevAction','term','statisticName', ...
     'statistic','DF1','DF2','pValue','nTrials','nAnimals'});
modelStatsT = normalizeMixedTable(modelStatsT,1:6,7:12);

coefficientT = cell2table(coefficientRows,'VariableNames', ...
    {'cell_type','period','metric','prevAction','coefficient','estimate', ...
     'SE','statistic','DF','pValue','lowerCI','upperCI','nTrials'});
coefficientT = normalizeMixedTable(coefficientT,1:5,6:13);

lightEffectsT = cell2table(lightRows,'VariableNames', ...
    {'cell_type','metric','prevAction','condition', ...
     'lightOnMinusOff_logOdds','oddsRatio','FStat','DF1','DF2','pValue'});
lightEffectsT = normalizeMixedTable(lightEffectsT,1:4,5:10);
end


function [animalT, groupT, modelStatsT, coefficientT] = ...
        sustainedAnalysis(T, phaseOrder, cellTypes, conditions, actions)
animalRows = cell(0,8);

for ct = 1:numel(cellTypes)
    for a = 1:numel(actions)
        for c = 1:numel(conditions)
            animals = unique(T.animalID(T.cell_type == cellTypes(ct) & ...
                T.condition == conditions(c)));
            for i = 1:numel(animals)
                for p = 1:numel(phaseOrder)
                    idx = T.cell_type == cellTypes(ct) & ...
                          T.condition == conditions(c) & ...
                          T.animalID == animals(i) & ...
                          T.phase == phaseOrder(p) & ...
                          T.prevAction == actions(a) & ...
                          T.prevReward == 0;
                    if sum(idx) == 0
                        continue;
                    end
                    animalRows(end+1,:) = {cellTypes(ct),actions(a), ...
                        conditions(c),phaseOrder(p),animals(i),sum(idx), ...
                        mean(T.switch(idx),'omitnan'),"all model-ready"}; %#ok<AGROW>
                end
            end
        end
    end
end

animalT = cell2table(animalRows,'VariableNames', ...
    {'cell_type','prevAction','condition','phase','animalID','nTrials', ...
     'pLoseSwitch','plotCohort'});
animalT = normalizeMixedTable(animalT,[1:5 8],6:7);

groupRows = cell(0,8);
for ct = 1:numel(cellTypes)
    for a = 1:numel(actions)
        for c = 1:numel(conditions)
            for p = 1:numel(phaseOrder)
                idx = animalT.cell_type == cellTypes(ct) & ...
                      animalT.prevAction == actions(a) & ...
                      animalT.condition == conditions(c) & ...
                      animalT.phase == phaseOrder(p);
                values = animalT.pLoseSwitch(idx);
                n = sum(isfinite(values));
                mu = NaN;
                sem = NaN;
                if n > 0
                    mu = mean(values,'omitnan');
                    sem = std(values,'omitnan') ./ sqrt(n);
                end
                groupRows(end+1,:) = {cellTypes(ct),actions(a), ...
                    conditions(c),phaseOrder(p),mu,sem,n, ...
                    sum(animalT.nTrials(idx))}; %#ok<AGROW>
            end
        end
    end
end

groupT = cell2table(groupRows,'VariableNames', ...
    {'cell_type','prevAction','condition','phase','mean_pLoseSwitch', ...
     'sem_pLoseSwitch','nAnimals','nTrials'});
groupT = normalizeMixedTable(groupT,1:4,5:8);

% The source statistics require adjacent model-ready rows within a session.
Tstats = sortrows(T,{'animalID','iOrig','jOrig'});
isConsecutive = false(height(Tstats),1);
if height(Tstats) > 1
    isConsecutive(2:end) = ...
        Tstats.animalID(2:end) == Tstats.animalID(1:end-1) & ...
        Tstats.iOrig(2:end) == Tstats.iOrig(1:end-1) & ...
        Tstats.jOrig(2:end) == Tstats.jOrig(1:end-1) + 1;
end
Tstats = Tstats(isConsecutive & Tstats.prevReward == 0,:);

modelStatRows = cell(0,12);
coefficientRows = cell(0,14);
for ct = 1:numel(cellTypes)
    for a = 1:numel(actions)
        idx = Tstats.cell_type == cellTypes(ct) & ...
              Tstats.prevAction == actions(a);
        Tm = Tstats(idx,:);
        Tm.condition = categorical(Tm.condition,conditions);
        Tm.phase = categorical(Tm.phase,phaseOrder,'Ordinal',true);
        Tm.animalID = categorical(Tm.animalID);
        Tm.sessionID = categorical(Tm.sessionID);
        Tm.response = Tm.switch;

        mdl = fitglme(Tm, ...
            ['response ~ condition * phase + ' ...
             '(1|animalID) + (1|sessionID)'], ...
            'Distribution','Binomial','Link','logit');
        modelStatRows = appendAnovaRows(modelStatRows,anova(mdl), ...
            cellTypes(ct),"ALL", "loseSwitch",actions(a),height(Tm), ...
            numel(categories(Tm.animalID)));

        baseRows = appendCoefficientRows(cell(0,13),mdl,cellTypes(ct), ...
            "ALL","loseSwitch",actions(a),height(Tm));
        for r = 1:size(baseRows,1)
            coefficientRows(end+1,:) = [baseRows(r,:), ...
                {coefficientPhase(string(baseRows{r,5}),phaseOrder)}]; %#ok<AGROW>
        end
    end
end

modelStatsT = cell2table(modelStatRows,'VariableNames', ...
    {'cell_type','period','metric','prevAction','term','statisticName', ...
     'statistic','DF1','DF2','pValue','nTrials','nAnimals'});
modelStatsT = normalizeMixedTable(modelStatsT,1:6,7:12);

coefficientT = cell2table(coefficientRows,'VariableNames', ...
    {'cell_type','period','metric','prevAction','coefficient','estimate', ...
     'SE','statistic','DF','pValue','lowerCI','upperCI','nTrials','phase'});
coefficientT = normalizeMixedTable(coefficientT,[1:5 14],6:13);
end


function fig = plotFigureS13(acuteAnimal, acuteStats, sustainedAnimal, ...
        sustainedGroup, sustainedStats, sustainedCoef, phaseOrder, ...
        cellTypes, conditions, actions)
fig = figure('Color','w','Position',[80 60 1400 900]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

metricByColumn = ["winStay","winStay","loseSwitch","loseSwitch"];
actionByColumn = ["Push","Pull","Push","Pull"];
panelIndex = 0;

for ct = 1:numel(cellTypes)
    for col = 1:4
        panelIndex = panelIndex + 1;
        ax = nexttile;
        plotAcutePanel(ax,acuteAnimal,acuteStats,cellTypes(ct), ...
            metricByColumn(col),actionByColumn(col),conditions);
        addPanelLetter(ax,panelIndex);
        if col == 1
            text(ax,-0.27,0.5,upper(cellTypes(ct)), ...
                'Units','normalized','Rotation',90, ...
                'HorizontalAlignment','center','FontWeight','bold');
        end
    end
end

bottomCellType = ["sst","sst","th","th"];
bottomAction = ["Push","Pull","Push","Pull"];
for col = 1:4
    panelIndex = panelIndex + 1;
    ax = nexttile;
    plotSustainedPanel(ax,sustainedAnimal,sustainedGroup,sustainedStats, ...
        sustainedCoef,bottomCellType(col),bottomAction(col), ...
        phaseOrder,conditions);
    addPanelLetter(ax,panelIndex);
    if col == 1 || col == 3
        text(ax,-0.27,0.5,upper(bottomCellType(col)), ...
            'Units','normalized','Rotation',90, ...
            'HorizontalAlignment','center','FontWeight','bold');
    end
end
end


function plotAcutePanel(ax,animalT,statsT,cellType,metric,action,conditions)
colors = struct('ctrl',[0.68 0.68 0.68], 'ib',[1.00 0.42 0.46]);
lineColors = struct('ctrl',[0.25 0.25 0.25], 'ib',[1.00 0.28 0.32]);
preX = [0.86 1.14];
cluster = [1.90 2.55];
hold(ax,'on');

for c = 1:numel(conditions)
    cond = conditions(c);
    P = animalT(animalT.cell_type == cellType & ...
        animalT.metric == metric & animalT.prevAction == action & ...
        animalT.condition == cond & animalT.period == "PRE",:);
    drawBarAndPoints(ax,preX(c),P.value,colors.(cond),0.35);

    O = animalT(animalT.cell_type == cellType & ...
        animalT.metric == metric & animalT.prevAction == action & ...
        animalT.condition == cond & animalT.period == "OPTO",:);
    animals = unique(O.animalID);
    xOff = cluster(c)-0.14;
    xOn = cluster(c)+0.14;
    off = nan(numel(animals),1);
    on = nan(numel(animals),1);
    for i = 1:numel(animals)
        offRow = O.animalID == animals(i) & O.light == "Light Off";
        onRow = O.animalID == animals(i) & O.light == "Light On";
        if any(offRow)
            off(i) = O.value(find(offRow,1));
        end
        if any(onRow)
            on(i) = O.value(find(onRow,1));
        end
    end
    complete = isfinite(off) & isfinite(on);
    off = off(complete);
    on = on(complete);
    drawBarAndPoints(ax,xOff,off,colors.(cond),0.35);
    drawBarAndPoints(ax,xOn,on,colors.(cond),0.85);
    for i = 1:numel(off)
        plot(ax,[xOff xOn],[off(i) on(i)],'-', ...
            'Color',lineColors.(cond),'LineWidth',0.55, ...
            'HandleVisibility','off');
    end
end

if metric == "winStay"
    ylim(ax,[0.4 1]);
    ylabel(ax,'P(stay | win)');
    titleRoot = "Win-stay";
else
    ylim(ax,[0 1]);
    ylabel(ax,'P(switch | loss)');
    titleRoot = "Lose-switch";
end
xlim(ax,[0.58 2.86]);
xticks(ax,[1 cluster]);
xticklabels(ax,["PRE","Control","Inhibition"]);
xtickangle(ax,0);
yl = ylim(ax);
yLight = yl(1) + 0.025*range(yl);
text(ax,cluster(1)-0.14,yLight,'OFF','HorizontalAlignment','center', ...
    'FontSize',6.5,'Color',[0.35 0.35 0.35]);
text(ax,cluster(1)+0.14,yLight,'ON','HorizontalAlignment','center', ...
    'FontSize',6.5,'Color',[0.35 0.35 0.35]);
text(ax,cluster(2)-0.14,yLight,'OFF','HorizontalAlignment','center', ...
    'FontSize',6.5,'Color',[0.70 0.10 0.12]);
text(ax,cluster(2)+0.14,yLight,'ON','HorizontalAlignment','center', ...
    'FontSize',6.5,'Color',[0.70 0.10 0.12]);
title(ax,titleRoot + " | " + action);
box(ax,'off');
set(ax,'FontName','Arial','FontSize',8.5,'LineWidth',0.8);

preP = findTermP(statsT,cellType,"PRE",metric,action,"condition");
optoP = findTermP(statsT,cellType,"OPTO",metric,action,"condition");
yText = yl(2)-0.035*range(yl);
text(ax,1,yText,pLabel(preP),'HorizontalAlignment','center', ...
    'FontSize',8,'Interpreter','none');
text(ax,mean(cluster),yText,pLabel(optoP),'HorizontalAlignment','center', ...
    'FontSize',8,'FontWeight',significanceWeight(optoP), ...
    'Interpreter','none');
end


function plotSustainedPanel(ax,animalT,groupT,statsT,coefT, ...
        cellType,action,phaseOrder,conditions)
colors = struct('ctrl',[0 0 0], 'ib',[0.93 0.10 0.12]);
lineColors = struct('ctrl',[0.76 0.76 0.76], 'ib',[1.00 0.76 0.76]);
hold(ax,'on');

for c = 1:numel(conditions)
    cond = conditions(c);
    animals = unique(animalT.animalID(animalT.cell_type == cellType & ...
        animalT.prevAction == action & animalT.condition == cond));
    for i = 1:numel(animals)
        A = animalT(animalT.cell_type == cellType & ...
            animalT.prevAction == action & animalT.condition == cond & ...
            animalT.animalID == animals(i),:);
        [tf,x] = ismember(A.phase,phaseOrder);
        y = A.pLoseSwitch(tf);
        x = x(tf);
        [x,ord] = sort(x);
        y = y(ord);
        if numel(x) >= 2
            plot(ax,x,y,'-','Color',lineColors.(cond), ...
                'LineWidth',0.55,'HandleVisibility','off');
        end
    end

    G = groupT(groupT.cell_type == cellType & ...
        groupT.prevAction == action & groupT.condition == cond,:);
    [tf,x] = ismember(G.phase,phaseOrder);
    G = G(tf,:);
    x = x(tf);
    [x,ord] = sort(x);
    G = G(ord,:);
    errorbar(ax,x,G.mean_pLoseSwitch,G.sem_pLoseSwitch,'-o', ...
        'Color',colors.(cond),'MarkerFaceColor',colors.(cond), ...
        'MarkerEdgeColor','none','MarkerSize',4,'LineWidth',0.9, ...
        'CapSize',3);
end

xlim(ax,[0.7 5.3]);
ylim(ax,[0 1]);
xticks(ax,1:5);
xticklabels(ax,["Pre","Opto I","Opto II","Post I","Post II"]);
xtickangle(ax,45);
ylabel(ax,'P(switch | loss)');
if action == "Push"
    title(ax,'Push \rightarrow switch','Color',[0 0.65 0.70]);
else
    title(ax,'Pull \rightarrow switch','Color',[0.62 0.15 0.62]);
end
box(ax,'off');
set(ax,'FontName','Arial','FontSize',8.5,'LineWidth',0.8);

pInteraction = findTermP(statsT,cellType,"ALL","loseSwitch", ...
    action,"condition:phase");
text(ax,0.5,0.06,"p(interaction) = " + compactP(pInteraction), ...
    'Units','normalized','HorizontalAlignment','center', ...
    'FontSize',8,'Interpreter','none');

pByPhase = phaseContrastPValues(coefT,cellType,action,phaseOrder);
for p = 1:numel(phaseOrder)
    if isfinite(pByPhase(p)) && pByPhase(p) < 0.05
        text(ax,p,0.94,'*','HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',13);
    end
end
end


function drawBarAndPoints(ax,x,values,color,alphaValue)
values = values(isfinite(values));
bar(ax,x,mean(values,'omitnan'),0.24,'FaceColor',color, ...
    'FaceAlpha',alphaValue,'EdgeColor','none');
if isempty(values)
    return;
end
offset = zeros(size(values));
if numel(values) > 1
    offset = linspace(-0.035,0.035,numel(values))';
end
scatter(ax,x+offset,values,13,'MarkerFaceColor',color, ...
    'MarkerEdgeColor','none');
end


function addPanelLetter(ax,index)
letter = char(double('a') + index - 1);
text(ax,-0.16,1.07,letter,'Units','normalized', ...
    'FontWeight','bold','FontSize',12);
end


function p = findTermP(statsT,cellType,period,metric,action,targetTerm)
idx = statsT.cell_type == cellType & statsT.period == period & ...
      statsT.metric == metric & statsT.prevAction == action;
S = statsT(idx,:);
terms = lower(erase(S.term," "));
target = lower(erase(string(targetTerm)," "));
row = find(terms == target,1);
if isempty(row) && contains(target,":")
    parts = split(target,":");
    row = find(contains(terms,parts(1)) & contains(terms,parts(2)) & ...
        contains(terms,":"),1);
end
if isempty(row)
    p = NaN;
else
    p = S.pValue(row);
end
end


function p = phaseContrastPValues(coefT,cellType,action,phaseOrder)
p = nan(size(phaseOrder));
S = coefT(coefT.cell_type == cellType & ...
    coefT.prevAction == action,:);
for k = 1:numel(phaseOrder)
    idx = S.phase == phaseOrder(k);
    if sum(idx) == 1
        p(k) = S.pValue(idx);
    end
end
end


function label = pLabel(p)
label = "p = " + compactP(p);
if isfinite(p) && p < 0.05
    label = label + "*";
end
end


function label = compactP(p)
if ~isfinite(p)
    label = "NA";
elseif p < 0.0001
    label = "<0.0001";
elseif p < 0.01
    label = string(sprintf('%.4f',p));
elseif p < 0.995
    label = string(sprintf('%.2f',p));
else
    label = "1.0";
end
end


function weight = significanceWeight(p)
weight = 'normal';
if isfinite(p) && p < 0.05
    weight = 'bold';
end
end


function rows = appendAnovaRows(rows,stats,cellType,period,metric, ...
        action,nTrials,nAnimals)
vars = statVariableNames(stats);
termVar = findVariable(vars,'Term');
pVar = findVariable(vars,'pValue');
fVar = findVariable(vars,'FStat');
chiVar = findVariable(vars,'Chi2Stat');
df1Var = findVariable(vars,'DF1');
df2Var = findVariable(vars,'DF2');

terms = string(stats.(termVar));
for r = 1:numel(terms)
    statisticName = "FStat";
    statistic = NaN;
    if ~isempty(fVar)
        values = stats.(fVar);
        statistic = values(r);
    elseif ~isempty(chiVar)
        statisticName = "Chi2Stat";
        values = stats.(chiVar);
        statistic = values(r);
    end
    pValues = stats.(pVar);
    df1 = NaN;
    df2 = NaN;
    if ~isempty(df1Var)
        values = stats.(df1Var);
        df1 = values(r);
    end
    if ~isempty(df2Var)
        values = stats.(df2Var);
        df2 = values(r);
    end
    rows(end+1,:) = {cellType,period,metric,action,terms(r), ...
        statisticName,statistic,df1,df2,pValues(r),nTrials,nAnimals}; %#ok<AGROW>
end
end


function rows = appendCoefficientRows(rows,mdl,cellType,period,metric, ...
        action,nTrials)
C = mdl.Coefficients;
names = coefficientNames(mdl,C);
for r = 1:numel(names)
    rows(end+1,:) = {cellType,period,metric,action,names(r), ...
        statValue(C,r,'Estimate'),statValue(C,r,'SE'), ...
        firstAvailableStat(C,r),statValue(C,r,'DF'), ...
        statValue(C,r,'pValue'),statValue(C,r,'Lower'), ...
        statValue(C,r,'Upper'),nTrials}; %#ok<AGROW>
end
end


function rows = lightSimpleEffectRows(mdl,cellType,metric,action)
coefNames = string(mdl.CoefficientNames);
idxLight = find(contains(coefNames,"prevLight") & ...
    ~contains(coefNames,":"));
idxInteraction = find(contains(coefNames,"condition") & ...
    contains(coefNames,"prevLight") & contains(coefNames,":"));
if numel(idxLight) ~= 1 || numel(idxInteraction) ~= 1
    error('FigureS13:CoefficientIdentification', ...
        'Could not identify light and condition-by-light coefficients.');
end

nCoef = mdl.NumCoefficients;
Hctrl = zeros(1,nCoef);
Hctrl(idxLight) = 1;
Hib = Hctrl;
Hib(idxInteraction) = 1;

[pCtrl,FCtrl,DF1Ctrl,DF2Ctrl] = coefTest(mdl,Hctrl,0);
[pIb,FIb,DF1Ib,DF2Ib] = coefTest(mdl,Hib,0);
beta = fixedEffects(mdl);
logOdds = [Hctrl*beta; Hib*beta];

rows = {cellType,metric,action,"ctrl",logOdds(1),exp(logOdds(1)), ...
            FCtrl,DF1Ctrl,DF2Ctrl,pCtrl; ...
        cellType,metric,action,"ib",logOdds(2),exp(logOdds(2)), ...
            FIb,DF1Ib,DF2Ib,pIb};
end


function phase = coefficientPhase(name,phaseOrder)
name = string(name);
phase = "";
if ~contains(lower(name),"condition_ib")
    return;
end
if strcmpi(name,"condition_ib")
    phase = "PRE";
    return;
end
for k = 2:numel(phaseOrder)
    if contains(upper(name),"PHASE_" + phaseOrder(k))
        phase = phaseOrder(k);
        return;
    end
end
end


function names = coefficientNames(mdl,C)
try
    names = string(mdl.CoefficientNames(:));
catch
    names = string(C.Properties.RowNames);
end
if isempty(names)
    names = "coefficient_" + string((1:height(C))');
end
end


function value = firstAvailableStat(stats,row)
value = statValue(stats,row,'tStat');
if ~isfinite(value)
    value = statValue(stats,row,'zStat');
end
end


function value = statValue(stats,row,variable)
vars = statVariableNames(stats);
idx = find(strcmpi(vars,variable),1);
if isempty(idx)
    value = NaN;
else
    values = stats.(vars{idx});
    value = values(row);
end
end


function variable = findVariable(vars,name)
idx = find(strcmpi(vars,name),1);
if isempty(idx)
    variable = [];
else
    variable = vars{idx};
end
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
error('FigureS13:UnknownStatsContainer', ...
    'Could not determine variable names in statistical output.');
end


function T = normalizeMixedTable(T,stringColumns,numericColumns)
for i = stringColumns
    variable = T.Properties.VariableNames{i};
    T.(variable) = string(T.(variable));
end
for i = numericColumns
    variable = T.Properties.VariableNames{i};
    T.(variable) = unpackScalarColumn(T.(variable));
end
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


function exportFigurePair(fig,outputFolder,stem)
pngFile = fullfile(outputFolder,stem + ".png");
pdfFile = fullfile(outputFolder,stem + ".pdf");
if exist('exportgraphics','file') == 2
    exportgraphics(fig,pngFile,'Resolution',300);
    exportgraphics(fig,pdfFile,'ContentType','vector');
else
    set(fig,'PaperPositionMode','auto');
    print(fig,pngFile,'-dpng','-r300');
    print(fig,pdfFile,'-dpdf','-painters','-bestfit');
end
end


function assertFileExists(file)
if ~isfile(file)
    error('FigureS13:FileNotFound','File not found: %s',file);
end
end


function requireVariables(T,required,label)
missing = setdiff(required,string(T.Properties.VariableNames));
if ~isempty(missing)
    error('FigureS13:MissingVariables', ...
        '%s is missing required variables: %s', ...
        label,strjoin(missing,', '));
end
end


function assertUniqueKeys(T,keyVars,label)
keyT = T(:,keyVars);
if height(unique(keyT,'rows')) ~= height(keyT)
    error('FigureS13:DuplicateKeys', ...
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
