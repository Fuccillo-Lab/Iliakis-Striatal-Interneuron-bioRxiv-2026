function results = plot_figure4_behavior( ...
        trialFile, metadataFile, posteriorFile, animalFile, outputFolder)
%PLOT_FIGURE4_BEHAVIOR Reproduce the behavioral data in Figure 4C-F.
%
% RESULTS = PLOT_FIGURE4_BEHAVIOR(TRIALFILE, METADATAFILE, POSTERIORFILE,
% ANIMALFILE, OUTPUTFOLDER) plots animal-level P(pull on next trial) after
% rewarded pushes (C,E) and rewarded pulls (D,F), for SST and TH cohorts.
% The optional OUTPUTFOLDER saves the figure and CSV source/statistics tables.
%
% Inputs
% ------
% TRIALFILE: model-ready table from prepare_rlhmm_table.m. Requires
% animalID, iOrig, jOrig, Y, currChoice, currReward, prevChoice, prevReward.
% METADATAFILE: trial-keyed companion table from the same importer. Requires
% animalID, iOrig, jOrig, assignment to identify Pre/Opto/Post phases.
% POSTERIORFILE: rlhmm_posteriors.csv from the selected final fit. Only its
% animalID/iOrig/jOrig/Y keys are used, to retain the exact trial set from
% the historical inner join; state probabilities are not used in Figure 4.
% ANIMALFILE: animal roster with animal_id, cell_type, condition. The plotted
% cohorts are sst/th and their conditions must be labeled ctrl/ib.
%
% Historical choices
% ------------------
% The plotted probabilities first average trials within animal x phase,
% then give animals equal weight in each group mean/SEM. The plotted
% push-win-to-pull probability equals 1 minus the original push win-stay;
% pull-win-to-pull equals the original pull win-stay. For the interaction
% tests, the original script independently rebuilt previous-trial variables
% on adjacent retained rows with consecutive jOrig and fit binomial GLMEs
% to win-stay (condition * bin + random intercepts for animal and session).
% Here those same tests model pull on the next trial directly. Complementing
% the binary push-win response reverses coefficient signs but leaves the
% interaction p-value unchanged; the pull-win response is already pull.
% This intentional difference in pair selection is preserved. Significance
% stars mark p < 0.05 for the fitted condition coefficient at PRE or the
% condition-by-bin interaction coefficient at each later bin; the latter
% tests a change in the group contrast relative to PRE.
%
% Requires MATLAB's Statistics and Machine Learning Toolbox (fitglme).

    if nargin < 5
        outputFolder = '';
    end
    inputs = {trialFile, metadataFile, posteriorFile, animalFile};
    for i = 1:numel(inputs)
        if ~isfile(inputs{i})
            error('Input file does not exist: %s', inputs{i});
        end
    end

    B = readtable(trialFile, 'TextType', 'string');
    M = readtable(metadataFile, 'TextType', 'string');
    P = readtable(posteriorFile, 'TextType', 'string');
    A = readtable(animalFile, 'TextType', 'string');
    key = {'animalID','iOrig','jOrig'};
    requireVariables(B, [key, {'Y','currChoice','currReward','prevChoice','prevReward'}], 'trial table');
    requireVariables(M, [key, {'assignment'}], 'plotting metadata');
    requireVariables(P, [key, {'Y'}], 'posterior table');
    requireVariables(A, {'animal_id','cell_type','condition'}, 'animal roster');
    assertUnique(B, key, 'trial table');
    assertUnique(M, key, 'plotting metadata');
    assertUnique(P, key, 'posterior table');
    assertUnique(A, {'animal_id'}, 'animal roster');

    % The original analysis joined modelTfinal to posterior output by trial
    % key AND choice, then joined animal metadata. Keep the same selection.
    T = innerjoin(B, M(:, [key, {'assignment'}]), 'Keys', key);
    if height(T) ~= height(B)
        error('Plotting metadata must cover every model-ready trial exactly once.');
    end
    nBeforePosterior = height(T);
    T = innerjoin(T, P(:, [key, {'Y'}]), 'Keys', [key, {'Y'}]);
    fprintf('Figure 4C-F: %d/%d model-ready trials matched posterior keys.\n', ...
        height(T), nBeforePosterior);
    roster = A(:, {'animal_id','cell_type','condition'});
    roster = renamevars(roster, 'animal_id', 'animalID');
    T = innerjoin(T, roster, 'Keys', 'animalID');
    T = sortrows(T, key);

    binOrder = ["PRE","OPTO1","OPTO2","POST1","POST2"];
    T.bin = strings(height(T),1);
    assignment = string(T.assignment);
    T.bin(ismember(assignment, ["PRE1","PRE2","PRE"])) = "PRE";
    T.bin(ismember(assignment, ["O1","O2","O3","OPTO1"])) = "OPTO1";
    T.bin(ismember(assignment, ["O4","O5","O6","bIntra-Opto"])) = "OPTO2";
    T.bin(ismember(assignment, ["POST1","POST2"])) = "POST1";
    T.bin(ismember(assignment, ["POST3","POST4"])) = "POST2";
    T.bin = categorical(T.bin, binOrder, 'Ordinal', true);
    T.cell_type = lower(strtrim(string(T.cell_type)));
    T.condition = lower(strtrim(string(T.condition)));
    T = T(~isundefined(T.bin) & ...
        ismember(T.cell_type, ["sst","th"]) & ...
        ismember(T.condition, ["ctrl","ib"]), :);
    fprintf('Figure 4C-F: %d selected trials, %d animals after phase/cohort filters.\n', ...
        height(T), numel(unique(T.animalID)));

    valid = ismember(T.prevChoice, [-1,1]) & ...
        ismember(T.currChoice, [-1,1]) & ismember(T.prevReward, [0,1]);
    plotT = T(valid,:);
    plotT.pushWinToPull = nan(height(plotT),1);
    plotT.pullWinToPull = nan(height(plotT),1);
    pushWin = plotT.prevChoice == 1 & plotT.prevReward == 1;
    pullWin = plotT.prevChoice == -1 & plotT.prevReward == 1;
    plotT.pushWinToPull(pushWin) = double(plotT.currChoice(pushWin) == -1);
    plotT.pullWinToPull(pullWin) = double(plotT.currChoice(pullWin) == -1);

    [G, animalID, cellType, bin, condition] = findgroups( ...
        plotT.animalID, plotT.cell_type, plotT.bin, plotT.condition);
    animalBin = table(animalID, cellType, bin, condition);
    measures = {'pushWinToPull','pullWinToPull'};
    for m = 1:numel(measures)
        value = plotT.(measures{m});
        animalBin.(measures{m}) = splitapply(@(x) mean(x,'omitnan'), value, G);
        animalBin.(['n_' measures{m}]) = splitapply(@(x) sum(~isnan(x)), value, G);
    end

    specs = {
        'C', 'sst', 'pushWinToPull', 'Push win to pull';
        'D', 'sst', 'pullWinToPull', 'Pull win to pull';
        'E', 'th',  'pushWinToPull', 'Push win to pull';
        'F', 'th',  'pullWinToPull', 'Pull win to pull'};
    groupRows = table();
    statRows = table();
    coefficientRows = table();
    fig = figure('Color','w', 'Name','Figure 4C-F: optogenetic behavior', ...
        'Position',[80 200 1300 440]);
    layout = tiledlayout(fig,1,4,'TileSpacing','compact','Padding','compact');

    for k = 1:size(specs,1)
        panel = string(specs{k,1});
        cellTypeName = string(specs{k,2});
        measure = char(specs{k,3});
        S = animalBin(animalBin.cellType == cellTypeName,:);
        ax = nexttile(layout); hold(ax,'on');
        for c = ["ctrl","ib"]
            col = conditionColor(c);
            faint = 0.20 .* col + 0.80 .* [1 1 1];
            theseAnimals = unique(S.animalID(S.condition == c));
            for i = 1:numel(theseAnimals)
                R = sortrows(S(S.animalID == theseAnimals(i),:), 'bin');
                plot(ax,double(R.bin),R.(measure),'-','Color',faint,'LineWidth',0.8);
            end
            for b = 1:numel(binOrder)
                vals = S.(measure)(S.condition == c & string(S.bin) == binOrder(b));
                vals = vals(isfinite(vals));
                row = table(panel,cellTypeName,c,binOrder(b),numel(vals), ...
                    mean(vals,'omitnan'),std(vals,'omitnan')/sqrt(numel(vals)), ...
                    'VariableNames',{'panel','cell_type','condition','bin', ...
                    'nAnimals','meanPpull','semPpull'});
                groupRows = [groupRows;row]; %#ok<AGROW>
            end
            R = groupRows(groupRows.panel == panel & groupRows.condition == c,:);
            errorbar(ax,1:5,R.meanPpull,R.semPpull,'-o','Color',col, ...
                'MarkerFaceColor',col,'LineWidth',2,'CapSize',6);
        end
        if ismember(panel,["C","E"])
            ylim(ax,[0 0.6]);
        else
            ylim(ax,[0.4 1]);
        end
        xlim(ax,[0.8 5.2]);
        xticks(ax,1:5);
        xticklabels(ax,["Pre","Opto I","Opto II","Post I","Post II"]);
        xtickangle(ax,40);
        ylabel(ax,'P(pull on next trial)');
        title(ax,sprintf('%s  %s',panel,upper(cellTypeName)));
        box(ax,'off');

        [pInteraction,nTrials,nAnimals,coefficients] = interactionTest( ...
            T,cellTypeName,measure,binOrder);
        fprintf('Panel %s: %d consecutive win trials, %d animals, interaction p = %.6g.\n', ...
            panel,nTrials,nAnimals,pInteraction);
        statRows = [statRows; table(panel,cellTypeName,string(measure), ...
            nTrials,nAnimals,pInteraction, 'VariableNames', ...
            {'panel','cell_type','measure','nConsecutiveTrials', ...
            'nAnimals','pInteraction'})]; %#ok<AGROW>
        coefficients.panel = repmat(panel,height(coefficients),1);
        coefficients.cell_type = repmat(cellTypeName,height(coefficients),1);
        coefficients = movevars(coefficients,{'panel','cell_type'},'Before',1);
        coefficientRows = [coefficientRows;coefficients]; %#ok<AGROW>
        limits = ylim(ax);
        for b = 1:height(coefficients)
            if ~coefficients.significant(b), continue; end
            meanAtBin = groupRows.meanPpull(groupRows.panel == panel & ...
                groupRows.bin == binOrder(b));
            starY = min(limits(2) - 0.025 * diff(limits), ...
                max(meanAtBin) + 0.08 * diff(limits));
            text(ax,b,starY,'*','HorizontalAlignment','center', ...
                'FontSize',16,'FontWeight','bold');
        end
        text(ax,0.04,0.08,sprintf('p_{interaction} = %.3g',pInteraction), ...
            'Units','normalized','FontSize',9);
    end

    results = struct('figure',fig,'animal_bin',animalBin, ...
        'group_summary',groupRows,'interaction_stats',statRows, ...
        'coefficient_stats',coefficientRows);
    if ~isempty(outputFolder)
        if ~isfolder(outputFolder), mkdir(outputFolder); end
        exportgraphics(fig,fullfile(outputFolder,'figure4_panelsC_to_F.png'),'Resolution',300);
        exportgraphics(fig,fullfile(outputFolder,'figure4_panelsC_to_F.pdf'),'ContentType','vector');
        writetable(animalBin,fullfile(outputFolder,'figure4_panelsC_to_F_animal_bin.csv'));
        writetable(groupRows,fullfile(outputFolder,'figure4_panelsC_to_F_group_summary.csv'));
        writetable(statRows,fullfile(outputFolder,'figure4_panelsC_to_F_interaction_stats.csv'));
        writetable(coefficientRows,fullfile(outputFolder,'figure4_panelsC_to_F_coefficient_stats.csv'));
    end
end


function [p,nTrials,nAnimals,coefficientRows] = interactionTest( ...
        T,cellTypeName,measure,binOrder)
% Fit the source's strictly consecutive-trial GLME with a pull response.
    S = T(T.cell_type == cellTypeName,:);
    S = sortrows(S,{'animalID','iOrig','jOrig'});
    n = height(S);
    prevChoice = nan(n,1);
    prevReward = nan(n,1);
    pull = nan(n,1);
    for i = 2:n
        sameSession = S.animalID(i) == S.animalID(i-1) && S.iOrig(i) == S.iOrig(i-1);
        consecutive = S.jOrig(i) == S.jOrig(i-1) + 1;
        if sameSession && consecutive && ...
                ismember(S.currChoice(i),[-1,1]) && ...
                ismember(S.currChoice(i-1),[-1,1]) && ...
                ismember(S.currReward(i-1),[0,1])
            prevChoice(i) = S.currChoice(i-1);
            prevReward(i) = S.currReward(i-1);
            pull(i) = double(S.currChoice(i) == -1);
        end
    end
    action = 1;
    if strcmp(measure,'pullWinToPull'), action = -1; end
    keep = prevChoice == action & prevReward == 1 & isfinite(pull);
    S = S(keep,:);
    if isempty(S), error('No consecutive rewarded %s trials for %s.',measure,cellTypeName); end
    L = table();
    L.animalID = categorical(S.animalID);
    L.sessionID = categorical(string(S.animalID) + "_" + string(S.iOrig));
    L.bin = categorical(string(S.bin),binOrder,'Ordinal',true);
    L.condition = categorical(string(S.condition),["ctrl","ib"]);
    L.pull = pull(keep);
    mdl = fitglme(L,'pull ~ condition * bin + (1|animalID) + (1|sessionID)', ...
        'Distribution','Binomial','Link','logit', ...
        'DummyVarCoding','reference');
    stats = anova(mdl);
    if ismember('Term',stats.Properties.VarNames)
        terms = string(stats.Term);
    else
        terms = string(stats.Properties.RowNames);
    end
    row = contains(lower(terms),'condition') & contains(lower(terms),'bin');
    if nnz(row) ~= 1 || ~ismember('pValue',stats.Properties.VarNames)
        error('Cannot identify the condition-by-bin interaction in GLME ANOVA.');
    end
    p = stats.pValue(row);
    nTrials = height(L);
    nAnimals = numel(unique(S.animalID));

    % Under MATLAB's reference coding (ctrl, PRE), the condition coefficient
    % is the group contrast at PRE. A condition:bin coefficient tests the
    % CHANGE in that contrast from PRE, as in the supplied source model.
    fixed = mdl.Coefficients;
    names = string(fixed.Name);
    coefficientRows = table();
    for b = 1:numel(binOrder)
        if b == 1
            expected = "condition_ib";
            comparison = "group contrast at PRE";
        else
            expected = "bin_" + binOrder(b);
            comparison = "interaction vs PRE";
        end
        if b == 1
            match = names == expected;
        else
            match = names == ("condition_ib:" + expected) | ...
                names == (expected + ":condition_ib");
        end
        if nnz(match) ~= 1
            error('Cannot identify the %s coefficient in GLME output.',binOrder(b));
        end
        bin = binOrder(b);
        coefficient = names(match);
        estimate = fixed.Estimate(match);
        pCoefficient = fixed.pValue(match);
        significant = pCoefficient < 0.05;
        coefficientRows = [coefficientRows; table(bin,coefficient,comparison, ...
            estimate,pCoefficient,significant)]; %#ok<AGROW>
    end
end


function col = conditionColor(condition)
    if condition == "ctrl", col = [0 0 0]; else, col = [1 0 0]; end
end


function requireVariables(T, variables, description)
    missing = setdiff(variables,T.Properties.VariableNames);
    if ~isempty(missing)
        error('%s missing required columns: %s',description,strjoin(missing,', '));
    end
end


function assertUnique(T, variables, description)
    if height(unique(T(:,variables),'rows')) ~= height(T)
        error('%s has duplicate key rows: %s',description,strjoin(variables,', '));
    end
end
