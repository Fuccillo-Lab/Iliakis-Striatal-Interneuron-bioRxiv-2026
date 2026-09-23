function results = plot_rlhmm_figure5( ...
        plottingMetadataFile, posteriorFile, parameterFile, ...
        transitionFile, animalFile, outputFolder)
%PLOT_RLHMM_FIGURE5 Reproduce the RL-GLM-HMM analyses in Figure 5.
%
% RESULTS = PLOT_RLHMM_FIGURE5(PLOTTINGMETADATAFILE, POSTERIORFILE,
% PARAMETERFILE, TRANSITIONFILE, ANIMALFILE, OUTPUTFOLDER) reads the final
% selected fit and produces clean source figures and statistics tables for
% Figure 5A-O.
%
% Required inputs
% ---------------
% PLOTTINGMETADATAFILE
%   Companion CSV written by prepare_rlhmm_table.m. It contains the trial
%   keys plus assignment, highSide, and the historical precomputed deltaQ.
% POSTERIORFILE
%   rlhmm_posteriors.csv from the selected fit.
% PARAMETERFILE
%   rl_params.csv from the selected fit.
% TRANSITIONFILE
%   transition_matrix.csv from the selected fit.
% ANIMALFILE
%   Animal roster with animal_id, cell_type, condition, and sex.
% OUTPUTFOLDER
%   Optional. When supplied, figures are exported as PNG and PDF and the
%   analysis summaries are written as CSV files.
%
% Historical-reproduction decisions
% -----------------------------------
% 1. Panel F deliberately uses the legacy precomputed deltaQ found in the
%    companion plotting table. This reproduces the submitted paper; it is
%    not a state-specific Q trace reconstructed from the final model.
% 2. Panels D, E, G, H, and I use hard maximum-a-posteriori (MAP) state
%    assignments. Panel F's population curves and panels J-O use posterior
%    state probabilities.
% 3. In panel B, beta_pos and beta_neg are multiplied by 8 to display the
%    fitted effect across an 8-unit value difference. The saved statistics
%    table contains both raw and display-scaled coefficients.
% 4. Panel F preserves the submitted convention: population curves are
%    posterior-weighted, while thin animal curves use hard MAP states.
% 5. The panel G/H familywise correction includes all eight planned tests:
%    four stay/switch measures by two state contrasts.
%
% Software
% --------
% Base MATLAB is sufficient for reading, summarizing, and plotting. The
% Statistics and Machine Learning Toolbox is required for signrank and
% fitlme. No Python environment or upstream GLM-HMM source is required once
% the final fit CSV files have been generated.

    if nargin < 6
        outputFolder = '';
    end

    inputFiles = {plottingMetadataFile, posteriorFile, parameterFile, ...
        transitionFile, animalFile};
    for iFile = 1:numel(inputFiles)
        if ~isfile(inputFiles{iFile})
            error('Input file does not exist: %s', inputFiles{iFile});
        end
    end

    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    colors = [ ...
        75 105 177; ...   % State 1: #4b69b1
        60 192 196; ...   % State 2: #3cc0c4
        153 57 149] ./ 255; % State 3: #993995
    states = (1:3)';

    %% Read and validate selected-fit artifacts
    M = readtable(plottingMetadataFile);
    P = readtable(posteriorFile);
    params = sortrows(readtable(parameterFile), 'state');
    transition = readTransitionMatrix(transitionFile);
    animals = readtable(animalFile);

    trialKeys = {'animalID','iOrig','jOrig'};
    outcomeKeys = {'currChoice','currReward','Y'};
    probabilityVars = {'stateProb_1','stateProb_2','stateProb_3'};

    requireVariables(M, [trialKeys, outcomeKeys, ...
        {'assignment','highSide','deltaQ'}], 'plotting metadata');
    requireVariables(P, [trialKeys, outcomeKeys, probabilityVars, ...
        {'stateMAP'}], 'posterior table');
    requireVariables(params, {'state','alpha','bias','beta_pos','beta_neg', ...
        'prevChoice_prevReward','prevChoice_prevUnreward','animalBias', ...
        'recentChoiceFrac_z','prevChoice_switchFrac_z'}, 'parameter table');
    requireVariables(animals, {'animal_id','cell_type','condition','sex'}, ...
        'animal roster');

    if ~isequal(size(transition), [3 3]) || ...
            any(~isfinite(transition), 'all') || ...
            any(abs(sum(transition, 2) - 1) > 1e-8)
        error('Transition matrix must be finite, 3-by-3, and row-stochastic.');
    end
    if height(params) ~= 3 || ~isequal(params.state(:), states)
        error('Parameter table must contain exactly States 1, 2, and 3.');
    end

    assertUniqueTrialKeys(M, trialKeys, 'plotting metadata');
    assertUniqueTrialKeys(P, trialKeys, 'posterior table');

    M = sortrows(M, trialKeys);
    P = sortrows(P, trialKeys);
    if height(M) ~= height(P) || ...
            ~isequal(M{:, trialKeys}, P{:, trialKeys})
        error('Plotting metadata and posterior rows do not have identical trial keys.');
    end
    for iVar = 1:numel(outcomeKeys)
        variable = outcomeKeys{iVar};
        if ~isequaln(M.(variable), P.(variable))
            error('Metadata/posterior mismatch in variable %s.', variable);
        end
    end

    posteriorSum = sum(P{:, probabilityVars}, 2);
    if any(abs(posteriorSum - 1) > 1e-8)
        error('Posterior state probabilities do not sum to one.');
    end
    if any(~ismember(P.stateMAP, states))
        error('stateMAP contains values outside States 1-3.');
    end

    T = [M, P(:, [probabilityVars, {'stateMAP'}])];
    T.assignment = string(T.assignment);
    fprintf('Validated %d trials, %d sessions, and %d animals.\n', ...
        height(T), height(unique(T(:, {'animalID','iOrig'}), 'rows')), ...
        numel(unique(T.animalID)));

    %% Panel A: selected transition matrix
    figA = figure('Color', 'w', 'Name', 'Figure 5A: transition matrix');
    imagesc(transition, [0 1]);
    axis square;
    colormap(figA, flipud(gray(256)));
    colorbar;
    xticks(1:3); yticks(1:3);
    xticklabels(compose('To state %d', states));
    yticklabels(compose('From state %d', states));
    xlabel('Next state'); ylabel('Current state');
    title('Selected-model transition probabilities');
    for row = 1:3
        for col = 1:3
            text(col, row, sprintf('%.4f', transition(row,col)), ...
                'HorizontalAlignment', 'center', ...
                'Color', contrastTextColor(transition(row,col)), ...
                'FontWeight', 'bold');
        end
    end
    styleAxes(gca);

    transitionTable = array2table(transition, ...
        'VariableNames', compose('to_state_%d', states));
    transitionTable.from_state = states;
    transitionTable = movevars(transitionTable, 'from_state', 'Before', 1);

    %% Panel B: fitted weights
    featureVars = [ ...
        "bias"
        "animalBias"
        "beta_pos"
        "beta_neg"
        "prevChoice_prevReward"
        "prevChoice_prevUnreward"
        "recentChoiceFrac_z"
        "prevChoice_switchFrac_z"];
    featureLabels = [ ...
        "Intercept"
        "Animal bias"
        "\DeltaQ, push better"
        "\DeltaQ, pull better"
        "Rewarded prior"
        "Unrewarded prior"
        "Choice history"
        "Recent switches"];

    rawWeights = nan(numel(featureVars), 3);
    for iFeature = 1:numel(featureVars)
        rawWeights(iFeature,:) = params.(featureVars(iFeature))';
    end
    displayWeights = rawWeights;
    displayWeights(ismember(featureVars, ["beta_pos","beta_neg"]), :) = ...
        8 .* displayWeights(ismember(featureVars, ["beta_pos","beta_neg"]), :);

    figB = figure('Color', 'w', 'Name', 'Figure 5B: fitted weights');
    hold on;
    x = 1:numel(featureVars);
    for state = 1:3
        plot(x, displayWeights(:,state), '-', 'LineWidth', 2, ...
            'Color', colors(state,:), 'DisplayName', sprintf('State %d', state));
    end
    yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xticks(x); xticklabels(featureLabels); xtickangle(45);
    ylabel('Weight');
    title('Fitted RL-GLM-HMM weights');
    legend('Location', 'best', 'Box', 'off');
    styleAxes(gca);

    weightStats = table;
    for state = 1:3
        stateRows = table( ...
            repmat(state, numel(featureVars), 1), featureVars, ...
            rawWeights(:,state), displayWeights(:,state), ...
            'VariableNames', {'state','feature','raw_weight','display_weight'});
        weightStats = [weightStats; stateRows]; %#ok<AGROW>
    end
    weightStats.display_scaling = repmat("none", height(weightStats), 1);
    isDeltaQWeight = ismember(weightStats.feature, ["beta_pos","beta_neg"]);
    weightStats.display_scaling(isDeltaQWeight) = "x8";

    %% Panel C: learning rates
    figC = figure('Color', 'w', 'Name', 'Figure 5C: learning rates');
    bars = bar(states, params.alpha, 0.7, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55);
    bars.CData = colors;
    xticks(states); xticklabels(compose('State %d', states));
    ylabel('Learning rate (\alpha)');
    ylim([0 max(1, 1.15 * max(params.alpha))]);
    title('Learning rate by state');
    styleAxes(gca);
    alphaStats = params(:, {'state','alpha'});

    %% Shared hard-state animal summaries
    animalIDs = unique(T.animalID);
    nAnimals = numel(animalIDs);
    stateMap = T.stateMAP;

    %% Panel D: accuracy by MAP state
    highProbabilityChoice = 3 - (2 .* T.highSide);
    accuracy = double(highProbabilityChoice == T.currChoice);
    accuracyByAnimal = animalStateMeans( ...
        T.animalID, stateMap, accuracy, animalIDs, states) .* 100;
    [accuracyStats, accuracyPHolm] = compareStatesToOne( ...
        accuracyByAnimal, states, [2 1; 3 1], ...
        'medianDelta_percentagePoints');

    figD = plotStateBars(accuracyByAnimal, colors, ...
        'Accuracy (% high-probability choices)', [35 100], true);
    set(figD, 'Name', 'Figure 5D: accuracy');
    yline(gca, 50, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    title(gca, 'Choice accuracy by state');
    annotateStateStars(gca, [2 3], accuracyPHolm, 98);

    %% Panel E: push fraction by MAP state and overall
    isPush = double(T.currChoice == 1);
    pushByAnimal = animalStateMeans( ...
        T.animalID, stateMap, isPush, animalIDs, states);
    overallPush = nan(nAnimals, 1);
    for iAnimal = 1:nAnimals
        overallPush(iAnimal) = mean(isPush(T.animalID == animalIDs(iAnimal)));
    end
    [pushStats, pushPHolm] = compareStateToOverall( ...
        pushByAnimal, overallPush, states);

    figE = plotPushFraction(pushByAnimal, overallPush, colors);
    title(gca, 'Push fraction by state');
    annotateTextByP(gca, 1:3, pushPHolm, 0.965);

    %% Panel F: P(push) across legacy deltaQ
    binCenters = -8:2:8;
    binEdges = -9:2:9;
    minTrialsPerBin = 20;
    nBins = numel(binCenters);
    populationPush = nan(3, nBins);
    effectiveN = nan(3, nBins);
    animalPush = nan(nAnimals, 3, nBins);

    for state = 1:3
        weights = T.(sprintf('stateProb_%d', state));
        for iBin = 1:nBins
            inBin = T.deltaQ >= binEdges(iBin) & ...
                T.deltaQ < binEdges(iBin + 1);
            effectiveN(state,iBin) = sum(weights(inBin), 'omitnan');
            if effectiveN(state,iBin) >= minTrialsPerBin
                populationPush(state,iBin) = ...
                    sum(weights(inBin) .* isPush(inBin), 'omitnan') ./ ...
                    effectiveN(state,iBin);
            end

            for iAnimal = 1:nAnimals
                idx = T.animalID == animalIDs(iAnimal) & ...
                    stateMap == state & inBin;
                if sum(idx) >= minTrialsPerBin
                    animalPush(iAnimal,state,iBin) = mean(isPush(idx));
                end
            end
        end
    end

    figF = figure('Color', 'w', 'Name', 'Figure 5F: legacy deltaQ');
    hold on;
    stateHandles = gobjects(3,1);
    for state = 1:3
        paleColor = 0.82 .* ones(1,3) + 0.18 .* colors(state,:);
        for iAnimal = 1:nAnimals
            yAnimal = squeeze(animalPush(iAnimal,state,:));
            if sum(isfinite(yAnimal)) >= 3
                plot(binCenters, yAnimal, '-', 'Color', paleColor, ...
                    'LineWidth', 0.4, 'HandleVisibility', 'off');
            end
        end
        stateHandles(state) = plot(binCenters, populationPush(state,:), '-o', ...
            'Color', colors(state,:), 'MarkerFaceColor', colors(state,:), ...
            'MarkerEdgeColor', 'none', 'LineWidth', 2, ...
            'DisplayName', sprintf('State %d', state));
    end
    xlabel('Legacy \DeltaQ'); ylabel('P(push)');
    xlim([-8 8]); ylim([0 1]); xticks([-8 -4 0 4 8]);
    legend(stateHandles, 'Location', 'best', 'Box', 'off');
    title('Choice function by state');
    styleAxes(gca);

    deltaQStats = table;
    for state = 1:3
        deltaQStats = [deltaQStats; table( ... %#ok<AGROW>
            repmat(state,nBins,1), binCenters(:), ...
            effectiveN(state,:)', populationPush(state,:)', ...
            'VariableNames', {'state','deltaQ_bin_center', ...
            'effective_posterior_n','posterior_weighted_p_push'})];
    end

    %% Panels G/H: action-specific win-stay; eight-test planned family
    staySummary = makeStaySummary(T, animalIDs, states);
    metricDefs = { ...
        'p_push_winStay',  'Push win-stay',   false; ...
        'p_pull_winStay',  'Pull win-stay',   false; ...
        'p_push_loseStay', 'Push lose-switch', true; ...
        'p_pull_loseStay', 'Pull lose-switch', true};
    comparisons = [2 1; 3 1];
    stayStats = plannedStayTests(staySummary, metricDefs, comparisons);

    pushWin = summaryMetricMatrix(staySummary, ...
        'p_push_winStay', animalIDs, states, false);
    pullWin = summaryMetricMatrix(staySummary, ...
        'p_pull_winStay', animalIDs, states, false);

    figG = plotStateBars(pushWin, colors, 'P(stay | push, win)', [0 1], false);
    set(figG, 'Name', 'Figure 5G: push win-stay');
    title(gca, 'Push win-stay');
    pG = stayPValues(stayStats, "Push win-stay");
    annotateStateStars(gca, [2 3], pG, 0.97);

    figH = plotStateBars(pullWin, colors, 'P(stay | pull, win)', [0 1], false);
    set(figH, 'Name', 'Figure 5H: pull win-stay');
    title(gca, 'Pull win-stay');
    pH = stayPValues(stayStats, "Pull win-stay");
    annotateStateStars(gca, [2 3], pH, 0.97);

    %% Panel I: percentage of trials in each MAP state
    occupancy = nan(nAnimals, 3);
    for iAnimal = 1:nAnimals
        idxAnimal = T.animalID == animalIDs(iAnimal);
        nAnimalTrials = sum(idxAnimal);
        for state = 1:3
            occupancy(iAnimal,state) = ...
                100 .* sum(stateMap(idxAnimal) == state) ./ nAnimalTrials;
        end
    end
    [occupancyStats, occupancyPHolm] = compareStatesToOne( ...
        occupancy, states, comparisons, 'medianDelta_percentagePoints');

    figI = plotStateBars(occupancy, colors, ...
        '% of trials in state', [0 100], true);
    set(figI, 'Name', 'Figure 5I: state occupancy');
    title(gca, 'MAP-state occupancy');
    annotateStateStars(gca, [2 3], occupancyPHolm, 97);

    %% Panels J-O: posterior occupancy by phase and experimental group
    roster = animals(:, {'animal_id','cell_type','condition','sex'});
    roster.Properties.VariableNames{'animal_id'} = 'animalID';
    if height(unique(roster(:, {'animalID'}), 'rows')) ~= height(roster)
        error('Animal roster contains duplicate animal_id rows.');
    end
    nTrialsBeforeRoster = height(T);
    T = innerjoin(T, roster, 'Keys', 'animalID');
    if height(T) ~= nTrialsBeforeRoster
        error('Not every modeled animal had exactly one matching roster row.');
    end
    T.cell_type = lower(string(T.cell_type));
    T.condition = lower(string(T.condition));
    T.sex = lower(string(T.sex));

    T.bin = strings(height(T),1);
    T.bin(ismember(T.assignment, ["PRE1","PRE2","PRE"])) = "PRE";
    T.bin(ismember(T.assignment, ["O1","O2","O3","OPTO1"])) = "OPTO1";
    T.bin(ismember(T.assignment, ["O4","O5","O6","bIntra-Opto"])) = "OPTO2";
    T.bin(ismember(T.assignment, ["POST1","POST2"])) = "POST1";
    T.bin(ismember(T.assignment, ["POST3","POST4"])) = "POST2";

    binOrder = ["PRE","OPTO1","OPTO2","POST1","POST2"];
    conditionOrder = ["ctrl","ib"];
    cellTypes = ["sst","th"];
    posteriorSummary = summarizePosteriorByPhase( ...
        T, binOrder, conditionOrder);

    lmeStats = table;
    figPhase = gobjects(numel(cellTypes),1);
    for iCell = 1:numel(cellTypes)
        cellType = cellTypes(iCell);
        cellSummary = posteriorSummary( ...
            posteriorSummary.cell_type == cellType, :);
        [figPhase(iCell), cellStats] = plotAndFitPhasePanels( ...
            cellSummary, cellType, binOrder, conditionOrder, colors);
        lmeStats = [lmeStats; cellStats]; %#ok<AGROW>
    end

    %% Results and optional export
    results = struct;
    results.transition = transitionTable;
    results.weights = weightStats;
    results.learning_rates = alphaStats;
    results.accuracy = accuracyStats;
    results.push_fraction = pushStats;
    results.deltaQ = deltaQStats;
    results.stay_switch = stayStats;
    results.occupancy = occupancyStats;
    results.posterior_phase_summary = posteriorSummary;
    results.posterior_phase_lme = lmeStats;

    if ~isempty(outputFolder)
        saveFigurePair(figA, outputFolder, 'panelA_transition_matrix');
        saveFigurePair(figB, outputFolder, 'panelB_weights');
        saveFigurePair(figC, outputFolder, 'panelC_learning_rates');
        saveFigurePair(figD, outputFolder, 'panelD_accuracy');
        saveFigurePair(figE, outputFolder, 'panelE_push_fraction');
        saveFigurePair(figF, outputFolder, 'panelF_legacy_deltaQ');
        saveFigurePair(figG, outputFolder, 'panelG_push_win_stay');
        saveFigurePair(figH, outputFolder, 'panelH_pull_win_stay');
        saveFigurePair(figI, outputFolder, 'panelI_state_occupancy');
        saveFigurePair(figPhase(1), outputFolder, 'panelsJ_to_L_SST');
        saveFigurePair(figPhase(2), outputFolder, 'panelsM_to_O_TH');

        writeResultTable(transitionTable, outputFolder, 'panelA_transition_matrix.csv');
        writeResultTable(weightStats, outputFolder, 'panelB_weights.csv');
        writeResultTable(alphaStats, outputFolder, 'panelC_learning_rates.csv');
        writeResultTable(accuracyStats, outputFolder, 'panelD_accuracy_stats.csv');
        writeResultTable(pushStats, outputFolder, 'panelE_push_stats.csv');
        writeResultTable(deltaQStats, outputFolder, 'panelF_deltaQ_summary.csv');
        writeResultTable(stayStats, outputFolder, 'panelsG_H_stay_stats.csv');
        writeResultTable(occupancyStats, outputFolder, 'panelI_occupancy_stats.csv');
        writeResultTable(posteriorSummary, outputFolder, ...
            'panelsJ_to_O_posterior_summary.csv');
        writeResultTable(lmeStats, outputFolder, ...
            'panelsJ_to_O_lme_stats.csv');
        fprintf('Saved Figure 5 source plots and statistics to: %s\n', outputFolder);
    end
end


function requireVariables(T, required, description)
    missing = setdiff(required, T.Properties.VariableNames, 'stable');
    if ~isempty(missing)
        error('%s missing required variables: %s', ...
            description, strjoin(missing, ', '));
    end
end


function assertUniqueTrialKeys(T, keys, description)
    [~, uniqueRows] = unique(T(:, keys), 'rows', 'stable');
    if numel(uniqueRows) ~= height(T)
        error('%s contains duplicate animal-session-trial keys.', description);
    end
end


function transition = readTransitionMatrix(filename)
%READTRANSITIONMATRIX Accept headerless matrices and pandas-style CSVs.
% A pandas DataFrame written with index=False commonly has a first line
% containing the numeric-looking column names "0,1,2". Depending on the
% MATLAB release/import heuristics, readmatrix may treat that line as data.
    transition = readmatrix(filename);

    if isequal(size(transition), [4 3]) && ...
            all(abs(transition(1,:) - [0 1 2]) < 1e-12)
        transition = transition(2:end,:);
    elseif isequal(size(transition), [3 4]) && ...
            all(abs(transition(:,1)' - [0 1 2]) < 1e-12)
        transition = transition(:,2:end);
    end
end


function color = contrastTextColor(value)
    if value < 0.55
        color = [0 0 0];
    else
        color = [1 1 1];
    end
end


function styleAxes(ax)
    set(ax, 'Box', 'off', 'TickDir', 'out', ...
        'FontSize', 11, 'LineWidth', 1.1);
end


function matrix = animalStateMeans(animalID, stateMAP, values, animals, states)
    matrix = nan(numel(animals), numel(states));
    for iAnimal = 1:numel(animals)
        for iState = 1:numel(states)
            idx = animalID == animals(iAnimal) & stateMAP == states(iState);
            if any(idx)
                matrix(iAnimal,iState) = mean(values(idx), 'omitnan');
            end
        end
    end
end


function fig = plotStateBars(values, colors, yLabel, yLimits, connectAnimals)
    fig = figure('Color', 'w');
    hold on;
    means = mean(values, 1, 'omitnan');
    bars = bar(1:3, means, 0.75, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55);
    bars.CData = colors;
    if connectAnimals
        for iAnimal = 1:size(values,1)
            plot(1:3, values(iAnimal,:), '-', ...
                'Color', [0.78 0.78 0.78], 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end
    end
    for state = 1:3
        valid = isfinite(values(:,state));
        scatter(repmat(state,sum(valid),1), values(valid,state), 28, ...
            colors(state,:), 'filled', 'MarkerFaceAlpha', 0.32, ...
            'MarkerEdgeColor', 'none', 'jitter', 'on', ...
            'jitterAmount', 0.08);
    end
    xticks(1:3); xticklabels(compose('State %d', 1:3));
    ylabel(yLabel); ylim(yLimits); xlim([0.5 3.5]);
    styleAxes(gca);
end


function [statsTable, pHolm] = compareStatesToOne( ...
        values, states, comparisons, deltaVariableName)
    nTests = size(comparisons,1);
    comparison = strings(nTests,1);
    medianDelta = nan(nTests,1);
    W = nan(nTests,1);
    nPairs = nan(nTests,1);
    nRanked = nan(nTests,1);
    pRaw = nan(nTests,1);

    for iTest = 1:nTests
        stateA = comparisons(iTest,1);
        stateB = comparisons(iTest,2);
        valueA = values(:, states == stateA);
        valueB = values(:, states == stateB);
        [pRaw(iTest), W(iTest), nPairs(iTest), nRanked(iTest), ...
            medianDelta(iTest)] = pairedSignrank(valueA, valueB);
        comparison(iTest) = sprintf('State %d vs State %d', stateA, stateB);
    end
    pHolm = holmAdjust(pRaw);
    statsTable = table(comparison, medianDelta, W, nPairs, nRanked, ...
        pRaw, pHolm, 'VariableNames', {'comparison', deltaVariableName, ...
        'W','N_pairs','N_ranked','p_raw','p_Holm'});
end


function [statsTable, pHolm] = compareStateToOverall(values, overall, states)
    nTests = numel(states);
    comparison = strings(nTests,1);
    medianDelta = nan(nTests,1);
    W = nan(nTests,1);
    nPairs = nan(nTests,1);
    nRanked = nan(nTests,1);
    pRaw = nan(nTests,1);
    for state = 1:nTests
        [pRaw(state), W(state), nPairs(state), nRanked(state), ...
            medianDelta(state)] = pairedSignrank(values(:,state), overall);
        comparison(state) = sprintf('State %d vs All', states(state));
    end
    pHolm = holmAdjust(pRaw);
    statsTable = table(comparison, medianDelta, W, nPairs, nRanked, ...
        pRaw, pHolm, 'VariableNames', {'comparison','medianDelta', ...
        'W','N_pairs','N_ranked','p_raw','p_Holm'});
end


function [p, W, nPairs, nRanked, medianDelta] = pairedSignrank(a, b)
    keep = isfinite(a) & isfinite(b);
    a = a(keep);
    b = b(keep);
    delta = a - b;
    nPairs = numel(delta);
    nRanked = sum(delta ~= 0);
    medianDelta = median(delta, 'omitnan');
    p = NaN;
    W = NaN;
    if nRanked > 0
        [p, ~, signedRankStats] = signrank(a, b);
        W = signedRankStats.signedrank;
    end
end


function pAdjusted = holmAdjust(pRaw)
    pAdjusted = nan(size(pRaw));
    validIndices = find(isfinite(pRaw));
    if isempty(validIndices)
        return;
    end
    [pSorted, order] = sort(pRaw(validIndices));
    m = numel(pSorted);
    adjustedSorted = (m - (1:m)' + 1) .* pSorted(:);
    adjustedSorted = min(cummax(adjustedSorted), 1);
    pAdjusted(validIndices(order)) = adjustedSorted;
end


function fig = plotPushFraction(pushByAnimal, overallPush, colors)
    fig = figure('Color', 'w', 'Name', 'Figure 5E: push fraction');
    hold on;
    x = [1 2 3 4.5];
    means = [mean(pushByAnimal,1,'omitnan'), mean(overallPush,'omitnan')];
    bars = bar(x, means, 0.75, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55);
    bars.CData = [colors; 0.65 0.65 0.65];
    for state = 1:3
        valid = isfinite(pushByAnimal(:,state));
        scatter(repmat(x(state),sum(valid),1), pushByAnimal(valid,state), ...
            28, colors(state,:), 'filled', 'MarkerFaceAlpha', 0.35, ...
            'MarkerEdgeColor', 'none', 'jitter', 'on', 'jitterAmount', 0.08);
    end
    valid = isfinite(overallPush);
    scatter(repmat(x(4),sum(valid),1), overallPush(valid), 28, ...
        [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha', 0.35, ...
        'MarkerEdgeColor', 'none', 'jitter', 'on', 'jitterAmount', 0.08);
    yline(0.5, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xticks(x); xticklabels(["State 1","State 2","State 3","All"]);
    ylabel('Push fraction'); ylim([0 1]); xlim([0.5 5]);
    styleAxes(gca);
end


function annotateStateStars(ax, xPositions, pValues, yPosition)
    for i = 1:numel(xPositions)
        if isfinite(pValues(i)) && pValues(i) < 0.05
            text(ax, xPositions(i), yPosition, '*', ...
                'HorizontalAlignment', 'center', 'FontSize', 16, ...
                'FontWeight', 'bold');
        end
    end
end


function annotateTextByP(ax, xPositions, pValues, yPosition)
    for i = 1:numel(xPositions)
        if ~isfinite(pValues(i))
            label = '';
        elseif pValues(i) < 0.05
            label = '*';
        else
            label = 'n.s.';
        end
        text(ax, xPositions(i), yPosition, label, ...
            'HorizontalAlignment', 'center', 'FontSize', 12, ...
            'FontWeight', 'bold');
    end
end


function summary = makeStaySummary(T, animalIDs, states)
    T = sortrows(T, {'animalID','iOrig','jOrig'});
    n = height(T);
    sameConsecutiveTrial = [false; ...
        T.animalID(2:end) == T.animalID(1:end-1) & ...
        T.iOrig(2:end) == T.iOrig(1:end-1) & ...
        T.jOrig(2:end) == T.jOrig(1:end-1) + 1];
    previousChoice = [NaN; T.currChoice(1:end-1)];
    previousReward = [NaN; T.currReward(1:end-1)];
    stayed = double(T.currChoice == previousChoice);
    previousChoice(~sameConsecutiveTrial) = NaN;
    previousReward(~sameConsecutiveTrial) = NaN;
    stayed(~sameConsecutiveTrial) = NaN;

    nRows = numel(animalIDs) * numel(states);
    summary = table(zeros(nRows,1), zeros(nRows,1), ...
        nan(nRows,1), nan(nRows,1), nan(nRows,1), nan(nRows,1), ...
        'VariableNames', {'animalID','state','p_push_winStay', ...
        'p_pull_winStay','p_push_loseStay','p_pull_loseStay'});
    row = 0;
    for iAnimal = 1:numel(animalIDs)
        for state = 1:numel(states)
            row = row + 1;
            summary.animalID(row) = animalIDs(iAnimal);
            summary.state(row) = states(state);
            base = T.animalID == animalIDs(iAnimal) & T.stateMAP == states(state);
            summary.p_push_winStay(row) = conditionalMean( ...
                stayed, base & previousChoice == 1 & previousReward == 1);
            summary.p_pull_winStay(row) = conditionalMean( ...
                stayed, base & previousChoice == -1 & previousReward == 1);
            summary.p_push_loseStay(row) = conditionalMean( ...
                stayed, base & previousChoice == 1 & previousReward == 0);
            summary.p_pull_loseStay(row) = conditionalMean( ...
                stayed, base & previousChoice == -1 & previousReward == 0);
        end
    end
end


function value = conditionalMean(values, idx)
    if any(idx)
        value = mean(values(idx), 'omitnan');
    else
        value = NaN;
    end
end


function statsTable = plannedStayTests(summary, metricDefs, comparisons)
    nTests = size(metricDefs,1) * size(comparisons,1);
    metric = strings(nTests,1);
    comparison = strings(nTests,1);
    medianDelta = nan(nTests,1);
    W = nan(nTests,1);
    nPairs = nan(nTests,1);
    nRanked = nan(nTests,1);
    pRaw = nan(nTests,1);
    row = 0;
    for iMetric = 1:size(metricDefs,1)
        variable = metricDefs{iMetric,1};
        label = string(metricDefs{iMetric,2});
        invert = metricDefs{iMetric,3};
        for iComparison = 1:size(comparisons,1)
            row = row + 1;
            stateA = comparisons(iComparison,1);
            stateB = comparisons(iComparison,2);
            A = summary(summary.state == stateA, {'animalID',variable});
            B = summary(summary.state == stateB, {'animalID',variable});
            A.Properties.VariableNames{variable} = 'valueA';
            B.Properties.VariableNames{variable} = 'valueB';
            joined = innerjoin(A, B, 'Keys', 'animalID');
            valueA = joined.valueA;
            valueB = joined.valueB;
            if invert
                valueA = 1 - valueA;
                valueB = 1 - valueB;
            end
            [pRaw(row), W(row), nPairs(row), nRanked(row), ...
                medianDelta(row)] = pairedSignrank(valueA, valueB);
            metric(row) = label;
            comparison(row) = sprintf('State %d vs State %d', stateA, stateB);
        end
    end
    pHolm = holmAdjust(pRaw);
    statsTable = table(metric, comparison, medianDelta, W, nPairs, ...
        nRanked, pRaw, pHolm, 'VariableNames', {'metric','comparison', ...
        'medianDelta','W','N_pairs','N_ranked','p_raw','p_Holm'});
end


function matrix = summaryMetricMatrix(summary, variable, animals, states, invert)
    matrix = nan(numel(animals), numel(states));
    for iAnimal = 1:numel(animals)
        for state = 1:numel(states)
            idx = summary.animalID == animals(iAnimal) & ...
                summary.state == states(state);
            if any(idx)
                matrix(iAnimal,state) = summary.(variable)(idx);
            end
        end
    end
    if invert
        matrix = 1 - matrix;
    end
end


function p = stayPValues(statsTable, metric)
    rows = statsTable.metric == metric;
    subset = statsTable(rows,:);
    desired = ["State 2 vs State 1","State 3 vs State 1"];
    p = nan(2,1);
    for i = 1:2
        idx = subset.comparison == desired(i);
        if any(idx)
            p(i) = subset.p_Holm(idx);
        end
    end
end


function summary = summarizePosteriorByPhase(T, binOrder, conditionOrder)
    valid = ismember(T.bin, binOrder) & ...
        ismember(T.condition, conditionOrder) & ...
        ismember(T.cell_type, ["sst","th"]);
    T = T(valid,:);
    T.bin = categorical(T.bin, binOrder, 'Ordinal', true);
    T.condition = categorical(T.condition, conditionOrder);
    T.animalID = categorical(T.animalID);
    T.sex = categorical(T.sex);
    T.cell_type = categorical(T.cell_type, ["sst","th"]);

    [G, animalID, bin, condition, sex, cell_type] = findgroups( ...
        T.animalID, T.bin, T.condition, T.sex, T.cell_type);
    summary = table(animalID, bin, condition, sex, cell_type);
    for state = 1:3
        inputVariable = sprintf('stateProb_%d', state);
        outputVariable = sprintf('mean_state%d', state);
        summary.(outputVariable) = splitapply( ...
            @(x) mean(x, 'omitnan'), T.(inputVariable), G);
    end
end


function [fig, statsTable] = plotAndFitPhasePanels( ...
        summary, cellType, binOrder, conditionOrder, colors)
    conditionColors = [0 0 0; 1 0 0];
    fig = figure('Color', 'w', 'Name', ...
        sprintf('Figure 5 phase panels: %s', upper(cellType)), ...
        'Position', [100 100 1450 420]);
    tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
    statsTable = table;

    for state = 1:3
        ax = nexttile;
        hold(ax, 'on');
        yVariable = sprintf('mean_state%d', state);
        animalLevels = unique(summary.animalID);
        for iAnimal = 1:numel(animalLevels)
            animalRows = summary(summary.animalID == animalLevels(iAnimal),:);
            if isempty(animalRows)
                continue;
            end
            conditionIndex = find(conditionOrder == string(animalRows.condition(1)), 1);
            paleColor = 0.78 .* ones(1,3) + ...
                0.22 .* conditionColors(conditionIndex,:);
            [~, order] = sort(double(animalRows.bin));
            animalRows = animalRows(order,:);
            plot(ax, double(animalRows.bin), animalRows.(yVariable), '-', ...
                'Color', paleColor, 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end

        for iCondition = 1:numel(conditionOrder)
            means = nan(numel(binOrder),1);
            sems = nan(numel(binOrder),1);
            for iBin = 1:numel(binOrder)
                idx = string(summary.condition) == conditionOrder(iCondition) & ...
                    string(summary.bin) == binOrder(iBin);
                values = summary.(yVariable)(idx);
                means(iBin) = mean(values, 'omitnan');
                sems(iBin) = std(values, [], 'omitnan') ./ ...
                    sqrt(sum(isfinite(values)));
            end
            errorbar(ax, 1:numel(binOrder), means, sems, '-o', ...
                'Color', conditionColors(iCondition,:), ...
                'MarkerFaceColor', conditionColors(iCondition,:), ...
                'MarkerEdgeColor', 'none', 'LineWidth', 1.5, ...
                'MarkerSize', 5, 'DisplayName', conditionOrder(iCondition));
        end

        summary.bin = categorical(string(summary.bin), binOrder, 'Ordinal', true);
        summary.condition = categorical(string(summary.condition), conditionOrder);
        summary.animalID = categorical(summary.animalID);
        formula = sprintf('%s ~ condition * bin + (1 | animalID)', yVariable);
        model = fitlme(summary, formula, 'FitMethod', 'ML');
        anovaTable = anova(model);
        overallP = interactionPValue(anovaTable);
        phaseP = phaseInteractionPValues(model, binOrder(2:end));

        stateStats = table( ...
            repmat(cellType, numel(binOrder), 1), ...
            repmat(state, numel(binOrder), 1), ...
            ["overall"; binOrder(2:end)'], ...
            [overallP; phaseP], ...
            'VariableNames', {'cell_type','state','test','p_value'});
        statsTable = [statsTable; stateStats]; %#ok<AGROW>

        text(ax, 0.98, 0.94, sprintf('p_{interaction} = %.3g', overallP), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right');
        for iPhase = 1:numel(phaseP)
            if isfinite(phaseP(iPhase)) && phaseP(iPhase) < 0.05
                text(ax, iPhase + 1, 0.96, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 16, ...
                    'FontWeight', 'bold');
            end
        end
        xticks(ax, 1:numel(binOrder));
        xticklabels(ax, ["Pre","Opto I","Opto II","Post I","Post II"]);
        xtickangle(ax, 45);
        ylim(ax, [0 1]); xlim(ax, [0.7 5.3]);
        ylabel(ax, sprintf('P(State %d)', state));
        title(ax, sprintf('%s | State %d', upper(cellType), state), ...
            'Color', colors(state,:));
        styleAxes(ax);
        if state == 3
            legend(ax, {'Control','Inhibition'}, 'Location', 'best', 'Box', 'off');
        end
    end
end


function p = interactionPValue(anovaTable)
    termNames = lower(string(anovaTable.Term));
    idx = contains(termNames, 'condition') & contains(termNames, 'bin') & ...
        contains(termNames, ':');
    if ~any(idx)
        error('Could not locate the condition-by-bin interaction in LME ANOVA.');
    end
    p = anovaTable.pValue(find(idx, 1));
end


function p = phaseInteractionPValues(model, phaseNames)
    coefficientNames = lower(string(model.Coefficients.Name));
    p = nan(numel(phaseNames),1);
    for iPhase = 1:numel(phaseNames)
        phaseName = lower(phaseNames(iPhase));
        idx = contains(coefficientNames, 'condition') & ...
            contains(coefficientNames, 'ib') & ...
            contains(coefficientNames, 'bin') & ...
            contains(coefficientNames, phaseName) & ...
            contains(coefficientNames, ':');
        if sum(idx) == 1
            p(iPhase) = model.Coefficients.pValue(idx);
        elseif sum(idx) > 1
            error('Multiple LME coefficients matched phase %s.', phaseNames(iPhase));
        else
            warning('No condition-by-%s LME coefficient was found.', phaseNames(iPhase));
        end
    end
end


function saveFigurePair(fig, outputFolder, baseName)
    exportgraphics(fig, fullfile(outputFolder, [baseName '.png']), ...
        'Resolution', 300);
    exportgraphics(fig, fullfile(outputFolder, [baseName '.pdf']), ...
        'ContentType', 'vector');
end


function writeResultTable(T, outputFolder, filename)
    writetable(T, fullfile(outputFolder, filename));
end
