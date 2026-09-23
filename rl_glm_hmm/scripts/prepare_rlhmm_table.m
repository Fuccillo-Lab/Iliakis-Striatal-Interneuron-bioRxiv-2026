function T_export = prepare_rlhmm_table(dataFolder, calendarFile, outputFile, excludeBoxSwitch)
%PREPARE_RLHMM_TABLE Build the trial table used by the paper's RL-GLM-HMM.
%
% T_EXPORT = PREPARE_RLHMM_TABLE(DATAFOLDER, CALENDARFILE, OUTPUTFILE)
% reads all trial CSV files in DATAFOLDER, joins session metadata from
% CALENDARFILE, constructs the historical predictors and complete-case
% trial set, and writes OUTPUTFILE.
%
% T_EXPORT = PREPARE_RLHMM_TABLE(..., EXCLUDEBOXSWITCH) controls the
% historical box-stability exclusion (default: true). When true, an animal
% is excluded if it appears in more than one box anywhere in the imported
% data. This is the rule used for the paper's model input.
%
% The exported table intentionally omits deltaQ_z and deltaQ_all_z. Neither
% was used by the final RL-GLM-HMM. Accordingly, this function does not call
% the legacy globalQF.m script or require its optimization dependencies.
%
% Historical compatibility notes
% ------------------------------
% 1. The complete-case filter retains the original Xfinal variables,
%    including motor variables and deltaQ, so it reproduces the historical
%    trial-selection rule even though those columns are not exported.
% 2. animalBias is calculated before the later latency/history/complete-case
%    exclusions, as in the original analysis.
% 3. recentChoiceFrac_z is divided by its within-animal SD without mean
%    centering. switchFrac_z is conventionally centered and scaled.
%
% Example
% -------
% T = prepare_rlhmm_table( ...
%     "C:\Data\March2026\CsvOutput", ...
%     "C:\Data\March2026\calendar.csv", ...
%     "C:\Data\March2026\modelTfinal_for_rlhmm.csv", ...
%     true);

    if nargin < 4 || isempty(excludeBoxSwitch)
        excludeBoxSwitch = true;
    end

    dataFolder = char(dataFolder);
    calendarFile = char(calendarFile);
    outputFile = char(outputFile);

    if ~isfolder(dataFolder)
        error('Trial-data folder does not exist: %s', dataFolder);
    end
    if ~isfile(calendarFile)
        error('Calendar file does not exist: %s', calendarFile);
    end

    %% Read calendar and trial CSVs
    calendarOpts = detectImportOptions(calendarFile);
    if ismember('date', calendarOpts.VariableNames)
        calendarOpts = setvartype(calendarOpts, 'date', 'string');
    end
    calendar = readtable(calendarFile, calendarOpts);
    requireVariables(calendar, {'animal','date','assignment','box'}, 'calendar');

    calendar.dateKey = datetime(string(calendar.date), ...
        'InputFormat', 'dd-MM-yy', 'Format', 'yyyy-MM-dd');
    calendar.animalID = calendar.animal;

    files = dir(fullfile(dataFolder, '*.csv'));
    if isempty(files)
        error('No CSV files found in: %s', dataFolder);
    end

    fprintf('Reading %d trial CSV files...\n', numel(files));
    tic
    dataCell = cell(numel(files), 1);

    for iFile = 1:numel(files)
        filePath = fullfile(dataFolder, files(iFile).name);
        Ti = readtable(filePath);
        Ti.source_file = repmat(string(files(iFile).name), height(Ti), 1);
        dataCell{iFile} = Ti;

        reportEvery = max(1, floor(numel(files) / 20));
        if mod(iFile, reportEvery) == 0 || iFile == numel(files)
            fprintf('  Read %d/%d files\n', iFile, numel(files));
        end
    end

    allData = vertcat(dataCell{:});
    clear dataCell Ti
    fprintf('CSV load done: %d trials, %.2f sec\n', height(allData), toc);

    requiredTrialVariables = { ...
        'animalID','date','iOrig','jOrig','pushPull','trialOutcome', ...
        'choiceLatency','peakDispl','meanVel','peakVel','dirConsistency', ...
        'pathLength','numBouts','optoTrial','rewardDirection','phase', ...
        'optoType','deltaQ'};
    requireVariables(allData, requiredTrialVariables, 'trial CSVs');

    if isdatetime(allData.date)
        allData.dateKey = allData.date;
        allData.dateKey.Format = 'yyyy-MM-dd';
    else
        allData.dateKey = datetime(string(allData.date), ...
            'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
    end

    %% Join calendar metadata
    sessionLookup = calendar(:, {'animalID','dateKey','assignment','box'});
    [Gcalendar, keyAnimal, keyDate] = findgroups( ...
        sessionLookup.animalID, sessionLookup.dateKey);
    nPerSession = splitapply(@numel, sessionLookup.assignment, Gcalendar);

    if any(nPerSession > 1)
        duplicateSessions = table( ...
            keyAnimal(nPerSession > 1), ...
            keyDate(nPerSession > 1), ...
            nPerSession(nPerSession > 1), ...
            'VariableNames', {'animalID','date','nRows'});
        disp(duplicateSessions);
        error('Calendar contains duplicate animalID-by-date rows.');
    end

    allData = outerjoin(allData, sessionLookup, ...
        'Keys', {'animalID','dateKey'}, ...
        'MergeKeys', true, ...
        'Type', 'left');

    fprintf('Trials missing assignment: %d\n', sum(ismissing(allData.assignment)));
    fprintf('Trials missing box: %d\n', sum(ismissing(allData.box)));

    allData.dateKey = [];
    allData.assignment = string(allData.assignment);
    T = allData;

    %% Historical box-stability exclusion
    if excludeBoxSwitch
        [Ganimal, animalIDs] = findgroups(T.animalID);
        nUniqueBoxes = splitapply( ...
            @(x) numel(unique(x(~isnan(x)))), T.box, Ganimal);
        excludedAnimals = animalIDs(nUniqueBoxes > 1);

        if ~isempty(excludedAnimals)
            fprintf('Excluding %d animals that appeared in >1 box:\n', ...
                numel(excludedAnimals));
            disp(table(excludedAnimals, 'VariableNames', {'animalID'}));
            T = T(~ismember(T.animalID, excludedAnimals), :);
        else
            fprintf('No animals excluded for box switching.\n');
        end
    end

    %% Animal-level choice bias (computed before later exclusions)
    biasEligible = ismember(T.trialOutcome, [1 2]) & ...
                   ismember(T.pushPull, [1 2]);
    biasData = T(biasEligible, :);
    biasAnimals = unique(biasData.animalID);
    animalBiasValues = nan(numel(biasAnimals), 1);

    for iAnimal = 1:numel(biasAnimals)
        idx = biasData.animalID == biasAnimals(iAnimal);
        choicesForBias = biasData.pushPull(idx);
        animalBiasValues(iAnimal) = ...
            mean(choicesForBias == 1) - mean(choicesForBias == 2);
    end

    T.animalBias = nan(height(T), 1);
    [hasBias, biasLocation] = ismember(T.animalID, biasAnimals);
    T.animalBias(hasBias) = animalBiasValues(biasLocation(hasBias));

    %% Trial coding and lagged variables
    T = sortrows(T, {'animalID','iOrig','jOrig'});
    n = height(T);

    choice = nan(n, 1);
    choice(T.pushPull == 1) = 1;   % push
    choice(T.pushPull == 2) = -1;  % pull
    reward = double(T.trialOutcome == 1);

    valid = ismember(T.pushPull, [1 2]) & ...
            T.choiceLatency < 11000 & T.choiceLatency >= 0;

    sameAnimalSession = [false; ...
        T.animalID(2:end) == T.animalID(1:end-1) & ...
        T.iOrig(2:end) == T.iOrig(1:end-1)];

    sameConsecutiveTrial = [false; ...
        T.animalID(2:end) == T.animalID(1:end-1) & ...
        T.iOrig(2:end) == T.iOrig(1:end-1) & ...
        T.jOrig(2:end) == T.jOrig(1:end-1) + 1];

    prevChoiceRaw = [nan; choice(1:end-1)];
    prevRewardRaw = [nan; reward(1:end-1)];

    prevChoice = nan(n, 1);
    prevReward = nan(n, 1);
    currChoice = nan(n, 1);
    currReward = nan(n, 1);

    prevChoice(sameAnimalSession) = prevChoiceRaw(sameAnimalSession);
    prevReward(sameAnimalSession) = prevRewardRaw(sameAnimalSession);
    currChoice(sameAnimalSession) = choice(sameAnimalSession);
    currReward(sameAnimalSession) = reward(sameAnimalSession);
    prevUnreward = abs(1 - prevReward);

    %% Historical complete-case variables
    % These variables are retained here solely to reproduce the original
    % rmmissing/modelTfinal row set. They are not all exported or fitted.
    animalBias = T.animalBias;
    latency = T.choiceLatency;
    displ = abs(T.peakDispl);
    meanVel = abs(T.meanVel);
    peakVel = abs(T.peakVel);
    dirConsistency = abs(T.dirConsistency);
    pathLength = abs(T.pathLength);
    numBouts = T.numBouts;

    animalBias(~valid) = NaN;
    latency(~valid) = NaN;
    displ(~valid) = NaN;
    meanVel(~valid) = NaN;
    peakVel(~valid) = NaN;
    dirConsistency(~valid) = NaN;
    pathLength(~valid) = NaN;
    numBouts(~valid) = NaN;

    %% Optogenetic history variables
    T.previousOptoTrial = zeros(n, 1);
    idxPrevOpto = find(sameConsecutiveTrial);
    T.previousOptoTrial(idxPrevOpto) = T.optoTrial(idxPrevOpto - 1);

    lambda = 0.9;
    T.trialOptoLeak = zeros(n, 1);
    [Gsession, ~, ~] = findgroups(T.animalID, T.iOrig);
    nSessions = max(Gsession);

    fprintf('Computing trialOptoLeak...\n');
    tic
    for iSession = 1:nSessions
        idx = find(Gsession == iSession);
        trialNumber = T.jOrig(idx);
        opto = T.optoTrial(idx);
        leak = zeros(numel(idx), 1);

        for ii = 2:numel(idx)
            nMissing = trialNumber(ii) - trialNumber(ii-1) - 1;
            if nMissing < 0
                error('Trials are not sorted within session.');
            end

            leak(ii) = lambda^(nMissing + 1) * leak(ii-1);
            if trialNumber(ii) == trialNumber(ii-1) + 1
                leak(ii) = leak(ii) + opto(ii-1);
            end
        end

        T.trialOptoLeak(idx) = leak;
    end
    fprintf('trialOptoLeak done: %.2f sec\n', toc);

    gamma = 0.6;
    [sessionIndex, sessionAnimal, sessionIOrig] = ...
        findgroups(T.animalID, T.iOrig);
    sessionOptoDensity = splitapply(@mean, T.optoTrial, sessionIndex);
    T.sessionOptoDensity = sessionOptoDensity(sessionIndex);

    sessionOptoLeak = zeros(nSessions, 1);
    for iSession = 2:nSessions
        if sessionAnimal(iSession) == sessionAnimal(iSession-1)
            sessionGap = sessionIOrig(iSession) - sessionIOrig(iSession-1);
            if sessionGap < 1
                error('Sessions are not sorted by animalID and iOrig.');
            end
            sessionOptoLeak(iSession) = ...
                gamma^sessionGap * sessionOptoLeak(iSession-1) + ...
                sessionOptoDensity(iSession-1);
        end
    end
    T.sessionOptoLeak = sessionOptoLeak(sessionIndex);

    T.trialOptoLeak_z = zscoreWithinGroupFast( ...
        T.trialOptoLeak, T.animalID, true);
    T.sessionOptoLeak_z = zscoreWithinGroupFast( ...
        T.sessionOptoLeak, T.animalID, true);

    %% Choice-history variables
    recentWindow = 5;
    switchWindow = 12;
    recentChoiceFrac = laggedRecentChoiceFracFast( ...
        choice, T.animalID, T.iOrig, recentWindow);
    switchFrac = laggedSwitchFracFast( ...
        choice, T.animalID, T.iOrig, switchWindow);

    %% Reproduce historical complete-case trial selection
    animalID = T.animalID;
    iOrig = T.iOrig;
    jOrig = T.jOrig;
    assignment = T.assignment;
    isOptoTrial = T.optoTrial;
    phase = T.phase;
    isOptoSession = T.optoType;
    highSide = T.rewardDirection;
    deltaQ = T.deltaQ;
    boxNo = T.box;
    previousOptoTrial = T.previousOptoTrial;
    trialOptoLeak = T.trialOptoLeak;
    sessionOptoLeak = T.sessionOptoLeak;
    trialOptoLeak_z = T.trialOptoLeak_z;
    sessionOptoLeak_z = T.sessionOptoLeak_z;

    Xfinal = table( ...
        animalID, iOrig, jOrig, assignment, isOptoTrial, phase, ...
        isOptoSession, prevChoice, prevReward, prevUnreward, deltaQ, ...
        highSide, animalBias, recentChoiceFrac, boxNo, switchFrac, ...
        currChoice, currReward, latency, displ, meanVel, peakVel, ...
        dirConsistency, pathLength, numBouts, previousOptoTrial, ...
        trialOptoLeak, sessionOptoLeak, trialOptoLeak_z, ...
        sessionOptoLeak_z);

    Y = choice;
    Xfinal = Xfinal(valid, :);
    Y = Y(valid);

    modelTfinal = Xfinal;
    modelTfinal.Y = (Y == 1);

    nBeforeCompleteCase = height(modelTfinal);
    modelTfinal = rmmissing(modelTfinal);
    fprintf('Complete-case filter removed %d/%d otherwise-valid rows.\n', ...
        nBeforeCompleteCase - height(modelTfinal), nBeforeCompleteCase);

    %% Scale fitted predictors after complete-case selection
    modelTfinal = sortrows(modelTfinal, {'animalID','iOrig','jOrig'});
    modelTfinal.switchFrac_z = zscoreWithinGroupFast( ...
        modelTfinal.switchFrac, modelTfinal.animalID, true);
    modelTfinal.recentChoiceFrac_z = zscoreWithinGroupFast( ...
        modelTfinal.recentChoiceFrac, modelTfinal.animalID, false);

    %% Build fitted and auxiliary interaction terms
    modelTfinal.prevChoice_switchFrac_z = ...
        modelTfinal.prevChoice .* modelTfinal.switchFrac_z;
    modelTfinal.prevChoice_prevReward = ...
        modelTfinal.prevChoice .* modelTfinal.prevReward;
    modelTfinal.prevChoice_prevUnreward = ...
        modelTfinal.prevChoice .* modelTfinal.prevUnreward;
    modelTfinal.isOptoSession_previousOptoTrial = ...
        modelTfinal.isOptoSession .* modelTfinal.previousOptoTrial;
    modelTfinal.isOptoSession_trialOptoLeak_z = ...
        modelTfinal.isOptoSession .* modelTfinal.trialOptoLeak_z;
    modelTfinal.prevChoice_sessionOptoLeak_z = ...
        modelTfinal.prevChoice .* modelTfinal.sessionOptoLeak_z;
    modelTfinal.isOptoSession_prevChoice_previousOptoTrial = ...
        modelTfinal.isOptoSession .* modelTfinal.prevChoice .* ...
        modelTfinal.previousOptoTrial;
    modelTfinal.isOptoSession_prevChoice_trialOptoLeak_z = ...
        modelTfinal.isOptoSession .* modelTfinal.prevChoice .* ...
        modelTfinal.trialOptoLeak_z;
    modelTfinal.intercept = ones(height(modelTfinal), 1);

    %% Export (deltaQ_z/deltaQ_all_z deliberately omitted)
    exportVariables = { ...
        'animalID', ...
        'iOrig', ...
        'jOrig', ...
        'currChoice', ...
        'currReward', ...
        'intercept', ...
        'prevChoice', ...
        'prevReward', ...
        'animalBias', ...
        'boxNo', ...
        'recentChoiceFrac_z', ...
        'isOptoSession', ...
        'previousOptoTrial', ...
        'trialOptoLeak', ...
        'sessionOptoLeak', ...
        'trialOptoLeak_z', ...
        'sessionOptoLeak_z', ...
        'switchFrac_z', ...
        'prevChoice_switchFrac_z', ...
        'prevChoice_prevReward', ...
        'prevChoice_prevUnreward', ...
        'isOptoSession_previousOptoTrial', ...
        'isOptoSession_trialOptoLeak_z', ...
        'prevChoice_sessionOptoLeak_z', ...
        'isOptoSession_prevChoice_previousOptoTrial', ...
        'isOptoSession_prevChoice_trialOptoLeak_z', ...
        'Y'};

    T_export = modelTfinal(:, exportVariables);

    %% Sanity checks and write
    requiredForFitter = { ...
        'animalID','iOrig','jOrig','currChoice','currReward', ...
        'prevChoice_prevReward','prevChoice_prevUnreward', ...
        'animalBias','recentChoiceFrac_z','prevChoice_switchFrac_z'};
    requireVariables(T_export, requiredForFitter, 'RL-GLM-HMM export');

    if any(any(ismissing(T_export)))
        error('Export unexpectedly contains missing values.');
    end

    [~, uniqueKeyRows] = unique( ...
        T_export(:, {'animalID','iOrig','jOrig'}), 'rows', 'stable');
    if numel(uniqueKeyRows) ~= height(T_export)
        error('Export contains duplicate animalID-iOrig-jOrig trial keys.');
    end

    outputFolder = fileparts(outputFile);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    writetable(T_export, outputFile);

    nAnimals = numel(unique(T_export.animalID));
    nSessionsExported = height(unique( ...
        T_export(:, {'animalID','iOrig'}), 'rows'));
    fprintf('\nWrote RL-GLM-HMM table: %s\n', outputFile);
    fprintf('  %d trials | %d sessions | %d animals\n', ...
        height(T_export), nSessionsExported, nAnimals);
end


function requireVariables(T, required, description)
%REQUIREVARIABLES Error clearly when an input table lacks required columns.
    missing = setdiff(required, T.Properties.VariableNames, 'stable');
    if ~isempty(missing)
        error('%s missing required variables: %s', ...
            description, strjoin(missing, ', '));
    end
end


function out = zscoreWithinGroupFast(x, groupID, centerData)
%ZSCOREWITHINGROUPFAST Scale a vector separately within each group.
% centerData=true:  (x - mean) / SD
% centerData=false: x / SD
    if nargin < 3
        centerData = true;
    end

    out = nan(size(x));
    [G, groupKeys] = findgroups(groupID);

    for iGroup = 1:numel(groupKeys)
        idx = (G == iGroup);
        xx = x(idx);
        mu = mean(xx, 'omitnan');
        sigma = std(xx, 'omitnan');

        if isnan(sigma) || sigma == 0
            out(idx) = zeros(sum(idx), 1);
        elseif centerData
            out(idx) = (xx - mu) ./ sigma;
        else
            out(idx) = xx ./ sigma;
        end
    end
end


function recentChoiceFrac = laggedRecentChoiceFracFast( ...
        choice, animalID, iOrig, window)
%LAGGEDRECENTCHOICEFRACFAST Prior-window push-minus-pull fraction.
% The current trial is excluded. NaN choices within the preceding WINDOW
% trial rows are omitted from the denominator.
    n = numel(choice);
    recentChoiceFrac = nan(n, 1);
    [G, ~, ~] = findgroups(animalID, iOrig);
    nGroups = max(G);

    for iGroup = 1:nGroups
        idx = find(G == iGroup);
        c = choice(idx);

        isPush = double(c == 1);
        isPull = double(c == -1);
        isValidChoice = double(~isnan(c));

        pushLag = [0; isPush(1:end-1)];
        pullLag = [0; isPull(1:end-1)];
        validLag = [0; isValidChoice(1:end-1)];

        pushSum = movsum(pushLag, [window-1 0], 'Endpoints', 'shrink');
        pullSum = movsum(pullLag, [window-1 0], 'Endpoints', 'shrink');
        validSum = movsum(validLag, [window-1 0], 'Endpoints', 'shrink');

        out = nan(size(c));
        hasHistory = validSum > 0;
        out(hasHistory) = ...
            (pushSum(hasHistory) - pullSum(hasHistory)) ./ ...
            validSum(hasHistory);
        recentChoiceFrac(idx) = out;
    end
end


function switchFrac = laggedSwitchFracFast(choice, animalID, iOrig, window)
%LAGGEDSWITCHFRACFAST Switch rate in the preceding WINDOW trial rows.
% The current trial is excluded. NaN choices are removed before switches
% are counted.
    n = numel(choice);
    switchFrac = nan(n, 1);
    [G, ~, ~] = findgroups(animalID, iOrig);
    nGroups = max(G);

    for iGroup = 1:nGroups
        idx = find(G == iGroup);
        c = choice(idx);
        out = nan(size(c));

        for ii = 1:numel(c)
            firstHistoryRow = max(1, ii - window);
            historyChoices = c(firstHistoryRow:ii-1);
            historyChoices = historyChoices(~isnan(historyChoices));

            if numel(historyChoices) >= 2
                out(ii) = sum(diff(historyChoices) ~= 0) / ...
                    (numel(historyChoices) - 1);
            end
        end

        switchFrac(idx) = out;
    end
end
