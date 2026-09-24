function results = analyze_slice_plasticity(sstFile, thFile, outputFolder)
% ANALYZE_SLICE_PLASTICITY Reproduce Figure 4H-I AMPA/NMDA analyses.
%
% Each input CSV must contain animal_id, condition (ctrl or ib), and
% ampa_nmda. Cells are averaged within animal before two-sided Welch tests.
% OUTPUTFOLDER is optional; when supplied, the figure and source/statistics
% tables are saved there.

if nargin < 3
    outputFolder = '';
end

sst = readPlasticityCsv(sstFile, "SST");
th = readPlasticityCsv(thFile, "TH");
cellData = [sst; th];

% Animal is the inferential unit.
[G, cellType, condition, animalID] = findgroups( ...
    cellData.cell_type, cellData.condition, cellData.animal_id);
animalSummary = table(cellType, condition, animalID, ...
    splitapply(@mean, cellData.ampa_nmda, G), ...
    splitapply(@numel, cellData.ampa_nmda, G), ...
    'VariableNames', {'cell_type','condition','animal_id', ...
    'mean_ampa_nmda','n_cells'});

cellTypes = ["SST","TH"];
conditions = ["ctrl","ib"];

nCtrl = zeros(2,1);
nIb = zeros(2,1);
meanCtrl = zeros(2,1);
meanIb = zeros(2,1);
differenceIbMinusCtrl = zeros(2,1);
tCtrlMinusIb = zeros(2,1);
df = zeros(2,1);
pValue = zeros(2,1);

fig = figure('Color','w');
tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

colors.ctrl = [0.45 0.45 0.45];
colors.ib = [0.85 0.25 0.25];

previousRng = rng;
rng(1); % reproducible jitter

% The exploratory script used [1 10], which clipped two TH cells near 14.
yLimits = [1, max(10, ceil(max(cellData.ampa_nmda) + 0.25))];

for i = 1:numel(cellTypes)
    thisType = cellTypes(i);
    C = cellData(cellData.cell_type == thisType,:);
    A = animalSummary(animalSummary.cell_type == thisType,:);

    ctrl = A.mean_ampa_nmda(A.condition == "ctrl");
    ib = A.mean_ampa_nmda(A.condition == "ib");
    if numel(ctrl) < 2 || numel(ib) < 2
        error('%s requires at least two animals per condition.', thisType);
    end

    [~, pValue(i), ~, stats] = ttest2(ctrl, ib, 'Vartype','unequal');
    nCtrl(i) = numel(ctrl);
    nIb(i) = numel(ib);
    meanCtrl(i) = mean(ctrl);
    meanIb(i) = mean(ib);
    differenceIbMinusCtrl(i) = mean(ib) - mean(ctrl);
    tCtrlMinusIb(i) = stats.tstat;
    df(i) = stats.df;

    ax = nexttile;
    hold(ax, 'on');

    for c = 1:numel(conditions)
        thisCondition = conditions(c);
        cells = C(C.condition == thisCondition,:);
        animals = A(A.condition == thisCondition,:);
        color = colors.(char(thisCondition));

        scatter(ax, c + 0.10*randn(height(cells),1), cells.ampa_nmda, ...
            28, 'o', 'MarkerFaceColor',color, 'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.22);
        scatter(ax, c + 0.035*randn(height(animals),1), ...
            animals.mean_ampa_nmda, 60, 'o', ...
            'MarkerFaceColor',color, 'MarkerEdgeColor','none');

        m = mean(animals.mean_ampa_nmda);
        sem = std(animals.mean_ampa_nmda) / sqrt(height(animals));
        errorbar(ax, c, m, sem, 'k', 'LineStyle','none', ...
            'LineWidth',0.75, 'CapSize',12);
        plot(ax, c, m, 'k_', 'MarkerSize',18, 'LineWidth',1.5);
    end

    xlim(ax, [0.5 2.5]);
    ylim(ax, yLimits);
    xticks(ax, [1 2]);
    xticklabels(ax, {'Ctrl','IB'});
    ylabel(ax, 'AMPA/NMDA ratio');
    title(ax, thisType);
    box(ax, 'off');
    set(ax, 'FontName','Arial', 'FontSize',12, 'LineWidth',0.75);
end

rng(previousRng);

welchTests = table(cellTypes', nCtrl, nIb, meanCtrl, meanIb, ...
    differenceIbMinusCtrl, tCtrlMinusIb, df, pValue, ...
    'VariableNames', {'cell_type','n_ctrl_animals','n_ib_animals', ...
    'mean_ctrl','mean_ib','difference_ib_minus_ctrl', ...
    't_statistic_ctrl_minus_ib','df','p_value'});

disp(welchTests)

results = struct( ...
    'figure', fig, ...
    'cell_data', cellData, ...
    'animal_summary', animalSummary, ...
    'welch_tests', welchTests);

if ~isempty(outputFolder)
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    exportgraphics(fig, fullfile(outputFolder, ...
        'figure4_panelsH_to_I.png'), 'Resolution',300);
    exportgraphics(fig, fullfile(outputFolder, ...
        'figure4_panelsH_to_I.pdf'), 'ContentType','vector');
    writetable(cellData, fullfile(outputFolder, ...
        'figure4_panelsH_to_I_cell_data.csv'));
    writetable(animalSummary, fullfile(outputFolder, ...
        'figure4_panelsH_to_I_animal_summary.csv'));
    writetable(welchTests, fullfile(outputFolder, ...
        'figure4_panelsH_to_I_welch_tests.csv'));
end
end

function T = readPlasticityCsv(filename, cellType)
if ~isfile(filename)
    error('%s input file does not exist: %s', cellType, filename);
end

T = readtable(filename, 'TextType','string');
required = {'animal_id','condition','ampa_nmda'};
missing = setdiff(required, T.Properties.VariableNames);
if ~isempty(missing)
    error('%s input is missing columns: %s', ...
        cellType, strjoin(missing, ', '));
end
T = T(:, required);

T.animal_id = strtrim(string(T.animal_id));
T.condition = lower(strtrim(string(T.condition)));
if any(~ismember(T.condition, ["ctrl","ib"]))
    error('%s condition values must be ctrl or ib.', cellType);
end
if ~isnumeric(T.ampa_nmda) || any(~isfinite(T.ampa_nmda))
    error('%s ampa_nmda values must be finite numbers.', cellType);
end

T.cell_type = repmat(cellType, height(T), 1);
T = movevars(T, 'cell_type', 'Before', 1);
end
