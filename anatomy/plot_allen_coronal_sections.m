function [fig,sites,animalSummary] = plot_allen_coronal_sections( ...
    T,templateFile,varargin)
%PLOT_ALLEN_CORONAL_SECTIONS Plot targeting coordinates on Allen CCFv3.
%
% [fig,sites,animalSummary] = plot_allen_coronal_sections( ...
%     T,templateFile,Name,Value,...)
%
% T is a prefiltered animal table with bilateral histology coordinates.
% templateFile is a local copy of average_template_10.nrrd. Coordinates are
% assigned to coronal AP panels, converted into Allen voxels, and overlaid
% on the corresponding average-template slices.
%
% Required columns (default schema)
% ---------------------------------
% animalID
% good_targeting
% left_histo_fiber_ap,  left_histo_fiber_dv,  left_histo_fiber_ml
% right_histo_fiber_ap, right_histo_fiber_dv, right_histo_fiber_ml
%
% Principal name-value options
% ----------------------------
% LeftInclusionVariable  : optional hemisphere-specific inclusion column
% RightInclusionVariable : optional hemisphere-specific inclusion column
% ShowAnimals            : "all", "included", or "excluded"
% PanelSpacingMm         : AP interval, anchored at +0.75 mm (default 0.30)
% PanelAssignment        : "pair-mean" (default) or "site"
% ConnectPairs           : connect bilateral sites assigned to one panel
% CropMm                 : [] or [MLmin MLmax DVmin DVmax] in mm
% AnnotationFile         : optional Allen annotation NRRD for boundaries
% SaveFigure             : export PDF and 600-dpi PNG (default true)
% OutputStem             : output path without an extension
% FigureTitle            : optional title above the panel grid

assert(istable(T),'The first input must be a MATLAB table.');
assert(~isempty(T),'The supplied table contains no rows.');

templateFile = string(templateFile);
assert(isscalar(templateFile) && isfile(templateFile), ...
    'Allen template not found: %s',templateFile);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'AnnotationFile',"",@isTextScalar);
addParameter(p,'AnimalVariable',"animalID",@isTextScalar);
addParameter(p,'InclusionVariable',"good_targeting",@isTextScalar);
addParameter(p,'LeftInclusionVariable',"",@isTextScalar);
addParameter(p,'RightInclusionVariable',"",@isTextScalar);
addParameter(p,'IncludedValues',["y","yes","true","1","include","included"], ...
    @(x)isstring(x) || ischar(x) || iscellstr(x));
addParameter(p,'LeftPrefix',"left_histo_fiber_",@isTextScalar);
addParameter(p,'RightPrefix',"right_histo_fiber_",@isTextScalar);
addParameter(p,'ShowAnimals',"all", ...
    @(x)any(strcmpi(string(x),["all","included","excluded"])));
addParameter(p,'PanelAnchorMm',0.75,@isNumericScalar);
addParameter(p,'PanelSpacingMm',0.30, ...
    @(x)isNumericScalar(x) && x>0);
addParameter(p,'PanelAssignment',"pair-mean", ...
    @(x)any(strcmpi(string(x),["pair-mean","site"])));
addParameter(p,'ConnectPairs',true,@isLogicalScalar);
addParameter(p,'CropMm',[], ...
    @(x)isnumeric(x) && (isempty(x) || numel(x)==4));
addParameter(p,'DisplayRange',[0 400], ...
    @(x)isnumeric(x) && numel(x)==2 && x(2)>x(1));
addParameter(p,'IncludedColor',[59 84 162]./255,@isRgbTriplet);
addParameter(p,'ExcludedColor',[0.75 0.75 0.75],@isRgbTriplet);
addParameter(p,'ShowLegend',false,@isLogicalScalar);
addParameter(p,'SaveFigure',true,@isLogicalScalar);
addParameter(p,'OutputStem',"",@isTextScalar);
addParameter(p,'FigureTitle',"",@isTextScalar);
parse(p,varargin{:});
opt = p.Results;

annotationFile = string(opt.AnnotationFile);
showAnimals = lower(string(opt.ShowAnimals));
panelAnchor = double(opt.PanelAnchorMm);
panelSpacing = double(opt.PanelSpacingMm);

%% Reshape bilateral coordinates into one row per site

variableNames = string(T.Properties.VariableNames);
animalVariable = string(opt.AnimalVariable);
inclusionVariable = string(opt.InclusionVariable);
leftInclusionVariable = string(opt.LeftInclusionVariable);
rightInclusionVariable = string(opt.RightInclusionVariable);
leftColumns = string(opt.LeftPrefix)+["ap","dv","ml"];
rightColumns = string(opt.RightPrefix)+["ap","dv","ml"];

requiredVariables = [animalVariable,leftColumns,rightColumns];
usesSharedInclusion = ...
    (strlength(leftInclusionVariable)==0 || ...
     strlength(rightInclusionVariable)==0) && ...
    strlength(inclusionVariable)>0;
if usesSharedInclusion
    requiredVariables(end+1) = inclusionVariable;
end
if strlength(leftInclusionVariable)>0
    requiredVariables(end+1) = leftInclusionVariable;
end
if strlength(rightInclusionVariable)>0
    requiredVariables(end+1) = rightInclusionVariable;
end

missingVariables = setdiff(requiredVariables,variableNames);
assert(isempty(missingVariables), ...
    'Input table is missing required variable(s): %s', ...
    strjoin(missingVariables,', '));

if usesSharedInclusion
    sharedIncluded = normalizeInclusion( ...
        T.(inclusionVariable),opt.IncludedValues);
else
    sharedIncluded = true(height(T),1);
end

leftIncluded = sharedIncluded;
rightIncluded = sharedIncluded;
if strlength(leftInclusionVariable)>0
    leftIncluded = normalizeInclusion( ...
        T.(leftInclusionVariable),opt.IncludedValues);
end
if strlength(rightInclusionVariable)>0
    rightIncluded = normalizeInclusion( ...
        T.(rightInclusionVariable),opt.IncludedValues);
end

nAnimals = height(T);
sourceRow = (1:nAnimals)';
leftSites = table( ...
    T.(animalVariable),repmat("left",nAnimals,1), ...
    numericColumn(T.(leftColumns(1))), ...
    numericColumn(T.(leftColumns(2))), ...
    numericColumn(T.(leftColumns(3))), ...
    leftIncluded,sourceRow, ...
    'VariableNames',{'animalID','hemisphere','AP','DV','ML', ...
    'included','sourceRow'});
rightSites = table( ...
    T.(animalVariable),repmat("right",nAnimals,1), ...
    numericColumn(T.(rightColumns(1))), ...
    numericColumn(T.(rightColumns(2))), ...
    numericColumn(T.(rightColumns(3))), ...
    rightIncluded,sourceRow, ...
    'VariableNames',{'animalID','hemisphere','AP','DV','ML', ...
    'included','sourceRow'});

sites = [leftSites;rightSites];
validCoordinate = isfinite(sites.AP) & isfinite(sites.DV) & isfinite(sites.ML);
if any(~validCoordinate)
    warning('Ignoring %d site row(s) with missing AP/DV/ML.', ...
        nnz(~validCoordinate));
end
sites = sites(validCoordinate,:);

switch showAnimals
    case "included"
        sites = sites(sites.included,:);
    case "excluded"
        sites = sites(~sites.included,:);
end
assert(~isempty(sites),'No sites remain after ShowAnimals=%s.',showAnimals);

sites.pairLabel = string(sites.animalID);
assert(~any(ismissing(sites.pairLabel)), ...
    'Every plotted site must have a non-missing animal identifier.');

%% Assign sites to AP panels

[animalGroup,animalIDs] = findgroups(sites.pairLabel);
animalMeanAP = splitapply(@mean,sites.AP,animalGroup);
animalNSites = splitapply(@numel,sites.AP,animalGroup);
animalIncluded = splitapply(@(x)all(x),sites.included,animalGroup);
animalPanelTick = round((animalMeanAP-panelAnchor)./panelSpacing);

switch lower(string(opt.PanelAssignment))
    case "pair-mean"
        sites.panelTick = animalPanelTick(animalGroup);
    case "site"
        sites.panelTick = round((sites.AP-panelAnchor)./panelSpacing);
end

sites.panelAP = panelAnchor + sites.panelTick.*panelSpacing;
sites.meanAnimalAP = animalMeanAP(animalGroup);
animalSummary = table( ...
    animalIDs,animalMeanAP, ...
    panelAnchor + animalPanelTick.*panelSpacing, ...
    animalNSites,animalIncluded, ...
    'VariableNames',{'animalID','meanAP','panelAP','nSites','included'});

% Always retain the +0.75-mm anchor panel for consistent figure layouts.
minimumTick = min([sites.panelTick;0]);
maximumTick = max([sites.panelTick;0]);
panelTicks = maximumTick:-1:minimumTick;
panelList = panelAnchor + panelTicks.*panelSpacing;

%% Load and orient Allen volumes

fprintf('Loading Allen average template: %s\n',templateFile);
templateVolume = nrrdread(templateFile);
[templateVolume,templateResolution] = orientAllenVolume(templateVolume);

% SHARP-Track/SHARCQ convention at 10-um resolution: bregma is
% [540 0 570] voxels. Expressing this in millimeters makes the conversion
% independent of the downloaded template resolution.
bregmaMmFromOrigin = [5.4 0 5.7]; % [AP DV ML]
sites.apVoxel = round((bregmaMmFromOrigin(1)-sites.AP)./templateResolution);
sites.dvVoxel = round((bregmaMmFromOrigin(2)+sites.DV)./templateResolution);
sites.mlVoxel = round((bregmaMmFromOrigin(3)+sites.ML)./templateResolution);

useAnnotation = false;
annotationVolume = [];
annotationResolution = NaN;
if strlength(annotationFile)>0
    if isfile(annotationFile)
        fprintf('Loading Allen annotation volume: %s\n',annotationFile);
        annotationVolume = nrrdread(annotationFile);
        [annotationVolume,annotationResolution] = ...
            orientAllenVolume(annotationVolume);
        useAnnotation = true;
    else
        warning('Annotation file not found; omitting boundaries: %s', ...
            annotationFile);
    end
end

%% Draw the coronal panel grid

nPanels = numel(panelList);
nColumns = min(3,nPanels);
nRows = ceil(nPanels/nColumns);
fig = figure('Color','w','Units','inches', ...
    'Position',[0.5 0.5 3.0*nColumns 2.65*nRows]);
tl = tiledlayout(fig,nRows,nColumns, ...
    'TileSpacing','compact','Padding','compact');

includedColor = double(opt.IncludedColor);
excludedColor = double(opt.ExcludedColor);
excludedLineColor = 0.96.*excludedColor;

for panelNumber = 1:nPanels
    thisTick = panelTicks(panelNumber);
    thisAP = panelList(panelNumber);
    ax = nexttile(tl,panelNumber);

    apIndex = round((bregmaMmFromOrigin(1)-thisAP)./templateResolution);
    assert(apIndex>=1 && apIndex<=size(templateVolume,1), ...
        'AP %.3f mm maps outside the template (index %d).',thisAP,apIndex);

    templateSlice = squeeze(templateVolume(apIndex,:,:));
    imagesc(ax,templateSlice,opt.DisplayRange);
    axis(ax,'image');
    axis(ax,'off');
    set(ax,'YDir','reverse');
    colormap(ax,gray(256));
    hold(ax,'on');

    if useAnnotation
        annotationIndex = round( ...
            (bregmaMmFromOrigin(1)-thisAP)./annotationResolution);
        if annotationIndex>=1 && annotationIndex<=size(annotationVolume,1)
            labelSlice = squeeze(annotationVolume(annotationIndex,:,:));
            if ~isequal(size(labelSlice),size(templateSlice))
                labelSlice = resizeLabelsNearest( ...
                    labelSlice,size(templateSlice));
            end
            boundaryMask = labelBoundaryMask(labelSlice);
            boundaryImage = image(ax,ones([size(boundaryMask),3]));
            boundaryImage.AlphaData = 0.32.*boundaryMask;
        end
    end

    if logical(opt.ConnectPairs)
        animalsHere = unique( ...
            sites.pairLabel(sites.panelTick==thisTick),'stable');
        for animalNumber = 1:numel(animalsHere)
            useSite = sites.panelTick==thisTick & ...
                sites.pairLabel==animalsHere(animalNumber);
            if nnz(useSite)>=2
                pair = sortrows(sites(useSite,:),"mlVoxel");
                if all(pair.included)
                    lineColor = 0.72.*includedColor;
                else
                    lineColor = excludedLineColor;
                end
                plot(ax,pair.mlVoxel,pair.dvVoxel,'-', ...
                    'Color',lineColor,'LineWidth',0.9);
            end
        end
    end

    includedHere = sites.panelTick==thisTick & sites.included;
    if any(includedHere)
        scatter(ax,sites.mlVoxel(includedHere),sites.dvVoxel(includedHere), ...
            20,includedColor,'MarkerEdgeColor','w','LineWidth',0.8);
    end

    excludedHere = sites.panelTick==thisTick & ~sites.included;
    if any(excludedHere)
        scatter(ax,sites.mlVoxel(excludedHere),sites.dvVoxel(excludedHere), ...
            22,excludedColor,'x','LineWidth',1.4);
    end

    if ~isempty(opt.CropMm)
        crop = double(opt.CropMm);
        mlLimits = (bregmaMmFromOrigin(3)+crop(1:2))./templateResolution;
        dvLimits = (bregmaMmFromOrigin(2)+crop(3:4))./templateResolution;
        xlim(ax,mlLimits);
        ylim(ax,dvLimits);
    end

    title(ax,sprintf('AP %+0.2f mm',thisAP), ...
        'FontName','Arial','FontSize',10,'FontWeight','normal');
end

if strlength(string(opt.FigureTitle))>0
    sgtitle(tl,string(opt.FigureTitle), ...
        'FontName','Arial','FontSize',12,'FontWeight','bold');
end

if logical(opt.ShowLegend)
    firstAx = nexttile(tl,1);
    hold(firstAx,'on');
    legendHandles = gobjects(0);
    legendLabels = strings(0);
    if any(sites.included)
        legendHandles(end+1) = scatter(firstAx,nan,nan,42, ...
            includedColor,'filled','MarkerEdgeColor','w','LineWidth',0.8); %#ok<AGROW>
        legendLabels(end+1) = "Included"; %#ok<AGROW>
    end
    if any(~sites.included)
        legendHandles(end+1) = scatter(firstAx,nan,nan,45, ...
            excludedColor,'x','LineWidth',1.4); %#ok<AGROW>
        legendLabels(end+1) = "Excluded"; %#ok<AGROW>
    end
    legend(firstAx,legendHandles,cellstr(legendLabels), ...
        'Location','southoutside','Orientation','horizontal','Box','off');
end

%% Report and export

fprintf('  Animals plotted: %d\n',height(animalSummary));
fprintf('  Sites plotted:   %d\n',height(sites));
fprintf('  Raw AP range:    %+0.2f to %+0.2f mm\n', ...
    min(sites.AP),max(sites.AP));
fprintf('  Display order:   %s mm\n', ...
    strjoin(compose('%+0.2f',panelList),', '));

if logical(opt.SaveFigure)
    outputStem = string(opt.OutputStem);
    assert(isscalar(outputStem) && strlength(outputStem)>0, ...
        'OutputStem is required when SaveFigure is true.');
    outputFolder = string(fileparts(char(outputStem)));
    if strlength(outputFolder)>0 && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    exportFigurePair(fig,outputStem);
end

end


function included = normalizeInclusion(values,includedValues)
if islogical(values)
    included = values;
elseif isnumeric(values)
    included = isfinite(values) & values~=0;
else
    values = lower(strtrim(string(values)));
    includedValues = lower(strtrim(string(includedValues)));
    included = ismember(values,includedValues);
end
included = included(:);
end


function values = numericColumn(values)
if isnumeric(values) || islogical(values)
    values = double(values);
else
    values = str2double(string(values));
end
values = values(:);
end


function [volume,resolutionMm] = orientAllenVolume(volume)
% Return Allen volume dimensions as AP x DV x ML.
sz = size(volume);
assert(numel(sz)==3,'Expected a three-dimensional Allen volume.');
expectedRatio = 13.2/8.0;
if abs((sz(1)/sz(2))-expectedRatio) < 0.03
    % Already AP x DV x ML.
elseif abs((sz(2)/sz(1))-expectedRatio) < 0.03
    volume = permute(volume,[2 1 3]);
else
    error('Unrecognized Allen volume dimensions: %s',mat2str(sz));
end
resolutionMm = 13.2/size(volume,1);
end


function resized = resizeLabelsNearest(labels,targetSize)
rowIndex = round(linspace(1,size(labels,1),targetSize(1)));
columnIndex = round(linspace(1,size(labels,2),targetSize(2)));
resized = labels(rowIndex,columnIndex);
end


function boundaryMask = labelBoundaryMask(labels)
boundaryMask = false(size(labels));
horizontalChange = labels(:,2:end)~=labels(:,1:end-1);
verticalChange = labels(2:end,:)~=labels(1:end-1,:);
boundaryMask(:,2:end) = boundaryMask(:,2:end) | horizontalChange;
boundaryMask(:,1:end-1) = boundaryMask(:,1:end-1) | horizontalChange;
boundaryMask(2:end,:) = boundaryMask(2:end,:) | verticalChange;
boundaryMask(1:end-1,:) = boundaryMask(1:end-1,:) | verticalChange;
brainNeighborhood = labels>0;
brainNeighborhood(:,2:end) = ...
    brainNeighborhood(:,2:end) | labels(:,1:end-1)>0;
brainNeighborhood(2:end,:) = ...
    brainNeighborhood(2:end,:) | labels(1:end-1,:)>0;
boundaryMask = boundaryMask & brainNeighborhood;
end


function exportFigurePair(fig,outputStem)
pdfFile = outputStem+".pdf";
pngFile = outputStem+".png";

if exist('exportgraphics','file')==2
    exportgraphics(fig,pdfFile,'ContentType','vector');
    exportgraphics(fig,pngFile,'Resolution',600);
else
    % Compatibility path for MATLAB releases predating exportgraphics.
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

fprintf('  Saved PDF:       %s\n',pdfFile);
fprintf('  Saved PNG:       %s\n',pngFile);
end


function tf = isTextScalar(x)
tf = (ischar(x) && (isrow(x) || isempty(x))) || ...
    (isstring(x) && isscalar(x));
end


function tf = isLogicalScalar(x)
tf = (islogical(x) || isnumeric(x)) && isscalar(x);
end


function tf = isNumericScalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x);
end


function tf = isRgbTriplet(x)
tf = isnumeric(x) && numel(x)==3 && all(isfinite(x(:))) && ...
    all(x(:)>=0 & x(:)<=1);
end
