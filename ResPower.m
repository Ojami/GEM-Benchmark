function ResPower(algorithmNames)
% Calculates similarity levels (proxy of resolution power) across
% NCI-60 cell line GEMs.
% 
% Oveis Jamialahmadi
% Jan 2018

clc
algorithmNames(ismember(algorithmNames,'pFBA')) = {'pFBAc'};
cellSpcCheck = questdlg('Do you want to test the cell-specific medium as well?(_c version)', ...
	'Cell-specific medium', ...
	'Yes','No','No');
if strcmp(cellSpcCheck,'Yes')
    % Add _c algorithms for cell-specific growth medium
    NumLen = numel(algorithmNames); % Ignore pFBA
    for i = 1:NumLen
        if ~strcmp(algorithmNames{i},'pFBAc')
            algorithmNames{end+1} = [algorithmNames{i},'c'];
        end
    end
    algorithmNames = sort(algorithmNames);
end
fprintf('Similarity check -Jaccard index-....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithms are going to be tested:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
%--------------------------------------------------------------------------
load CancerData4ResPower

alSize = (0);
for count = 1:numel(algorithmNames)
    fprintf('%d-Processing for algorithm: %s\n',count,algorithmNames{count})
    numModels = 1:60; % Default number for NCI-60 panel
    GEMnames = strcat(algorithmNames{count},num2str(numModels'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');
    if ~exist(GEMnames{end},'file') % for 59 cell lines
        GEMnames(end) = [];
    end
    alSize(count) = numel(GEMnames);
    if alSize(count) < 60
        shuff = shuff59;
        cancerNames = Cancer59;
        Idx = [1,6,7,14,15,21,22,27,28,37,38,44,45,52,53,54,55,59];
    else
        Idx = [1,6,7,15,16,22,23,28,29,38,39,45,46,53,54,55,56,60];
        shuff = shuff60;
        cancerNames = Cancer60;
    end 
    
    models = ({});
    for ct1 = 1:numel(GEMnames)
        load(GEMnames{ct1})
        models{ct1} = OutM;
        clear OutM
    end
    models = models(shuff);
    Pairwise_mat = zeros(numel(models),numel(models));
    for ct1 = 1:numel(models)
        for ct2 = 1:numel(models)
            % Jacard index
               Pairwise_mat(ct1,ct2) = ...
                sum(ismember(models{ct1}.rxns,models{ct2}.rxns))/...
                (numel(models{ct1}.rxns)+numel(models{ct2}.rxns)-...
                sum(ismember(models{ct1}.rxns,models{ct2}.rxns)));
        end
    end
    
    Pairwise_mat(isnan(Pairwise_mat)) = 0;
    M = zeros(9,9);
    m1 = 1;
    for i1 = 1:2:numel(Idx)
        m2 = 1;
        for i2 = 1:2:numel(Idx)
            M(m1,m2) = sum(sum(Pairwise_mat(Idx(i1):Idx(i1+1),Idx(i2):Idx(i2+1))))/...
                numel(Pairwise_mat(Idx(i1):Idx(i1+1),Idx(i2):Idx(i2+1)));
            m2 = m2 + 1; 
        end
        m1 = m1 + 1;
    end
    figure;
    colormap('jet')
    imagesc(M);
    colorbar
%     save(algorithmNames{count},'M')
    set(gca,'YTick',1:numel(Idx)/2)
    set(gca,'XTick',1:numel(Idx)/2)
    set(gca,'YTickLabel',unique(cancerNames,'stable'))
    set(gca,'XTickLabel',unique(cancerNames,'stable'))

	set(gcf,'name',algorithmNames{count},'numbertitle','off')
    set(gca, 'Ticklength', [0 0],'FontSize',11)
    savefig([algorithmNames{count},'_similarityLevel'])
end