function findEssGenes(algorithmNames,modelType,TRFBAcheck)
% Finds essential genes of generated GEMs
% Oveis Jamialahmadi
% Jan 2018

clc
if nargin < 3
    TRFBAcheck = false;
end
if nargin < 2
    modelType = 1;
end
if nargin < 1
    algorithmNames = {'pFBA','PRIME','GIMME','iMAT','INIT',...
        'mCADRE','FASTCORE','FASTCORMICS','CORDA'};
end

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

fprintf('Finding essential genes....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')
for count = 1:numel(algorithmNames)
    fprintf('%d-%s\n',count,algorithmNames{count})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');

alSize = (0);

if modelType == 1
    load Recon1_genes
    load consistentRecon1
    consistentRecon = consistentRecon1;
    clear consistentRecon1
elseif modelType == 2
    load Recon2_genes
    load consistentRecon2
    consistentRecon = model;
    clear model
% elseif modelType == 3 ######% For testing purposes
%     load Recon3_genes
%     Recon_genes = Recon3_genes;
%     clear Recon3_genes
%     load Recon3D
%     consistentRecon = model;
%     clear model
else
    error('Wrong modelType!')
end

for count = 1:numel(algorithmNames)
    if strcmp(algorithmNames{count},'TRFBA') || strcmp(algorithmNames{count},'TRFBAc')
        TRFBAcheck = true;
    else
        TRFBAcheck = false;
    end
    fprintf('------------------------------------------------\n')
    fprintf('%d-Processing for algorithm: %s\n',count,algorithmNames{count})
    numModels = 1:60; % Default number for NCI-60 panel
    GEMnames = strcat(algorithmNames{count},num2str(numModels'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');
    if ~exist(GEMnames{end},'file') % for 59 cell lines
        GEMnames(end) = [];
    end
    alSize(count) = numel(GEMnames);
    PotenGenes = ({});
    % Use parallel computing for TRFBA to speed-up the simulation
    if TRFBAcheck
        poolobj = parpool('local',3);
    end
    for cnt = 1:alSize(count)
        fprintf('Processing for cell line: %d of %d\n',cnt,alSize(count))
         load(GEMnames{cnt})
         model = OutM; clear OutM
        if TRFBAcheck
           grRatio = TRFBAessentiality(model,modelType,RegulatedGenes);
           clear RegulatedGenes
        else
            if numel(model.rules) ~= numel(model.rxns)
                % Modify the rules filed to avoid error in deleteModelGenes
                % function of COBRA toolbox
                idx = ismember(consistentRecon.rxns,model.rxns);
                model.rules = model.rules(idx);
            end
            model = checkObjective(model,modelType);
            [grRatio,~,~] = singleGeneDeletion (model,'FBA');
            %     list of model essential genes based on their effect on growth
            %     inhibition.
        end
        
        try
            Gr1loci = grRatio<=0.99; % Based on the criteria of Folger
            Gr1 = model.genes(Gr1loci);
            Gr2 = setdiff(Recon_genes,Gr1);
            PotenGenes{cnt,1} = Gr1; % Potential essential genes
            PotenGenes{cnt,2} = Gr2; % Other metabolic genes in Recon1
            PotenGenes{cnt,3} = grRatio;
        catch % No model is available
            PotenGenes{cnt,1} = 0;
            PotenGenes{cnt,2} = 0;
            PotenGenes{cnt,3} = 1;
        end
        save(['Potentialgenes',algorithmNames{count},'_',num2str(modelType)],'PotenGenes')
        clear model
    end
    if exist('poolobj','var')
        delete(poolobj)
    end
    
    fprintf('Enrichment analysis with cell line-specific essential genes...\n')
    MetabolicEssentiality_sp(PotenGenes,algorithmNames(count),modelType);
end

% =========================== Subfunctions ================================
function model = checkObjective(model,modelType)
% Check if biomass is available in the GEM and set it as obj func
if modelType == 1
    model.c(find(model.c)) = 0;
    if ~ismember('bioMassR',model.rxns)
        model = addReaction(model,'bioMassR', [' 0.118707 ala-L[c] +', ...
        '0.0009274 arg-L[c] + 0.0120562 asn-L[c] + 0.6927686 asp-L[c] +',...
        '0.0012056 cys-L[c] + 0.1910446 gln-L[c] + 0.1678596 glu-L[c] +',...
        '0.1391101 gly[c] + 0.0287494 his-L[c] + 0.0037096 ile-L[c] +',...
        '0.0111288 leu-L[c] + 0.009274 lys-L[c] + 0.0018548 met-L[c] +',...
        '0.0037096 phe-L[c] + 0.0179916 pro-L[c] + 0.03709604 ser-L[c] +',...
        '0.0204028 thr-L[c] + 0.0009274 trp-L[c] + 0.0055644 tyr-L[c] +',...
        '0.0120562 val-L[c] + 0.0022743 damp[c] + 0.0015162 dcmp[c] +',...
        '0.0015162 dgmp[c] + 0.0022743 dtmp[c] + 0.0092443 cmp[c] +',...
        '0.0104769 gmp[c] + 0.0055466 ump[c] + 0.0055466 amp[c] +',...
        '0.0024143 sphmyln_hs[c] + 0.0082708 chsterol[c] +',...
        '0.0065713 xolest_hs[c] + 0.00386 mag_hs[c] + 0.0031466 dag_hs[c] +',...
        '0.0030972 pail_hs[c] + 0.0200511 pe_hs[c] + 0.005891 ps_hs[c] +',...
        '0.0171707 pchol_hs[c] + 0.0012641 lpchol_hs[c] +',...
        '0.0005735 clpn_hs[c] + 0.016414 pa_hs[c] + 0.0274816 tag_hs[c] +',...
        '100 atp[c]  -> 100 adp[c] + 100 h[c] + 100 pi[c]']);
    end
    model = changeObjective(model,'bioMassR');
elseif modelType == 2
    if ~ismember('biomass_reaction',model.rxns)
        model = addReaction(model,'biomass_reaction',['20.6508 h2o[c] +',...
            ' 20.7045 atp[c] + 0.38587 glu_L[c] + 0.35261 asp_L[c] +',...
            ' 0.036117 gtp[c] + 0.27942 asn_L[c] + 0.50563 ala_L[c] + ',...
            '0.046571 cys_L[c] + 0.326 gln_L[c] + 0.53889 gly[c] + ',...
            '0.39253 ser_L[c] + 0.31269 thr_L[c] + 0.59211 lys_L[c] + ',...
            '0.35926 arg_L[c] + 0.15302 met_L[c] + 0.023315 pail_hs[c] + ',...
            '0.039036 ctp[c] + 0.15446 pchol_hs[c] + 0.055374 pe_hs[c] + ',...
            '0.020401 chsterol[c] + 0.002914 pglyc_hs[c] + 0.011658 clpn_hs[c] + ',...
            '0.053446 utp[c] + 0.009898 dgtp[n] + 0.009442 dctp[n] + ',...
            '0.013183 datp[n] + 0.013091 dttp[n] + 0.27519 g6p[c] + ',...
            '0.12641 his_L[c] + 0.15967 tyr_L[c] + 0.28608 ile_L[c] + ',...
            '0.54554 leu_L[c] + 0.013306 trp_L[c] + 0.25947 phe_L[c] + ',...
            '0.41248 pro_L[c] + 0.005829 ps_hs[c] + 0.017486 sphmyln_hs[c] + ',...
            '0.35261 val_L[c]  -> 20.6508 adp[c] + 20.6508 h[c] + 20.6508 pi[c]']);
    end
    model = changeObjective(model,'biomass_reaction');
else
    error('modelType must be 1, 2 or 3!')
end

function grRatio = TRFBAessentiality(model,modelType,RegulatedGenes)
if isempty(model)
    grRatio = [];
    return
end
s = optimizeCbModel(model,'max');
wildGR = s.f; clear s
grRatio = (1);
modelgenes = model.genes; modelmets = model.mets; modelrxnGeneMat = model.rxnGeneMat;
parfor count = 1:numel(modelgenes)
    fprintf('gene:%d of %d\n',count,numel(modelgenes))
    changeCobraSolver('gurobi5');
    modelT = model;
    if any(strcmp(RegulatedGenes,modelgenes{count}))
        nmets = numel(modelmets);
        GeneLoci = find(ismember(RegulatedGenes,modelgenes{count}));
        modelT = setModelb(modelT,nmets,GeneLoci);

    elseif any(modelrxnGeneMat(:,count))
        modelT = deleteModelGenes(modelT,modelgenes{count},0);
    end
    modelT = checkObjective(modelT,modelType);
    s = optimizeCbModel(modelT,'max');

    grRatio(count) = s.f/wildGR;
    s = (0);
end
grRatio(isnan(grRatio)) = 1; % To not select it as essential gene wrongly

function model = setModelb(model,nmets,GeneLoci)
model.b(GeneLoci+nmets) = 0;