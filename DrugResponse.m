function DrugResponse(dataSet,algorithmNames,modelType)
% Performs drug response simulation for generated GEMs using 9 algorithmNames 
% under test + pFBAc

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

fprintf('Drug response simulations....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be testes:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');
alSize = (0);

% Read Drug data from three different datasets in order--------------------
load(dataSet)
if ~exist('removeIdx','var')
    removeIdx = 0;
end
if ~exist('nexp','var')
    nexp = 0;
end
% Get drug targets for drugs listed in DrugBank database
% load drugTargets

load DBankData
if modelType == 1
    load Recon1_genes
elseif modelType == 2
    load Recon2_genes
else
    error('Wrong modelType!')
end
Drugname = unique(allDrugs); drugTargets = cell(1,numel(Drugname));
for cnt = 1:numel(Drugname)
    drugTargets{cnt} = findDBankTargets(Drugname{cnt},Recon_genes,DBankData);
end

for count = 1:numel(algorithmNames)
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
    
    % Define statistics variables:
    sigDrug_idx = (0); drgCtr = 1; % Significant drugs for any p-value <0.05
    mean_sigDrug_idx = (0); drgCtrm = 1; % Significant drugs for mean p-value <0.05
    sig_drugname = ({}); mean_sig_drugname = ({});
    stat = struct('r',0,'p',1,'fr',0,'fp',1);
    sig_stat = struct('r',0,'p',1,'fr',0,'fp',1);
    mean_sig_stat = struct('r',0,'p',1,'fr',0,'fp',1);
    
    for ct = 1:numel(drugTargets)
        SimulatedFluxesT = (0); 
        if isempty(drugTargets{ct})
            SimulatedFluxesT = 'None';
            fprintf('Drug:%d of %d was overlooked\n',ct,numel(drugTargets))
            continue
        end
        this_drug = drugTargets{ct}; 
        fprintf('Drug:%d of %d:%s\n',ct,numel(drugTargets),Drugname{ct})
        for i = 1:alSize(count)
%             fprintf('Drug:%d of %d:%s | model count = %d\n',ct,numel(drugTargets),Drugname{ct},i)
            load(GEMnames{i})
            model = OutM; clear OutM
            model = checkObjective(model,modelType);
            if strcmp('pFBA',algorithmNames{count}) % minimize L1-norm for pFBA
                Flx = optimizeCbModel(model,'max','one');
            else
                Flx = optimizeCbModel(model,'max');
            end
            
            bmass = find(model.c);
            model.ub(bmass) = 0.5*Flx.f;
            model.lb(bmass) = 0.5*Flx.f;
            model.c(bmass) = 0;
            where_genes = ismember(model.genes,this_drug);
            [where_rxns,~] = find(model.rxnGeneMat(:,where_genes));
            where_rxns = unique(where_rxns);
            if isempty(where_rxns)
                SimulatedFluxesT(i,1) = eps; 
            end
            if strcmp(algorithmNames{count},'TRFBA') || strcmp(algorithmNames{count},'TRFBAc')
                checkRxns = model.rxns(where_rxns);
                checkRxnsUt = regexprep(checkRxns,'\_f$','');
                checkRxnsUt = regexprep(checkRxnsUt,'\_b$','');
                checkRxnsU = unique(checkRxnsUt);
                for ch1 = 1:numel(checkRxnsU)
                    if any(strcmp(model.rxns,checkRxnsU{ch1})) % Irreversible rxn
                        model.c(strcmp(model.rxns,checkRxnsU{ch1})) = 1;
                        Flx = optimizeCbModel(model,'max');
                        if Flx.f < eps
                            % Avoid confusing with non available elements
                            SimulatedFluxesT(i,ch1) = eps; 
                        else
                            SimulatedFluxesT(i,ch1) = Flx.f;
                        end
                        model.c(find(model.c)) = 0;
                    else % Reversible rxns: fwd obj = 1, bkw obj = -1
                        RevRxnsTemp = checkRxns(ismember(checkRxnsUt,checkRxnsU{ch1}));
                        fwdRxn = RevRxnsTemp(~cellfun('isempty',regexp(RevRxnsTemp,'\_f$')));
                        wherex_f = ismember(model.rxns,fwdRxn);
                        model.c(wherex_f) = 1;
                        bkwRxn = RevRxnsTemp(~cellfun('isempty',regexp(RevRxnsTemp,'\_b$')));
                        wherex_b = ismember(model.rxns,bkwRxn);
                        model.c(wherex_b) = -1;
                        Flx = optimizeCbModel(model,'max');
                        if Flx.f < eps
                            % Avoid confusing with non available elements
                            SimulatedFluxesT(i,ch1) = eps; 
                        else
                            SimulatedFluxesT(i,ch1) = Flx.f;
                        end
                        model.c(find(model.c)) = 0;
                    
                    end   
                end
            else % Other algorithmNames----------------------------------------
                for i2 = 1:numel(where_rxns)
                    model.c(where_rxns(i2)) = 1;
                    if strcmp('pFBA',algorithmNames{count})
                        Flx = optimizeCbModel(model,'max','one');
                    else
                        Flx = optimizeCbModel(model,'max');
                    end
                    if Flx.f < eps
                        % Avoid confusing with non available elements
                        SimulatedFluxesT(i,i2) = eps; 
                    else
                        SimulatedFluxesT(i,i2) = Flx.f;
                    end
                    model.c(where_rxns(i2)) = 0;
                end
            end
            clear model
        end
        if all(~SimulatedFluxesT)
            SimulatedFluxesT = 'None';
        end
    % ======================= IC50 analysis ===============================
    [S,sigName,sigS,meanSigName,meanSigS] = IC50analyzer(SimulatedFluxesT,...
        allIC50,allDrugs,Drugname,ct,nexp,removeIdx,alSize(count));
    stat(ct) = S; sig_drugname{ct} = sigName; sig_stat(ct) = sigS; 
    mean_sig_drugname{ct} = meanSigName; mean_sig_stat(ct) = meanSigS;       
    end
     fprintf('-----------------------------------------------------------\n')
     %# Commented: Only for testing purposes
%     save(['DrugResponse_',dataSet,algorithmNames{count},'.mat'],'stat','sig_drugname','mean_sig_drugname',...
%             'mean_sig_stat','sig_stat')
    % Export significant data
    testR = {mean_sig_stat(:).fr}; testP = {mean_sig_stat(:).fp};
    testDrugs = mean_sig_drugname(~cellfun(@isempty,testR));
    testR = testR(~cellfun(@isempty,testR)); testR = abs([testR{:}]);
    testP = testP(~cellfun(@isempty,testP)); testP = [testP{:}];
    testDrugs = testDrugs(testR>0);
    testDrugs = [testDrugs{:}];
    testDrugs = testDrugs';
    testP = testP(testR>0)';
    testR = testR(testR>0)';
    exportData = [{'Drug','R','p_value'};[testDrugs,...
        strtrim(cellstr(num2str(testR))),strtrim(cellstr(num2str(testP)))]];
    save(['DrugResponse_',dataSet,algorithmNames{count},'.mat'],'exportData')
    clear testDrugs exportData testR testP
end
% ===========================Subfunctions=================================
function [stat,sig_drugname,sig_stat,mean_sig_drugname,mean_sig_stat] = IC50analyzer(SimulatedFluxesT,...
    allIC50,allDrugs,Drugname,count,nexp,removeIdx,alSize)
if strcmp(SimulatedFluxesT,'None')
    stat.r = 0; stat.p = 1;
    stat.fr = 0; stat.fp = 1;
    sig_stat.r = 0; sig_stat.p = 1;
    sig_stat.fr = 0; sig_stat.fp = 1;
    mean_sig_stat.r = 0; mean_sig_stat.p = 1;
    mean_sig_stat.fr = 0; mean_sig_stat.fp = 1;
    sig_drugname = 'NA';
    mean_sig_drugname = 'NA';
    return
end

if alSize < 60
    allIC50(52,:) = []; % NCI-H23
    if removeIdx % Remove different cell lines from NCI-60 panel
        removeIdx(removeIdx>52) = removeIdx(removeIdx>52) - 1;
        SimulatedFluxesT(removeIdx,:) = [];
    end
else
    if removeIdx % Remove different cell lines from NCI-60 panel
        SimulatedFluxesT(removeIdx,:) = [];
    end
end

drgCtr = 1; drgCtrm = 1;
rho = (0); p_val = (0);
thisDrugIdx = find(ismember(allDrugs,Drugname{count}));
f2 = 1; found_ic50 = (0); where_nan = ({}); found_loci = (0);
for f1 = 1:numel(thisDrugIdx) % find non-available ic50 data
    where_nan{f1} = find(isnan(allIC50(:,thisDrugIdx(f1))));
    if isempty(where_nan{f1}) % all data available, so choose this one
        found_ic50(1:size(allIC50,1),f2) = allIC50(:,thisDrugIdx(f1));
        found_loci(f2) = thisDrugIdx(f1);
        f2 = f2 + 1;
    end
end
if f2 == 1
    where_nan_numel = cellfun(@numel,where_nan);
    [~,nan_index] = min(where_nan_numel);
    fprintf('IC50 values for some cell lines(%d) are not availble::Drug:%s\n',min(where_nan_numel),Drugname{count})
    fprintf('Those cell lines will be removed for the mentioned drug!\n')
    found_ic50 = allIC50(:,thisDrugIdx(nan_index));
    found_ic50(where_nan{nan_index}) = [];
end

if nexp(~isnan(nexp))
    if size(found_ic50,2) > 1  % If there are duplicate IC50 values for a drug
        found_loci1 = nexp(found_loci); % chosse the one with highest number of experiments
        [~,loc1] = max(found_loci1);
        how_loc1 = find(found_loci1(loc1) == found_loci1);
        if numel(how_loc1) > 1
            found_ic50 = mean(found_ic50,2);
        else
            found_ic50 = found_ic50(:,loc1);
        end
    end
else
    found_ic50 = mean(found_ic50,2);
end
% Do comparison for each rxn in all models of a certain algorithm
for i2 = 1:size(SimulatedFluxesT,2)
    curr_rxn = SimulatedFluxesT(:,i2);
    if f2 == 1 % some cell lines' fluxes should be removed
        curr_rxn(where_nan{nan_index}) = [];
    end
    if size(curr_rxn,1) < size(curr_rxn,2)
        curr_rxn = curr_rxn';
    end
    [rho(i2),p_val(i2)] = corr(curr_rxn,found_ic50,'type','Spearman');
end
rhoNet = rho(~isnan(rho));
p_valNet = p_val(~isnan(rho));
stat.r = rhoNet;
stat.p = p_valNet;
stat.fr = mean(rhoNet); stat.fp = mean(p_valNet);
Significant_pval = p_valNet(p_valNet<0.05);
if ~isempty(Significant_pval)
    sigDrug_idx(drgCtr) = count;
    drgCtr = drgCtr + 1;
end
if mean(p_valNet) < 0.05 % Mean p-value for all rxns corresponding to a certain drug
    mean_sigDrug_idx(drgCtrm) = count;
    drgCtrm = drgCtrm + 1;
end
if drgCtr > 1
    sig_drugname = Drugname(sigDrug_idx);
    sig_stat = stat;
else
    sig_drugname = 'NA';  sig_stat.r = 0; sig_stat.p = 1;
    sig_stat.fr = 0; sig_stat.fp = 1;
end
if drgCtrm > 1
    mean_sig_drugname = Drugname(mean_sigDrug_idx);
    mean_sig_stat = stat;
else
    mean_sig_drugname = 'NA';  mean_sig_stat.r = 0; mean_sig_stat.p = 1;
    mean_sig_stat.fr = 0; mean_sig_stat.fp = 1;
end
        
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
    error('modelType must be 1 or 2!')
end
