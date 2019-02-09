function MetChecker(algorithmNames,modelType)
% Evaluates the ability of generated GEMs to predict metabolite uptake/secretion
% rates measured in the work of Jain et al.(2012), denoted by CORE data
% 
% Oveis Jamialahmadi,
% Jan 2018

clc
if nargin < 2
    modelType = 1;
end
if nargin < 1
    algorithmNames = {'pFBA','TRFBA','PRIME','GIMME','iMAT','INIT',...
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
    
fprintf('Simulating metabolite uptake/secretion rates....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');

% Metabolite check =================================================
load('MetCheckNeeds')
AllRes = ({});
alSize = (0);
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
   
    allExMets = ({});
	
    for k = 1:5%numel(sharedMets)
        fprintf('Metabolite: %d-%s (of %d)\n',k,sharedMets{k},numel(sharedMets))
        exMets = zeros(numel(GEMnames),1);
        for j = 1:numel(GEMnames) % Growth simulations
            load(GEMnames{j})
            model = checkObjective(OutM,modelType);
            if strcmp('pFBA',algorithmNames{count}) % minimize L1-norm for pFBA
                Sol = optimizeCbModel(model,'max','one');
            else
                Sol = optimizeCbModel(model,'max');
            end
            biomassFlux = Sol.f;
            clear Sol
            if strcmp(algorithmNames{count},'TRFBA') || strcmp(algorithmNames{count},'TRFBAc')
%                 model.rxns = regexprep(model.rxns,'\_f$','');
%                 model.rxns = regexprep(model.rxns,'\_r$','');
%                 model.rxns = regexprep(model.rxns,'\_b$','');
                findExchange = find(~cellfun('isempty',strfind(model.rxns,sharedRxns{k})));
                if isempty(findExchange)%Check for other rxn symbols
                     tempRxn1 = sharedRxns{k};
                     tempRxn1 = strrep(tempRxn1,'(','_');
                     tempRxn1 = strrep(tempRxn1,')','');
                     tempRxn1 = strrep(tempRxn1,'_L','__L');
                     findExchange = find(~cellfun('isempty',strfind(model.rxns,tempRxn1)));
                end
                    
                if numel(findExchange) > 2 % Something goes wrong 
                    error('Wrong check!')
                end

            else
                findExchange = find(ismember(model.rxns,sharedRxns{k}));
            end
            
            if ~isempty(findExchange) % can find the associated exchange rxn in the model
                model.lb(find(model.c)) = 0.9*biomassFlux;
                model.c(find(model.c)) = 0;
                if strcmp(algorithmNames{count},'TRFBA') || strcmp(algorithmNames{count},'TRFBAc')
                    checkRxns = model.rxns(findExchange);
                    checkRxnsIdx = ~cellfun('isempty',regexp(checkRxns,'\_f$'));
                    model.c(ismember(model.rxns,checkRxns(checkRxnsIdx))) = 1; % Forward
                    model.c(ismember(model.rxns,checkRxns(~checkRxnsIdx))) = -1; % Backward
                    
                else
                    model.c(findExchange) = 1;
                end
                if strcmp('pFBA',algorithmNames{count})
                    metFlx = optimizeCbModel(model,'max','one');
                else
                    metFlx = optimizeCbModel(model,'max');
                end
                exMets(j) = metFlx.f;
            else
                fprintf('Metabolite %s does not exist in the model:%d\n',sharedMets{k},j)
                exMets(j) = 1e10;
            end

        end
        
        allExMets{k} = exMets; 
        clear model 
    end
     
% Deal with simulated results ---------------------------------------------
    rho = nan(1,numel(allExMets)); p_val = nan(1,numel(allExMets));
    for ct = 1:numel(allExMets)
        exMets = allExMets{ct};
        exMets(45) = [];% remove 'MDA-N'
        if  size(exMets,1) <  size(exMets,2)
            exMets = exMets';
        end
        % remove 'No match' cell lines from data
        find10 = find(exMets==1e10);
        exMets(find10) = [];
        if isempty(exMets)
            continue
        end
        curmets = allMetabolites{ct};
        if alSize(count) < 60
            curmets(51) = []; % Rmove NCI-H23 for mCADRE/FASTCORE/FASTCORMICS
        end
        curmets(find10) = [];
        [rho(ct),p_val(ct)] = corr(exMets,curmets,'type','Spearman');
    end
    mod_rho = rho(~isnan(rho));
    mod_p_val = p_val(~isnan(rho));
    mod_sharedMets = sharedMets(~isnan(rho));
    mod_rho1 = mod_rho(mod_rho~=0);
    mod_p_val1 = mod_p_val(mod_rho~=0);
    mod_sharedMets1 = mod_sharedMets(mod_rho~=0);
    outdata1 = ({});
    for i1 = 1:numel(mod_rho1)
        outdata1{i1,1} = mod_rho1(i1);
        outdata1{i1,2} = mod_p_val1(i1);
        outdata1{i1,3} = mod_sharedMets1{i1};
    end
    if isempty(outdata1)
        AllRes{count,1} = algorithmNames{count};
        AllRes{count,2} = 'NA';
        continue
    end
    cell2table(outdata1,'VariableNames',{'Rho','P_value','Metabolite'})
    % Taken from PRIME:
    %Correcting for multiple hypothesis using FDR and significance level of
    alpha = 0.05;
    [pID,~] = FDR(mod_p_val1,alpha);
    if isempty(pID)
        fprintf('------------------------------------------------\n')
        fprintf('No significant metabolite pairs after correcting for muliple comparisons!');
        fprintf('------------------------------------------------\n')
        AllRes{count,1} = algorithmNames{count};
        AllRes{count,2} = 'NA';
        continue;
     end

    %Identifying the set of significantly correlated metabolites
    [p_sort,idx_sort] = sort(mod_p_val1);
    sig_p = find(p_sort<=pID);
    idx = idx_sort(1:numel(sig_p));
    mod_p_valC = mod_p_val1(idx);
    mod_rhoC = mod_rho1(idx);
    mod_sharedMetsC = mod_sharedMets1(idx);
    outdataC = ({});
    for i2 = 1:numel(mod_rhoC)
        outdataC{i2,1} = mod_rhoC(i2);
        outdataC{i2,2} = mod_p_valC(i2);
        outdataC{i2,3} = mod_sharedMetsC{i2};
    end
    cell2table(outdataC,'VariableNames',{'Rho','Corrected_P_value','Metabolite'})
    save(['Metcheck_',algorithmNames{count},'_',num2str(modelType)],'outdataC')
    AllRes{count,1} = algorithmNames{count};
    AllRes{count,2} = outdataC;
    clear outdataC
    fprintf('------------------------------------------------\n')
end
save('AlResults4Exometabolomics','AllRes') % Output for all selected algorithms 

%-----------------------------Subfunctions---------------------------------
% ============================= subfunctions ==============================
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
elseif modelType == 3
    if ~ismember('BIOMASS_reaction',model.rxns)
        bmassF =  ['20.6508 h2o[c] +',...
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
                '0.35261 val_L[c]  -> 20.6508 adp[c] + 20.6508 h[c] + 20.6508 pi[c]'];
        bmassF = strrep(bmassF,'_','__');
        bmassF = strrep(bmassF,'__hs','_hs');
        bmassF = strrep(bmassF,'[','_');
        bmassF = strrep(bmassF,']','');
        model = addReaction(model,'BIOMASS_reaction',bmassF);
    else
        model = changeObjective(model,'BIOMASS_reaction');
    end
else
    error('modelType must be 1 or 2!')
end

function [pID,pN] = FDR(p,q)
% This function has been taken from PRIME main paper
% FORMAT pt = FDR(p,q)
%
 % p   - vector of p-values
 % q   - False Discovery Rate level
 %
 % pID - p-value threshold based on independence or positive dependence
 % pN  - Nonparametric p-value threshold
 %______________________________________________________________________________
 % @(#)FDR.m    1.3 Tom Nichols 02/01/18
 
 p = sort(p(:));
 V = length(p);
 I = (1:V)';
 
 
 cVID = 1;
 cVN = sum(1./(1:V));
 
 pID = p(max(find(p<=I/V*q/cVID)));
 pN = p(max(find(p<=I/V*q/cVN)));