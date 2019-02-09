function growthEval(algorithmNames,modelType)
% Evaluates the relative error for growth prediction ability of generated
% GEMs using 9 algorithms under test + pFBAc
% 
% Oveis Jamialahmadi
% Jan 2018

clc
if nargin < 2
    modelType = 1;
end
if nargin < 1
    algorithmNames = {'pFBA','TRFBA','PRIME','GIMME','iMAT','INIT',...
        'mCADRE','FASTCORE','FASTCORMICS','CORDA'};
end

fprintf('Growth prediction evaluation for algorithms....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');
load NCImetainfo

% Add _c algorithms for cell-specific growth medium
algorithmNames(ismember(algorithmNames,'pFBA')) = {'pFBAc'};
NumLen = numel(algorithmNames); % Ignore pFBA
for i = 1:NumLen
    if ~strcmp(algorithmNames{i},'pFBAc')
        algorithmNames{end+1} = [algorithmNames{i},'c'];
    end
end
algorithmNames = sort(algorithmNames);
alSize = (0);
eval_error = 100.*ones(60,numel(algorithmNames));
pred_gr = zeros(60,numel(algorithmNames));
for count = 1:numel(algorithmNames)
    fprintf('------------------------------------------------\n')
    fprintf('%d-Processing for algorithm: %s\n',count,algorithmNames{count})
    expGR = Measured_gr;
    numModels = 1:60; % Default number for NCI-60 panel
    GEMnames = strcat(algorithmNames{count},num2str(numModels'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');
    if ~exist(GEMnames{end},'file') % for 59 cell lines
        GEMnames(end) = [];
        expGR(52) = [];
    end
    alSize(count) = numel(GEMnames);
    for j = 1:numel(GEMnames) % Growth simulations
        fprintf('GEM:%d of %d | %d-%s\n',j,numel(GEMnames),count,algorithmNames{count})
        load(GEMnames{j})
        if isempty(OutM)
            eval_error(j,count) = 100;
            pred_gr(j,count) = 0;
            continue
        end
        OutM = checkObjective(OutM,modelType);
        if strcmp('pFBA',algorithmNames{count}) % minimize L1-norm for pFBA
            Sol = optimizeCbModel(OutM,'max','one');
        else
            Sol = optimizeCbModel(OutM,'max');
        end
        if ~isempty(Sol.x) %able to simulate
            if Sol.f
                eval_error(j,count) = abs(Sol.f - expGR(j)) / expGR(j);
                pred_gr(j,count) = Sol.f;
            end
        end
        clear OutM
    end
    clear GEMnames
    fprintf('------------------------------------------------\n')
end

% Post-processing the evaluated errors
[~,rmvIndx] = find(eval_error == 100);
uniqueMem = unique(rmvIndx);
countMem = histc(rmvIndx,uniqueMem);
eval_error(:,uniqueMem(countMem == 60)) = []; % Remove algorithms not capable of predicting growth
algorithmNames(uniqueMem(countMem == 60)) = [];
alSize(uniqueMem(countMem == 60)) = [];
pred_gr(:,uniqueMem(countMem == 60)) = [];
modelCount = (0); corrR = (0); corrP = (0);
for count = 1:numel(algorithmNames)
    currEr = eval_error(:,count);
    currGr = pred_gr(:,count);
    expGR = Measured_gr;
    whereNA = currEr == 100; whereNA = whereNA(1:alSize(count));
    if ~sum(whereNA) % all available!
        if alSize(count) < 60
            currEr(end) = currEr(1); % Only for boxplot function of MATLAB
            expGR(52) = [];
            currGr(end) = [];
        end
        modelCount(count) = numel(whereNA);
        [corrR(count),corrP(count)] = corr(expGR,currGr,'Type','Spearman');
    else % Some model cannot predict growth
        validGR = currEr(~whereNA);
        
        currEr(whereNA) = validGR(1); % Only for boxplot function of MATLAB
        if alSize(count) < 60
            currEr(end) = validGR(1); % Only for boxplot function of MATLAB
            expGR(52) = [];
            currGr(end) = [];
            currGr = currGr(~whereNA);
            expGR = expGR(~whereNA);
        else
            currGr = currGr(~whereNA);
            expGR = expGR(~whereNA);
        end
        modelCount(count) = numel(validGR);
        [corrR(count),corrP(count)] = corr(expGR,currGr,'Type','Spearman');
    end
    eval_error(:,count) = currEr;
end
disp(table(corrR',corrP','VariableNames',{'R','p_value'},'RowNames',algorithmNames'))
fprintf('------------------------------------------------\n')
disp(table(modelCount',alSize','VariableNames',{'Number_of_Capable_GEMs','All_models'},'RowNames',algorithmNames'))
figure;
eval_error(eval_error==100) = NaN;
boxplot(eval_error,algorithmNames,'symbol','kO','labelorientation','inline','colors','k',...
    'outliersize',4)
ylabel('Relative error')
set(gca,'FontSize',10,'FontWeight','bold')
set(findobj(gca,'Type','text'),'FontSize',10,'FontWeight','bold')
set(findobj(gca,'Type','line'),'LineWidth',1.3)

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
else
    error('modelType must be 1 or 2!')
end