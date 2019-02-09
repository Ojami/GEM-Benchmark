function growthEval_CV(algorithmNames,modelType)
% Evaluates the relative error for growth prediction ability of generated
% GEMs using 9 algorithms under test + pFBAc
% Oveis Jamialahmadi
% Jan 2018

clc
folderNames = strcat(algorithmNames,'_CV');
checkExist = (0);
for i = 1:numel(folderNames)
    checkExist(i) = exist(folderNames{i},'file');
end
NotExistAlg = algorithmNames(~checkExist);
if ~isempty(NotExistAlg)
    fprintf('GEMs of the following algorithms are not present in CV folder:\n')
    for i = 1:numel(NotExistAlg)
        fprintf('%d- %s\n',i,NotExistAlg{i})
    end
    fprintf('------------------------------------------------\n')
end
algorithmNames = algorithmNames(logical(checkExist));
algorithmNames = strcat(algorithmNames,'CV');
fprintf('Growth prediction evaluation for algorithms....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')

for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');

load NCImetainfo
algorithmNames = sort(algorithmNames);

alSize = (0);
eval_error = 100.*ones(75,numel(algorithmNames));
pred_gr = zeros(75,numel(algorithmNames));
for count = 1:numel(algorithmNames)
%     if strcmp(algorithmNames{count},'TRFBACV')
%         continue
%     end
    fprintf('------------------------------------------------\n')
    fprintf('%d-Processing for algorithm: %s\n',count,algorithmNames{count})
    expGR = Measured_gr(2);

    GEMnames = strcat(algorithmNames{count},num2str([1:75]'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');

    alSize(count) = numel(GEMnames);
    for j = 1:numel(GEMnames) % Growth simulations
        fprintf('GEM:%d of %d | %d-%s\n',j,numel(GEMnames),count,algorithmNames{count})
        load(GEMnames{j})
        OutM = checkObjective(OutM,modelType);        
        
        if strcmp('pFBA',algorithmNames{count}) % minimize L1-norm for pFBA
            Sol = optimizeCbModel(OutM,'max','one');
        else
            Sol = optimizeCbModel(OutM,'max');
        end
        if ~isempty(Sol.x) %able to simulate
            if Sol.f
                eval_error(j,count) = abs(Sol.f - expGR) / expGR;
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
eval_error(:,uniqueMem(countMem == 75)) = []; % Remove algorithms not capable of predicting growth
pred_gr(:,uniqueMem(countMem == 75)) = []; % Remove algorithms not capable of predicting growth
algorithmNames(uniqueMem(countMem == 75)) = [];
alSize(uniqueMem(countMem == 75)) = [];

modelCount = (0);
for i = 1:numel(alSize)
    currEr = eval_error(:,i);
    currGr = pred_gr(:,i);
    currEr = currEr(1:alSize(i));
    currGr = currGr(1:alSize(i));
    whereNA = currEr == 100;
    if ~sum(whereNA)
        modelCount(i) = alSize(i);
    else
        modelCount(i) = sum(whereNA);
        currEr(whereNA) = [];
        currGr(whereNA) = [];
    end
    figure
    X = 1:numel(currEr);
    Y = currGr./max(currGr);
%     save(algorithmNames{i},'X','Y')
    plot(1:numel(currGr),currGr./max(currGr),'k.-','LineWidth',2,'Marker','O','MarkerFaceColor','k')
    xlabel('Model count','FontName','Times New Roman','FontSize',14);
    ylabel('Rleative error','FontName','Times New Roman','FontSize',14);
    title(algorithmNames{i})
    set(gca,'FontSize',10,'FontWeight','bold')
end

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

