function CVanalyzer(algorithmNames)
clc

algorithmNames = strcat(algorithmNames,'CV');
metainfo = strcat(algorithmNames,'_metainfo');
for ct = 1:numel(metainfo)
    fprintf('CV analyzer for algorithm: %s (%d of %d)\n',strrep(algorithmNames{ct},'CV',''),ct,numel(metainfo))
    numModels = 1:75; % Default number for NCI-60 panel
    GEMnames = strcat(algorithmNames{ct},num2str(numModels'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');

    load(metainfo{ct}); % Indices of validation sets
    recoveredRxns = (0); modelSize = recoveredRxns; valid = modelSize; ctr1 = 1;
    for ctr = 1:size(Testindx,2)
        fprintf('count:%d of %d\n',ctr,size(Testindx,2))
        load(GEMnames{ctr})
        CurrTestindx = Testindx(:,ctr);
        CurrTestindx(~CurrTestindx) = []; % Remove zero indices
        if ~isempty(OutM.rxns) % Check the ability of algorithm to generate a model
%             % Build the consistent part of the OutM
%             Aactive = fastcc(OutM,1e-4,0);
%             OutM = removeRxns(OutM,OutM.rxns(setdiff(1:numel(OutM.rxns),Aactive)));
            valid(ctr1) = numel(core(CurrTestindx));
            recoveredRxns(ctr1) = sum(ismember(OutM.rxns,core(CurrTestindx))); % Recovered rxns
            modelSize(ctr1) = numel(OutM.rxns);
            ctr1 = ctr1 + 1;
        end
        clear OutM
    end
    clear Testindx OutM core
    Res(ct).meanValid = mean(valid);
    Res(ct).stdValid = std(valid);
    Res(ct).meanRecovered = mean(recoveredRxns);
    Res(ct).stdRecovered = std(recoveredRxns);
    Res(ct).meanModel = mean(modelSize);
    Res(ct).stdModel = std(modelSize);
    Res(ct).Hyperg = hygecdf(floor(mean(recoveredRxns))-1,2473,floor(mean(valid)),floor(mean(modelSize)),'upper');
    fprintf('.............................................\n')
end

algorithmNames = strrep(algorithmNames,'CV','');
Res = squeeze(cell2mat(struct2cell(Res)))';
disp(table(Res(:,5),Res(:,6),Res(:,3),Res(:,4),Res(:,1),Res(:,2),Res(:,7)...
     ,'VariableNames',{'GEM_mean','GEM_SD','Recovered_mean','Recovered_SD',...
     'Validation_mean','Validation_SD','p_value'},'RowNames',algorithmNames'))