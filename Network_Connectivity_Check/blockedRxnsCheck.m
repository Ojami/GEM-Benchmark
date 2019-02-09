function blockedRxnsCheck(TRFBAcheck)
clc
% Identifies the fraction of blocked reactions for all the GEMs generated
% by each algorithm.
% TRFBAcheck: true if dealing with TRFBA or TRFBA-CORE GEMs, otherwise
%             false
% O. Jamialahmadi, Jan 2019.

if nargin < 1
    TRFBAcheck = false;
end

warning('The function is under developmenet, and will be added to the main GUI application')

currPath = which('blockedRxnsCheck.m');
currPath = strrep(currPath,'Network_Connectivity_Check\blockedRxnsCheck.m','');
addpath(genpath([currPath,'GEMs_Comparison'])); 

changeCobraSolver('gurobi5');
algorithmNames = {'pFBAc','PRIME','GIMME','iMAT','INIT',...
        'mCADRE','FASTCORE','FASTCORMICS','CORDA'};

NumLen = numel(algorithmNames); % Ignore pFBA
for i = 1:NumLen
    if ~strcmp(algorithmNames{i},'pFBAc')
        algorithmNames{end+1} = [algorithmNames{i},'c'];
    end
end
algorithmNames = sort(algorithmNames);

BlockedRxnsClosed = (0); BlockedRxnsOpen = (0);
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
    for inctr = 1:numel(GEMnames)
        fprintf('Model = %d (of %d)\n',inctr,numel(GEMnames))
        load(GEMnames{inctr})
        model = OutM; clear OutM
        if TRFBAcheck
            A = identifyFastBlockedRxns(model,model.rxns,0);
            A = modifyBlocked4TRFBA(A);
            BlockedRxnsClosed(inctr,count) = numel(A)./numel(model.rxns);
            clear A
            % Remove medium constraints------------------------------
            exRxns = model.rxns(strncmpi('ex_',model.rxns,3));
            model.ub(ismember(model.rxns,exRxns)) = max(abs(model.ub));
            %--------------------------------------------------------
             A = identifyFastBlockedRxns(model,model.rxns,0);
             A = modifyBlocked4TRFBA(A);
            BlockedRxnsOpen(inctr,count) = numel(A)./numel(model.rxns);
            clear A
        else
            A = fastcc(model, 1e-4,0);
            BlockedRxnsClosed(inctr,count) = (numel(model.rxns) - numel(A))./numel(model.rxns);
            clear A
            % Remove medium constraints------------------------------
            exRxns = model.rxns(strncmpi('ex_',model.rxns,3));
            model.lb(ismember(model.rxns,exRxns)) = -1*(max(abs(model.lb)));
            %--------------------------------------------------------
            A = fastcc(model, 1e-4,0);
            BlockedRxnsOpen(inctr,count) = (numel(model.rxns) - numel(A))./numel(model.rxns);
            clear A
        end
        % Use two approaches: 1- with and 2-without medium constraints
    
        clear model
        
    end
    fprintf('======================================\n')
end
rmpath(genpath([currPath,'GEMs_Comparison']));
save BlockedRxnsClosed BlockedRxnsClosed
save BlockedRxnsOpen BlockedRxnsOpen
end

function A = modifyBlocked4TRFBA(A)
frxns = regexp(A,'_f$');
frxns = ~cellfun('isempty',frxns);
frxns = A(frxns); frxns = regexprep(frxns,'_f$','');
brxns = regexp(A,'_b$');
brxns = ~cellfun('isempty',brxns);
brxns = A(brxns); brxns = regexprep(brxns,'_b$','');
brxns_diff = setdiff(brxns,frxns); brxns_diff = strcat(brxns_diff,'_b');
frxns_diff = setdiff(frxns,brxns); frxns_diff = strcat(frxns_diff,'_f');
A(ismember(A,union(brxns_diff,frxns_diff))) = [];
end