function RxnExpVals = GPRmapper(ParsedGPR,corrRxn,GeneNms,GeneExp,model)
% Converts GPRs from GPRparser to reaction expression values using min (for
% AND) and max (for OR) logic.

% -----------------------Map expression data using GRP rules---------------
AllGenesExp_temp = (0);
for i = 1:size(ParsedGPR,1)
    TheseGenes = ParsedGPR(i,:);
    TheseGenes(ismember(TheseGenes,'')) = [];
    TheseGenesExp = (0);
    for j = 1:numel(TheseGenes)
        findEachGene = GeneExp(ismember(GeneNms,TheseGenes{j}));
        if isempty(findEachGene) % Model genes are not in expression data
            TheseGenesExp(j) = 0;
        else
            findEachGene = mean(findEachGene);
            TheseGenesExp(j) = findEachGene;
        end
    end
    AllGenesExp_temp(i) = min(TheseGenesExp); % AND GPR associations
end
% Decide on OR GPR associations
AllGenesExp = (0); ctr = 1; AllRxns = ({});
corrRxn_temp = corrRxn;
while numel(corrRxn_temp)
    ThisRxn = corrRxn_temp{1};
    findThisRxn = find(ismember(corrRxn_temp,ThisRxn));
    findGeneExp = AllGenesExp_temp(findThisRxn);
    AllGenesExp(ctr) = max(findGeneExp); % OR
    AllRxns{ctr} = ThisRxn;
    corrRxn_temp(findThisRxn) = [];
    AllGenesExp_temp(findThisRxn) = [];
    ctr = ctr + 1;
end
AllRxns(~AllGenesExp) = [];
AllGenesExp(~AllGenesExp) = [];

RxnExpVals = zeros(numel(model.rxns),1);
[~,idx1] = intersect(model.rxns,AllRxns);
RxnExpVals(idx1) = AllGenesExp(1:numel(AllRxns));
