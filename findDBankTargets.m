function Drugtargets = findDBankTargets(Drugname,modelGenes,DBankData)

dbankID = DBankData.dbankID;
dbankGenes = DBankData.dbankGenes;
dbankIDs2genes = DBankData.dbankIDs2genes;
dbanknames = DBankData.dbanknames;
drugSyn = DBankData.drugSyn;
entrez_id = DBankData.entrez_id;
gene_symbol = DBankData.gene_symbol;

% Search for input drug in DrugBank =======================================
DrugID = dbankID(ismember(lower(dbanknames),lower(Drugname))); %Find Drug ID
if isempty(DrugID)
    drugSyn = strcat('|',{' '},drugSyn,' |');
    chckStr = ['| ',lower(Drugname),' |'];
    isSyn = strfind(lower(drugSyn),chckStr);
    isSyn = ~cellfun('isempty',isSyn);
    if ~any(isSyn)
        fprintf('- %s is not in DrugBank database!\n',Drugname)
        fprintf('---------------------------------------------\n')
        Drugtargets = '';
        return
    else
        DrugID = dbankID(isSyn);
    end
end
if ~ischar(DrugID) % For strfind input
    DrugID = DrugID{1};
end
dIDX = strfind(dbankIDs2genes,DrugID); dIDX = ~cellfun('isempty',dIDX);
if isempty(find(dIDX, 1))
    fprintf('%s has no targets in DrugBank!\n',Drugname)
    fprintf('---------------------------------------------\n')
    Drugtargets = '';
    return
end
target_genes = dbankGenes(dIDX); % Find gene names for all drug targets

% check if some Entrez ids could not be found
NotfoundEntrez = target_genes(~ismember(target_genes,gene_symbol));
if ~isempty(NotfoundEntrez)
    fprintf('Some Entrez IDs could not be found for drug:%s\n',Drugname)
    fprintf('Total targets for this drug:%d\n',numel(target_genes))
    for cnt = 1:numel(NotfoundEntrez)
        if isempty(NotfoundEntrez{cnt})
            fprintf('%s-%s\n',num2str(cnt),'N/A')
        else
            fprintf('%s-%s\n',num2str(cnt),NotfoundEntrez{cnt})
        end
    end
    fprintf('---------------------------------------------\n')
end

%Find Entrez ids
Entrez_targets = entrez_id(ismember(gene_symbol,target_genes));

% Find target genes in general metabolic model
Drugtargets = modelGenes(ismember(modelGenes,Entrez_targets));
if isempty(Drugtargets)
     fprintf('No target genes in general metabolic model for drug:%s\n',Drugname)
     fprintf('---------------------------------------------\n')
     Drugtargets = '';
     return
end
% check if some genes cannot be found in general metabolic model
NotfoundinGeneralModel = Entrez_targets(~ismember(Entrez_targets,modelGenes));
if ~isempty(NotfoundinGeneralModel)
    fprintf('Some genes are not in general metabolic model for drug:%s\n',Drugname)
    fprintf('Total targets for this drug:%d\n',numel(target_genes))
    convert2genesym = gene_symbol(ismember(entrez_id,NotfoundinGeneralModel));
    for cnt = 1:numel(convert2genesym)
        fprintf('%s-%s\n',num2str(cnt),convert2genesym{cnt})
    end
    fprintf('---------------------------------------------\n')
end
%==========================================================================