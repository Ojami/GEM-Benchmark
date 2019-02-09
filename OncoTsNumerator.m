function OncoTsNumerator(algorithmNames,modelType)
% Counts the number of TS, Oncogenes and Loss of functions (LOF) in GEMs
% NOTE: Hypergeometric test was set to 'lower cumulative' for LOF and TS datasets,
% because the lower number of TS/LOFs in the output GEM is of interest as
% opposed to the oncogenes.

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

fprintf('Nmuerating TS, Oncogenes and LOFs....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')

alSize = (0);

% load required datasets
load TS_Onco_LOF_list

if modelType == 1
    load Recon1activeGenes
    Recon_genes = Recon1activeGenes; clear Recon1activeGenes
elseif modelType == 2
    load Recon2activeGenes
    Recon_genes = Recon2activeGenes; clear Recon2activeGenes
else
    error('Wrong modelType!')
end
[numOnco,numTS,numLOF] = deal({});
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
    outRes = struct('onco',0,'ts',0,'lof',0,'active',0,'recon',0,'on1',0,'ts1',0,'lof1',0,...
            'ponco',0,'pts',0,'plof',0);
    for ct = 1:alSize(count)
        load(GEMnames{ct})
        curgenes = OutM.genes;
        curgmat = OutM.rxnGeneMat;
        % Identify active genes in each model
        ag = 1; activeGenes = ({});
        for i1 = 1:numel(curgenes) 
            if ~isempty(find(curgmat(:,i1), 1))
                activeGenes(ag) = curgenes(i1);
                ag = ag + 1;
            end
        end
        activeGenes = unique(activeGenes);
        oncoCount = intersect(activeGenes,Oncogenes);
        TSCount = intersect(activeGenes,TumSupp);
        LOFCount = intersect(activeGenes,LOFset);
        
        outRes(ct).onco = oncoCount;
        outRes(ct).ts = TSCount; 
        outRes(ct).lof = LOFCount;
        outRes(ct).active = activeGenes;
        outRes(ct).recon = numel(Recon_genes);
        outRes(ct).on1 = numel(intersect(Recon_genes,Oncogenes));
        outRes(ct).ts1 = numel(intersect(Recon_genes,TumSupp));
        outRes(ct).lof1 = numel(intersect(Recon_genes,LOFset));
        outRes(ct).ponco = 1-hygecdf(numel(oncoCount)-1,numel(Recon_genes),outRes(ct).on1,numel(activeGenes)); % Upper cumulative Q
        
        outRes(ct).pts = hygecdf(numel(TSCount),numel(Recon_genes),outRes(ct).ts1,numel(activeGenes)); % Lower cumulative P
        outRes(ct).plof = hygecdf(numel(LOFCount),numel(Recon_genes),outRes(ct).lof1,numel(activeGenes)); % Lower cumulative P
        
        clear oncoCount TSCount LOFCount OutM
    end
    [numOnco{count},numTS{count},numLOF{count}] = NumAnalzer(outRes);
end

% Dsip formatted results
oncoC = zeros(1,numel(numOnco)); TSC = oncoC; LOFC = oncoC;
for i = 1:numel(numOnco)
    if numOnco{i}.en
        oncoC(i) = 1;
    end
    if numTS{i}.en
        TSC(i) = 1;
    end
    if numLOF{i}.en
        LOFC(i) = 1;
    end
end
% dispResults(oncoC,numOnco,algorithmNames,alSize)
if sum(oncoC)
    fprintf('--------- Oncogenes --------------------\n')
    dispResults(oncoC,numOnco,algorithmNames,alSize)
end
if sum(TSC)
    fprintf('--------- Tumor suppressors --------------------\n')
    dispResults(TSC,numTS,algorithmNames,alSize)
end
if sum(LOFC)
    fprintf('--------- Loss of functions --------------------\n')
    dispResults(LOFC,numLOF,algorithmNames,alSize)
end

% ====================subfunctions=========================================
function dispResults(oncoC,numOnco,algorithmNames,alSize)
 numOnco = squeeze(cell2mat(struct2cell([numOnco{logical(oncoC)}])))';
 algorithmNames = algorithmNames(logical(oncoC));
 alSize = alSize(logical(oncoC));
 
 disp(table(numOnco(:,1),numOnco(:,2),numOnco(:,3),numOnco(:,4)./alSize'...
     ,'VariableNames',{'En','pval','Enstd','f'},'RowNames',algorithmNames'))


function [reso,rest,resl]= NumAnalzer(outRes)

co = 1; ct = 1; cl = 1; eno = (0); ent = (0); enl = (0);
pvalo = (0); pvalt = (0); pvall = (0);
for count = 1:numel(outRes)
    if outRes(count).ponco<=0.05
        eno(co) = (numel(outRes(count).onco)/numel(outRes(count).active))/(outRes(count).on1/outRes(count).recon);
        pvalo(co) = outRes(count).ponco;
        co = co + 1;
    end
    if outRes(count).pts<=0.05
        ent(ct) = (numel(outRes(count).ts)/numel(outRes(count).active))/(outRes(count).ts1/outRes(count).recon);
        pvalt(ct) = outRes(count).pts;
        ct = ct + 1;
    end
    if outRes(count).plof<=0.05
        enl(cl) = (numel(outRes(count).lof)/numel(outRes(count).active))/(outRes(count).lof1/outRes(count).recon);
        pvall(cl) = outRes(count).plof;
        cl = cl + 1;
    end
end

reso.en = mean(eno);
reso.p = mean(pvalo);
reso.enstd = std(eno);
reso.num = co - 1;

rest.en = mean(ent);
rest.p = mean(pvalt);
rest.enstd = std(ent);
rest.num = ct - 1;

resl.en = mean(enl);
resl.p = mean(pvall);
resl.enstd = std(enl);
resl.num = cl - 1;
