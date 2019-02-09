function MetabolicEssentiality_sp(PotenGenes,algorithmName,modelType)
% Identifying essential metabolic genes in generated models and compare
% them with gene dependency scores (<0) of CERES (Meyers et al. (2017)).

% Oveis Jamialahmadi
% Jan 2018

PotenGenes_main = PotenGenes;
if modelType == 1
    load Recon1_genes
elseif modelType == 2
    load Recon2_genes
else
    error('Wrong modelType!')
end


%#### CELL LINE BASED
% sp_cat = {'SKOV3','HOP62','MALME3M','NCIH322','HS578T','MDAMB435S','EKVX','HT29','OVCAR5','KM12','COLO201','KPL1','CAKI1','UACC257','SW620','UACC62','MDAMB231','A549','OVCAR8','SF295','U251MG','BT549','786O','NCIH460','SR786','RPMI8226'};
sp_cat = {'SKOV3','HOP62','MALME3M','HS578T','EKVX','HT29','OVCAR5','KM12','CAKI1',...
    'UACC257','SW620','UACC62','MDAMB231','A549','OVCAR8','SF295','U251MG','BT549',...
    '786O','NCIH460','SR786','RPMI8226'};
can_cat = sp_cat;
sp_cat = strcat('EGlist_',sp_cat);
% Approximate cells #### FOR FUTURE DEVELOPMENT PURPOSES ##################
% id_cat_60 = [3,5,7,9,13,14,16,18,20,21,25,27,30,31,33,36,37,41,42,43,44,50,54,55,56,58];
% id_cat_59 = [3,5,7,9,13,14,16,18,20,21,25,27,30,31,33,36,37,41,42,43,44,50,53,54,55,57];
% #########################################################################
% Exact cells: indices correspond to NCI60_HGU133A.mat
id_cat_60 = [3,5,7,13,16,18,20,21,30,31,33,36,37,41,42,43,44,50,54,55,56,58];
id_cat_59 = [3,5,7,13,16,18,20,21,30,31,33,36,37,41,42,43,44,50,53,54,55,57];


for sp = 1:numel(sp_cat)
    
    fprintf('Cell line: %s\n',can_cat{sp})
    load([sp_cat{sp},'.mat'])
    eval(['EGlist_general = ',sp_cat{sp},';'])    
    eval(['clear ',sp_cat{sp}])

    PotenGenes = PotenGenes_main;
    if size(PotenGenes,1) < 60
        PotenGenes = PotenGenes(id_cat_59(sp),:);
    else
        PotenGenes = PotenGenes(id_cat_60(sp),:);
    end
    %##################### FOR FUTURE DEVELOPMENT PURPOSES ################
%   ------------------- CORRELATION BASED --------------------------------
%     genelist = PotenGenes{1,4};
%     grRatio = PotenGenes{1,3};
% %     genelist(grRatio > (1-1e-6)) = []; grRatio(grRatio > (1-1e-6)) = []; % Remove ones
%     % #Modify scores
%     EG_scores(EG_scores>=-.5) = 1; EG_scores(EG_scores<-.5) = 0; 
%     grRatio(grRatio<=0.01) = 0; grRatio(grRatio>0.01) = 1;
%     [~,n1,n2] = intersect(EG_genes,genelist,'stable');
%     EG_genes = EG_genes(n1); EG_scores = EG_scores(n1);
%     genelist = genelist(n2); grRatio = grRatio(n2);
%     if size(grRatio,1) == 1
%         grRatio = grRatio';
%     end
%     clear n1 n2
%     if iscell(EG_scores)
%         EG_scores = [EG_scores{:}]';
%     end
%     PearR = (0); PearP = (1);
%     [PearR,PearP] = corr(grRatio,EG_scores,'Type','Spearman');

%     %####################################################################

    Allfound = ({}); 
    for ct = 1:size(PotenGenes,1)

        genelist = PotenGenes{ct,1}; % for each model
        genelist = unique(genelist);
        if isnumeric(genelist)
            genelist = {'XXXX'};
        end
        genelist1 = unique(Recon_genes); % All unique genes in background GEM

        found_cells = EGlist_general(ismember(EGlist_general,genelist));
        found_cells1 = EGlist_general(ismember(EGlist_general,genelist1));
        Allfound{ct}.frac1 = numel(found_cells)/numel(genelist);
        Allfound{ct}.ess1 = numel(found_cells);
        Allfound{ct}.gene1 = numel(genelist);
        Allfound{ct}.ess2 = numel(found_cells1);
        Allfound{ct}.gene2 = numel(genelist1);
        Allfound{ct}.frac2 = numel(found_cells1)/numel(genelist1); % Generic metabolic model
        Allfound{ct}.enrich = Allfound{ct}.frac1/ Allfound{ct}.frac2;
        % Compute hypergeometric test
        Allfound{ct}.p= hygecdf(numel(found_cells)-1,numel(genelist1),numel(found_cells1),numel(genelist),'upper');

    end
    CERESraw = Allfound;
    clear Allfound

    % Analyze and display results
    outData = MetEssAnalysis(CERESraw,algorithmName);
    fprintf('====================================\n')
end

fprintf('===============DONE=========================\n')
end
% ====================subfunctions=========================================
function CERESres = MetEssAnalysis(CERESraw,algorithmName)


CERESres = struct('frac',0,'p',1,'pstd',1,'en',0,'enstd',0);

    
c = 1; this_p = (1); this_frac = (0); 
for count = 1:numel(CERESraw)
    Allfound = CERESraw{count};
    if Allfound.p < 0.05
        this_p(c) = Allfound.p;
        this_frac(c) = Allfound.enrich;
        c=c+1;
    end
    clear Allfound
end
CERESres.frac = (c-1)/numel(CERESraw);
CERESres.p = mean(this_p);
CERESres.pstd = std(this_p);
CERESres.en = mean(this_frac);
CERESres.enstd =std(this_frac);

CERESres = squeeze(cell2mat(struct2cell(CERESres)))';

disp(table(CERESres(2),CERESres(4)...
 ,'VariableNames',{'pval','Enrichment'},'RowNames',algorithmName'))
end