function ResPower_Noisy(algorithmNames)
% Similarity level for GEMs generated for Noisy data
% Oveis Jamialahmadi
% Jan 2018

clc
load Noisyset % Set of noisy expression data generated for Cell line #2 in current panel

folderNames = strcat(algorithmNames,'_Noisy');
algorithmNames = strcat(algorithmNames,'Noisy');
fprintf('Resolution power for algorithmNames....\n')
fprintf('------------------------------------------------\n')
fprintf('Following algorithm are going to be tested:\n')
for i = 1:numel(algorithmNames)
    fprintf('%d-%s\n',i,algorithmNames{i})
end
fprintf('------------------------------------------------\n')
changeCobraSolver('gurobi5');

algorithmNames = sort(algorithmNames);
folderNames = sort(folderNames);

%--------------------------------------------------------------------------
for count = 1:numel(algorithmNames)
    fprintf('%d-%s (of %d)\n',count,algorithmNames{count},numel(algorithmNames))

    numModels = 1:numel(dir([pwd,'\GEMs_Noisy\',folderNames{count},'\*.mat']));
    GEMnames = strcat(algorithmNames{count},num2str(numModels'),'.mat');
    GEMnames = cellstr(GEMnames);
    GEMnames = strrep(GEMnames,' ','');
    Pairwise_mat = zeros(numel(GEMnames),numel(GEMnames));
    for ct1 = 1:numel(GEMnames)
        load(GEMnames{ct1})
        model1 = OutM; clear OutM
        for ct2 = 1:numel(GEMnames)
            load(GEMnames{ct2})
            model2 = OutM; clear OutM
            % Jacard index
            try
               Pairwise_mat(ct1,ct2) = ...
                sum(ismember(model1.rxns,model2.rxns))/...
                (numel(model1.rxns)+numel(model2.rxns)-...
                sum(ismember(model1.rxns,model2.rxns)));
            catch
                Pairwise_mat(ct1,ct2) = NaN;
            end
            clear model2
        end
        clear model1
    end
%     Pairwise_mat(isnan(Pairwise_mat)) = 0;
    [~,N2] = find(isnan(Pairwise_mat));
    N1 = (0); cnt = 1;
    while numel(N2) > size(Pairwise_mat,1)
        [~,N1(cnt)] = max(histc(N2,1:size(Pairwise_mat,1)));
        N2(N2==N1(cnt)) = [];
        cnt = cnt + 1;
    end
    Pairwise_mat(N1,:) = [];
    Pairwise_mat(:,N1) = [];
    hfig = figure;
    colormap(hfig,'jet')
    imagesc(Pairwise_mat);
    colorbar
    R = rho; % From noisy data
    R(N1) = [];
    GEMnames(N1) = [];
    set(gca,'YTick',1:2:numel(GEMnames))
    set(gca,'XTick',1:2:numel(GEMnames))
    R = strtrim(cellstr(num2str(R')));
    newR = regexp(R,'0.(\d{2})','match');
    R_empty = logical(cellfun('isempty',newR));
    R(~R_empty) = [newR{:}];
    set(gca,'YTickLabel',R(1:2:end))
    set(gca,'XTickLabel',R(1:2:end))
	set(gcf,'name',algorithmNames{count},'numbertitle','off')
    title(algorithmNames{count})
    set(gca, 'Ticklength', [0 0],'FontSize',17,'FontName','Cambria')
    try
        set(gca,'XTickLabelRotation',90)
    catch
    end
    xlabel(gca,'Spearman R','FontSize',19,'FontWeight','Bold','FontName','Cambria')
    ylabel(gca,'Spearman R','FontSize',19,'FontWeight','Bold','FontName','Cambria')
    
    %# Disabled post-processing steps
%     save(algorithmNames{count},'R','Pairwise_mat')
%     savefig(hfig,algorithmNames{count})
%     print(hfig,'-r500','-dpng',algorithmNames{count})
%     close(gcf)
%     clear hfig
    
end
