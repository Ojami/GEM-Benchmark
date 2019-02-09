function model = defineGrowthMedium(model,modelType,cellSpecificGrowth,cellCounter)
% General growth medium from modified RPMI-1640 of Folger et al.(2010)
% modelType: 1 for Recon1; 2 for Recon2
% cellSpecificGrowth: True is generating GEMs using Jain et al.(2012) data.

% Adjusting lb of exchange rxns: findEX_Rxns function of FASTCORMICS was
% employed here to identify all demand and exchange rxns;
if ~cellSpecificGrowth
    cellCounter = 1;
end
% inorganic_exrxns = {'EX_ca2','EX_cl','EX_fe2','EX_fe3','EX_h2o(','EX_h(',...
%     'EX_i(','EX_k(','EX_na1(','EX_o2(','EX_pi('}; % Only those in RPMI-1640
if modelType == 3
    inorganic_exrxns = {'HMR_9082','HMR_9150','HMR_9076','HMR_9047','HMR_9079','HMR_9081','HMR_9077','HMR_9048','HMR_9072'};
else
    inorganic_exrxns = {'EX_ca2','EX_cl','EX_fe2','EX_h2o(','EX_h(',...
        'EX_k(','EX_na1(','EX_o2(','EX_pi('}; % Only those in RPMI-1640
end
Temp_lb = [];
for ex1 = 1:numel(inorganic_exrxns) % Store current lbs for inorganic compounds in growth medium
    lbt = ~cellfun('isempty',strfind(model.rxns,inorganic_exrxns{ex1}));
    if sum(lbt) % Present in the model
        try
            Temp_lb(ex1) = model.lb(lbt);
            inorganic_exrxns(ex1) = model.rxns(lbt); % Build new exrxns strings based on input model
        catch % If more than one equivalent is found
            lbt = find(lbt);
            lbt = lbt(~cellfun('isempty',regexp(model.rxns(lbt),[inorganic_exrxns{ex1},'\W+'],'match')));
            Temp_lb(ex1) = model.lb(lbt);
            inorganic_exrxns(ex1) = model.rxns(lbt);
        end
    else
        Temp_lb(ex1) = 0;
    end
end

exRxns=findEX_Rxns(model); % A trimmed version of that in FASTCORMICS
exRxns = exRxns(strncmpi('ex_',exRxns,3)); % Not considering DM_ and Sink_ rxns
model.lb(ismember(model.rxns,exRxns)) = 0; % Set all exchange rxns' lb to 0

% Set lb of inorganic compounds in the growth medium to their original level
[n1,n2] = ismember(inorganic_exrxns,model.rxns);
model.lb(n2(n2>0)) = Temp_lb(n1); 

if cellSpecificGrowth % Use CORE data of Jain et al. (2012)
    OutRates = CORE2Uptakerates(cellCounter);
else
    OutRates = [-.5,-.05,-5,-.005,-.005,-.005,-.005,-.005,-.005,-.005,-.005,-.005,...
    -.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,...
    -.05,-.05,-.05,-.05,-.05,-.05]; 
end

if modelType == 1% Recon 1
    model = changeRxnBounds(model,{'EX_gln_L(e)','EX_gthrd(e)','EX_glc(e)',...
        'EX_thm(e)','EX_ribflv(e)','EX_pydx(e)','EX_ncam(e)','EX_inost(e)',...
        'EX_fol(e)','EX_pnto_R(e)','EX_chol(e)','EX_btn(e)',...
        'EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_Lcystin(e)',...
        'EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)',...
        'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_cys_L(e)',...
        'EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'}...
        ,OutRates,'l');
elseif modelType == 2 % Recon2
    load('exbound_Recon2')
    [~,Nidx1,Nidx2] = intersect(exbound_Recon2.r,model.rxns);
    model.lb(Nidx2) = exbound_Recon2.lb(Nidx1);
    model = changeRxnBounds(model,{'EX_gln_L(e)','EX_gthrd(e)','EX_glc(e)',...
        'EX_thm(e)','EX_ribflv(e)','EX_pydx(e)','EX_ncam(e)','EX_inost(e)',...
        'EX_fol(e)','EX_pnto_R(e)','EX_chol(e)','EX_btn(e)',...
        'EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_Lcystin(e)',...
        'EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)',...
        'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_cys_L(e)',...
        'EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'}...
        ,OutRates,'l');
elseif modelType == 3 % iHepatocyte
    OutRates(16) = []; % 'EX_Lcystin(e)' not in iHepatocyte
    model = changeRxnBounds(model,{'HMR_9063','HMR_9351','HMR_9034','HMR_9159',...
        'HMR_9143','HMR_9400','HMR_9378','HMR_9361','HMR_9146','HMR_9145',...
        'HMR_9083','HMR_9109','HMR_9066','HMR_9062','HMR_9070',...
        'HMR_9067','HMR_9038','HMR_9039','HMR_9040','HMR_9041','HMR_9042',...
        'HMR_9043','HMR_9065','HMR_9069','HMR_9044','HMR_9045','HMR_9064','HMR_9046'}...
        ,OutRates,'l');
elseif modelType == 4 % Recon3D
    load('exbound_Recon2')
    tempR = exbound_Recon2.r;
    tempB = exbound_Recon2.lb;
    tempR(1) = []; tempB(1) = [];
    
    tempR = strrep(tempR,'sink_','SK_');
    tempR = strrep(tempR,'(e)','_e');
    tempR = strrep(tempR,'(c)','_c');
    tempR = strrep(tempR,'(r)','_r');
    tempR = strrep(tempR,'_L','__L');
    tempR = strrep(tempR,'_R','__R');
    tempR(ismember(tempR,'EX_glc_e')) = {'EX_glc__D_e'};
    
%     %###### For Recon3D_2
%     tempR = strrep(tempR,'(e)','[e]');
%     tempR = strrep(tempR,'(c)','[c]');
%     tempR = strrep(tempR,'(r)','[r]');
%     tempR(ismember(tempR,'EX_glc_e')) = {'EX_glc_D[e]'};
    
    [~,Nidx1,Nidx2] = intersect(tempR,model.rxns);
    model.lb(Nidx2) = tempB(Nidx1);
else
    error('modelType input must be 1 or 2!')
end

% Set uptake of Alanine and Glutamate to zero (Refer to the supplementary
% materials)
model.lb(ismember(model.rxns,{'EX_ala_L[e]','EX_glu_L[e]'})) = 0;
model.lb(ismember(model.rxns,{'EX_ala_L(e)','EX_glu_L(e)'})) = 0;
model.lb(ismember(model.rxns,{'EX_ala__L_e','EX_glu__L_e'})) = 0;

% -------- Subfuctions ----------------------------------------------------
function [exRxns]=findEX_Rxns(model)
exRxnsInd=find(sum(abs(model.S),1)==1);

for i=1:numel(exRxnsInd)
    exmets= find(model.S(:,exRxnsInd(i)));
    if model.S(exmets,exRxnsInd(i))==1
        model.S(exmets,exRxnsInd(i))=-1;
    end
end
exRxns=model.rxns(exRxnsInd);

function OutRates = CORE2Uptakerates(CellNO)
% Growth medium uptake rates for RPMI-1640 based on the experimental data
% of Jain et al. 2012 (fmol/cell/h), Dolfie cell volum (pL) and typical 
% cell specific volume 4.3 mL/g [Ref. 10 of supplementary materials of
% Dolfi paper].
% Metabolite uptake rates (mmol/gDW/h) :
% C [fmol/cell/h]* 1e-12 [mmol/fmol] * 4.3 [mL/gDW] 
% -------------------------------------------------
%       cell volume [pL/cell] * 1e-9 [mL/pL]
%==========================================================================
% Read CORE data of Jain et al. 2012 --------------------------------------
% The coverage is 59/60: MDA-N and MDA-MB-468 are different.
load('defineGrowthMedium_needs')
% ---------------------------------------------------------------------------
c2mdata = {'glutamine','EX_gthrd(e)','glucose','thiamine','EX_ribflv(e)','EX_pydx(e)',...
    'EX_ncam(e)','EX_inost(e)','folate','pantothenate','choline','biotin',...
    'arginine','asparagine','aspartate','EX_Lcystin(e)','glycine','EX_his_L(e)',...
    'isoleucine','leucine','lysine','methionine','phenylalanine','EX_cys_L(e)',...
    'serine','threonine','tryptophan','tyrosine','valine'};
fill_missed = [-.05,-.005,-.005,-.005,-.005,-.05,-.05,-.05];
CellNo45 = [-.5,-.05,-5,-.005,-.005,-.005,-.005,-.005,-.005,-.005,-.005,-.005,...
    -.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,-.05,...
    -.05,-.05,-.05,-.05,-.05,-.05];
allmet = ({});
ncitemp = nci;
ncitemp(45) = [];
ctr = 1;
for cm = 1:numel(c2mdata)
    wheremets = find(strcmpi(mets,c2mdata{cm}));
    if ~isempty(wheremets)
        met = dat(wheremets,:); % for metabolite of interest
        met = met./vols; % normalize to cell volume
        % Convert to uptake rates:
        met = met.*1e-3.*4.3;
    else
        met = ones(1,59).*fill_missed(ctr);
        ctr = ctr + 1;
    end
    [~,n2]=ismember(ncitemp,corecell);
    n2(~n2) = [];
    met = met(n2); % Rearrange based on NCI-60 cells
    met = [met(1:44),CellNo45(cm),met(45:end)];
    allmet{cm} = met';
end
OutRates = (0);
if CellNO < 60
    for cm = 1:numel(c2mdata) % Extract the target cell line's uptake rates
        OutRates(cm) = allmet{cm}(CellNO);
    end
else % Mean of all 60 cell lines
    OutRates = [-0.178839189825914,-0.0500000000000000,-0.782745974685274,-1.88592049920245e-05,-0.00500000000000000,-0.00500000000000000,-0.00500000000000000,-0.00500000000000000,-7.76817209040688e-05,-9.46567803166800e-05,-0.00176700412403146,-8.33334625037071e-05,-0.0122648799699342,-0.00520232253017099,-0.00362304517326634,-0.0500000000000000,0.000508189295412696,-0.0500000000000000,-0.0128776702605091,-0.0156312273625244,-0.0160981535539064,-0.00548236814552869,-0.00597847578195898,-0.0500000000000000,-0.0256512691277860,-0.0107048352258026,-0.00218812025596096,-0.00778873365523522,-0.0116134859981881];
end
OutRates(OutRates>0) = 0;

