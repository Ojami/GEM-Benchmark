function [ParsedGPR,corrRxn] = GPRparser(model)
% Modified version of subfunction extractGPRs of COBRA toolbox, with a much
% less CPU-time

AllRules = model.grRules;
[AllGPR,corrRxn] = deal({});
for i = 1:numel(AllRules)
    if isempty(AllRules{i})
        continue
    end
    EachRule = AllRules{i};
    orFind = textscan(EachRule,'%s','delimiter','or'); orFind = orFind{1};
    orFind(cellfun('isempty',orFind)) = [];

    GPR_temp = regexp(orFind,'(\d*)[.](\d*)|(\d*)|(\w+)(\d*)[.](\d*)|(\w+)(\d+)','match');
    GPR_temp(cellfun('isempty',GPR_temp)) =[];
    corrRxn = [corrRxn;repmat(model.rxns(i),numel(GPR_temp),1)];
    AllGPR = [AllGPR;GPR_temp];
end

maxSize = max(cellfun(@numel,AllGPR));  
FillerFunc = @(x) [x repmat({'Del'},1,maxSize-numel(x))]; 
ParsedGPR = cellfun(FillerFunc,AllGPR,'UniformOutput',false); 
ParsedGPR = vertcat(ParsedGPR{:});
ParsedGPR = strrep(ParsedGPR,'Del','');
