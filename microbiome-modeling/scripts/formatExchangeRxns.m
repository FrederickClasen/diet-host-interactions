function newModel=formatExchangeRxns(model)
mModel = model;
% rename the exchange reactions without prefix
biomassExchange = mModel.rxns(contains(mModel.rxns,'EX_biomass'));
mModel = removeReactions(mModel,biomassExchange);
for j=1:length(mModel.rxns)
    if contains(mModel.rxns{j},'EX_')
        [~,remain] = strtok(mModel.rxns{j},'_');
        remain = strrep(remain,'_EX','EX');
        mModel.rxns{j} = remain;
    end
end

mModel.rxns = cellstr(mModel.rxns);

% add community biomass reaction and set objective for all biomass rxns and
% exchange biomass
rxnsToAdd.rxns = {'community_biomass'};
rxnsToAdd.equations = {'biomass[e] =>'};
rxnsToAdd.c = 1;
mModel = addRxns(mModel,rxnsToAdd,3);
biomassObjectives = mModel.rxns(contains(mModel.rxns,'biomass'));
mModel = setParam(mModel,'obj',biomassObjectives,1);
mModel.equations{length(mModel.rxns)} = 'biomass[e] =>';
%fba = optimizeCbModel(mModel);
%disp(fba.f);
newModel=mModel;
end