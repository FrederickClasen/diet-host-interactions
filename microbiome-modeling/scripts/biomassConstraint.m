function newModel=biomassConstraint(model,experiment,modelnames,modelshortnames)

fModel = model;
objectives = cellstr(fModel.rxns(contains(fModel.rxns,'Biomass')));
objectives{length(objectives) + 1} = 'community_biomass';
fModel = changeObjective(fModel,objectives,1);
% change the constraints for biomass to relative abundances
biomass_constraints = importdata("data/models/selected/"+strrep(experiment,'/','')+"_abundances.csv");
biomass_constraints.textdata(1,:) = [];
biomass_constraints.textdata(:,2) = [];
biomass_constraints.textdata(:,2) = [];
biomass_constraints.model = biomass_constraints.textdata;

for i=1:length(biomass_constraints.model)
     if modelnames{i} == biomass_constraints.model{i}
         bound = biomass_constraints.data(i,2);
         rxn = modelshortnames{i}+"_Biomass"+modelshortnames{i};
         fModel = setParam(fModel,'lb',rxn,bound);
         fModel = setParam(fModel,'ub',rxn,1000);
         %fModel = setParam(fModel,'ub',rxn,bound+(1000*bound));
     end
end

%fba = optimizeCbModel(fModel);
newModel = fModel;
% code to open exchange reactions one by one
% exchRxns = cellstr(fModel.rxns(contains(fModel.rxns,'EX_')));
% sModel = fModel;
% for i=1:length(exchRxns)
%     sModel = setParam(sModel,'lb',exchRxns(i),-1000);
%     fba = optimizeCbModel(sModel);
%     disp(i);
%     if abs(fba.f) > 0
%         disp(exchRxns(i))
%         break
%     end
% end
end