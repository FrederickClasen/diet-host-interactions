%% ADD PATHS

addpath(genpath('/Users/clasenf/OneDrive - The Francis Crick Institute/METEOR'));      % path to directory
addpath(genpath('/Users/clasenf/GEM/Raven/'));                                     % RAVEN toolbox
addpath(genpath('/Users/clasenf/GEM/COBRA/'));                                     % COBRA toolbox
changeCobraSolver('mosek','all');
setRavenSolver('mosek');

%% LOAD AGORA MODELS

experiment = "nonDEN_Fed/"; %;"nonDEN_Fasted/";"DEN_Fed/";"DEN_Fasted/"};
PATH = "data/models/selected/"+experiment;
files = dir(PATH);
models = struct();
for j=3:length(files)
    MODELNAMEPATH = files(j).name;
    f = PATH + MODELNAMEPATH;
    model = load(f,'model');
    model = model.model;
    MODELNAME = erase(MODELNAMEPATH,'.mat');
    parts = split(MODELNAME,'_');
    shortName = lower(strcat(parts{1}(1),parts{2}(1)));
    shortName = shortName + string(j);
    %assignin('base',MODELNAME,model)
    models(j-2).modelname = MODELNAME;
    models(j-2).shortname = char(shortName);
    models(j-2).model = model;
    %models(i-2).model.equations = constructEquations(models(i-2).model);
end

%% MAKE THE COMMUNITY MODEL 

for i=1:length(models)
    disp(i);
    tModel = models(i).model;
    
    % rename cytosolic compartment
    newComp = "["+char(models(i).shortname)+"_c]";
    tModel.mets = strrep(tModel.mets,'[c]',newComp);
    
    % remove exchange reactions - this removes the biomass exchange as well
    % we add a new biomass to the [e] compartment
    exchRxns = tModel.rxns(contains(tModel.rxns,'Biomass'));
    %exchRxns = getExchangeRxns(tModel);
    models(i).exchRxns = exchRxns;
    tModel = removeReactions(tModel,exchRxns,true,true);
    
    
    % add several model fields for merging
    tModel.id = models(i).modelname;
    tModel.rev = zeros(length(tModel.rxns),1);
    tModel.rev(find(tModel.lb == -1000)) = 1;
    tModel.comps = {char(models(i).shortname+"_c");'e'};
    tModel.compNames = {char(models(i).shortname+"_c");'e'};
    tModel.metComps = zeros(length(tModel.mets),1);
    for j=1:length(tModel.mets)
       if contains(tModel.mets(j),models(i).shortname)
           tModel.metComps(j) = 1;
       else  
           tModel.metComps(j) = 2;
       end
    end
    
    % add biomass to the [e] compartment
    %biomassMet = string(tModel.mets(find(contains(tModel.mets,'biomass'))));
    newBiomassEq = "Produced biomass["+char(models(i).shortname)+"_c]"+" => biomass[e]";
    rxnsToAdd.rxns = {char("Biomass"+char(models(i).shortname))};
    rxnsToAdd.equations = {char(newBiomassEq)};
    rxnsToAdd.c = 1;
    tModel = addRxns(tModel,rxnsToAdd,3,[],true);
    tModel.equations = constructEquations(tModel);
    
    tModel.rxns = cellstr(char(models(i).shortname)+"_" + tModel.rxns);
    % merge to mergedModel
    if i==1
        mergedModel = tModel;
    else
        mergedModel = mergeTwoModels(mergedModel,tModel);
    end
end

%% TEST IF ALL MODELS CAN DO BIOMASS

for i=1:length(models)
    tModel = models(i).model;
    biomassRxn = tModel.rxns(contains(tModel.rxns,'EX_biomass'));
    tModel = setParam(tModel,'obj',biomassRxn,1);
    fba = optimizeCbModel(tModel);
    disp(fba.f);
end 

%% GET ALL EXCHANGE METABOLITES

for i=1:length(models)
    I = getIndexes(models(i).model,findMetFromCompartment(models(i).model,'[e]'),'mets');
    eMets = models(i).model.metNames(I);
    models(i).exchangeMets = eMets;
end

allExchange = models(1).exchangeMets;
for i=1:length(models)
    allExchange = union(allExchange,models(i).exchangeMets);
end 

%% SIMULATE FLUX FOR ALL COMMUNITY MODELS

load('mergedModels.mat');
speciesKO = struct();
diets = {'WD';'CD';'WDCARBS';'WDLIPIDS';'WDAA';'WDAACARBS';'WDAALIPIDS';'WDCARBSLIPIDS'};
counter = 0;
for j=1:length(diets)
    diet = string(diets{j});
    for i=1:length(mergedModels)
        counter = counter + 1;
        
        disp('Fixing metComps ...');
        newModel = fixMetComps(mergedModels(i).model);
        
        disp('Formatting exchange reactions ...');
        newModel = formatExchangeRxns(newModel);
        
        disp('Adding diet constraint ...');
        newModel = dietConstraint(newModel,diet);
        
        disp('Adding biomass constraint ...');
        newModel = biomassConstraint(newModel,mergedModels(i).modelname,mergedModels(i).modelNames,mergedModels(i).shortnames);
        
        % this part for the species KO
        
%         speciesKO(counter).diet = diet;
%         speciesKO(counter).modelname = mergedModels(i).modelname;
%         shortnames = mergedModels(i).shortnames;
%         modelnames = mergedModels(i).modelNames;
%         modelKO = struct();
%         for k=1:length(shortnames)
%             shortname = shortnames{k};
%             rxns = cellstr(newModel.rxns(contains(newModel.rxns,shortname)));
%             tModel = setParam(newModel,'eq',rxns,0);
%             I = getIndexes(tModel,'EX_glyc(e)','rxns');
%             try
%                 fba = optimizeCbModel(tModel);
%                 modelKO(k).species = modelnames{k};
%                 %modelKO(k).rxns = tModel.rxns;
%                 %modelKO(k).equations = tModel.equations;
%                 if fba.stat
%                     modelKO(k).flux = fba.x;
%                     disp(fba.x(I))
%                 else
%                     modelKO(k).flux = zeros(length(tModel.rxns),1);
%                 end
%             catch
%                 modelKO(k).species = modelnames{k};
%                 %modelKO(k).rxns = tModel.rxns;
%                 %modelKO(k).equations = tModel.equations;
%                 modelKO(k).flux = zeros(length(tModel.rxns),1);
%                 disp('ERROR');
%             end
%         end
%         speciesKO(counter).KO = modelKO;
        
        fba=optimizeCbModel(newModel);
        %I = getIndexes(newModel,'EX_glyc(e)','rxns');
        %disp(fba.x(I));
        tempModel = newModel;
        T = table(tempModel.rxns,...
                  tempModel.equations);
        T2 = array2table(fba.x);
        T2.Properties.VariableNames = {'FLUX'};
        T = [T T2];
        colheads = {'ID','EQUATION','FLUX'};
        T.Properties.VariableNames = colheads;
        writetable(T,'data/models/Fluxes/Fluxes'+diet+'.xlsx','Sheet',mergedModels(i).modelname,'WriteVariableNames',true);
    end
end







