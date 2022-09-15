%% SIMULATION

clear;
experiment = 'SIM';
disp('----- loading multi-tissue model -----');
load('data/GEM/multiTissueModels.mat');
load('data/diet.mat');

gut_sinks = {'EX_nmn(e)';'EX_pnto_R(e)';'EX_4hbz(e)';'EX_gthox(e)';...
             'EX_12dgr180(e)';'EX_spmd(e)';'EX_sheme(e)';'EX_q8(e)';...
             'EX_pheme(e)';'EX_cgly(e)';'EX_2dmmq8(e)';'EX_26dap_M(e)';...
             'EX_zn2(e)';'EX_so4(e)';'EX_mg2(e)';'EX_mn2(e)';'EX_k(e)';'EX_cu2(e)';...
             'EX_cobalt2(e)';'EX_ca2(e)';'EX_cl(e)';'EX_thm(e)';'EX_fe2(e)';'EX_fe3(e)';...
             'EX_etoh(e)';'EX_fol(e)';'EX_pi(e)';'EX_ribflv(e)';...
             'EX_na1(e)';'EX_thm(e)';'EX_adocbl(e)';'EX_h2o(e)';'EX_zn2(e)';'EX_pydx(e)';...
             'EX_chol(e)';'EX_mqn7(e)';'EX_mqn8(e)';'EX_xan(e)';'EX_Lcystin(e)';'EX_4hpro(e)';...
             'EX_cys_L(e)';'EX_h(e)';'EX_cl(e)';'EX_adn(e)';'EX_amp(e)';'EX_dad_2(e)';...
             'EX_dcyt(e)';'EX_h(e)';'EX_hxan(e)';'EX_ins(e)';'EX_thymd(e)';'EX_ura(e)';...
             'EX_uri(e)';'EX_lac_L(e)';'EX_lac_D(e)';'EX_rib_D(e)';'EX_arab_L(e)';'EX_no3(e)';...
             'EX_so4(e)';'EX_but(e)'};
         
blood_sinks = {'EXC_BOTH_C00301[b]';'EXC_BOTH_C00120[b]';'EXC_BOTH_C00016[b]';...
               'EXC_BOTH_C14818[b]';'EXC_BOTH_C14819[b]';'EXC_BOTH_C00001[b]';'EXC_BOTH_C00318[b]';...
               'EXC_BOTH_C00003[b]';'EXC_BOTH_C00006[b]';'EXC_BOTH_C00255[b]';'EXC_BOTH_C00101[b]';...
               'EXC_BOTH_C00214[b]';'EXC_BOTH_Cx1003[b]'};
                      
diets = {'WD','CD'};
           
for i=1:length(models)
    
    disp(models(i).modelName);
    
    disp('----- constraining model -----');
    
    model = models(i).model;
    
    % set blood and gut uptake to zero
    blood_exchange = model.rxns(find(contains(string(model.subSystems),"Blood exchange")));
    model = setParam(model,'lb',blood_exchange,0);
    gut_exchange = model.rxns(find(contains(string(model.subSystems),"Exchange/demand reaction")));
    model = setParam(model,'lb',gut_exchange,0);
    
    % OBJECTIVE CONSTRAINTS
    % biomass objectives for microbiome
    microbiome_objectives = model.rxns(find(contains(string(model.rxns),"biomass")));
    hmr = not( cellfun( @isempty, regexp( microbiome_objectives, 'HMR' ) ) );
    microbiome_objectives( hmr ) = [];
    exc = not( cellfun( @isempty, regexp( microbiome_objectives, 'EXC' ) ) );
    microbiome_objectives( exc ) = [];
    %tissue objectives
    %tissue_objectives = model.rxns(find(contains(string(model.rxns),"MMRN_Biomass")));
    tissue_objectives = {'t_MMRN_Biomass'};

    objectives = [microbiome_objectives;tissue_objectives];
    model = setParam(model,'obj',objectives,1);
    
    % add sink reactions from predefined lists (above)
    sink_reactions = model.rxns(find(contains(string(model.rxns),"sink")));
    sinks = [sink_reactions;gut_sinks;blood_sinks];
    if find(contains(model.rxns,'EX_arab_D(e)'))
        sinks{end+1} = 'EX_arab_D(e)';
    end 
    model = setParam(model,'lb',sinks,-0.5);
    
    % EFLUX CONSTRAINTS
    % kidney and wat
    model = constrainRxns(model,'data/Eflux/Kidney/'+models(i).modelName+'.csv',2.5);
    model = constrainRxns(model,'data/Eflux/WAT/'+models(i).modelName+'.csv',2.5);
    
    % liver and tumour
    if contains(models(i).modelName,'nonDEN')
        model = constrainRxns(model,'data/Eflux/Liver/'+models(i).modelName+'.csv',2.5);
        tumour_rxns = model.rxns(find(contains(string(model.subSystems),"t_")));
        gut_to_blood = not( cellfun( @isempty, regexp( tumour_rxns, 'EX_M2B' ) ) );
        tumour_rxns( gut_to_blood ) = [];
        model = setParam(model,'eq',tumour_rxns,0);
        model = setParam(model,'eq',model.rxns(find(contains(string(model.subSystems),"l2t_Exchange reactions"))),0);
    else
        model = constrainRxns(model,'data/Eflux/Liver/'+models(i).modelName+'.csv',2.5);
        model = constrainRxns(model,'data/Eflux/Tumour/'+models(i).modelName+'.csv',2.5);
    end
    
    % RER CONSTRAINT
    if contains(models(i).modelName,'nonDEN')
        model = setParam(model,'ub','EXC_BOTH_C00007[b]',-15);
        model = setParam(model,'lb','EXC_BOTH_C00007[b]',-1000);
        model = setParam(model,'ub','EXC_BOTH_C00011[b]',1000);
        model = setParam(model,'lb','EXC_BOTH_C00011[b]',11);
    else
        model = setParam(model,'ub','EXC_BOTH_C00007[b]',-17);
        model = setParam(model,'lb','EXC_BOTH_C00007[b]',-1000);
        model = setParam(model,'ub','EXC_BOTH_C00011[b]',1000);
        model = setParam(model,'lb','EXC_BOTH_C00011[b]',13);
    end
    
    %CDP-DAG is the only metabolite needed from LDL so can add this or LDL
    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'EXC_BOTH_m46[b]',...
                      'l_EXC_BOTH_m46[c]',...
                      't_EXC_BOTH_m46[c]',...
                      'w_EXC_BOTH_m46[c]'};
    rxnsToAdd.equations = {'CDP-diacylglycerol-LD-PI pool[b] <=> ',...
                           'CDP-diacylglycerol-LD-PI pool[l_c] <=> CDP-diacylglycerol-LD-PI pool[b]',...
                           'CDP-diacylglycerol-LD-PI pool[t_c] <=> CDP-diacylglycerol-LD-PI pool[b]',...
                           'CDP-diacylglycerol-LD-PI pool[w_c] <=> CDP-diacylglycerol-LD-PI pool[b]'};
    rxnsToAdd.subSystems = {'Blood exchange',...
                           'l_Exchange reactions',...
                           't_Exchange reactions',...
                           'w_Exchange reactions'};
    model = addRxns(model,rxnsToAdd,3,[],true);
    model = setParam(model,'lb','EXC_BOTH_m46[b]',-0.5);
    
    % remove any tissue here
    kidney = cellstr(model.rxns(find(contains(string(model.subSystems),"k_"))));
    wat = cellstr(model.rxns(find(contains(string(model.subSystems),"w_"))));
    rem = [kidney;wat];
    model = removeRxns(model,rem);
    
    %model = removeReactions(model,{'w_HMR_5233','w_HMR_5238'},true,true,true);
    
    dna = cellstr(model.mets(find(contains(model.metNames,"DNA replication"))));
    protein = cellstr(model.mets(find(contains(model.metNames,"Protein biosynthesis"))));
    rna = cellstr(model.mets(find(contains(model.metNames,"RNA transcription"))));
    pep = cellstr(model.mets(find(contains(model.metNames,"Peptidoglycan polymer (n-1 subunits)"))));
    chol = cellstr(model.mets(find(contains(model.metNames,"cholesterol-ester pool"))));
    rem = [dna;protein;rna;pep;chol];
    model = removeMets(model,rem,false,true,true,true);
    model.equations = constructEquations(model);
    
    for j=1:length(diets)
        % DIET CONSTRAINT
        dmodel = setParam(model,'lb',diet.aa.rxns,-diet.aa.ub(:,j));
        dmodel = setParam(dmodel,'lb',diet.lipids.rxns,-diet.lipids.ub(:,j));
        dmodel = setParam(dmodel,'lb',diet.carbs.rxns,-diet.carbs.ub(:,j));
        
        %dmodel = setParam(dmodel,'lb',diet.carbs.rxns,-diet.carbs.ub(:,2));
        %dmodel = setParam(dmodel,'lb',diet.lipids.rxns,-diet.lipids.ub(:,2));

        dmodel.ub(abs(dmodel.ub) > 1000) = 1000;
        dmodel.lb(abs(dmodel.lb) > 1000) = -1000;

        disp('----- constraining done -----');
        disp('-----    running FBA    -----');
        
        fba = solveLP(dmodel,1);
        
        fluxVector = zeros(length(dmodel.rxns),3);
        fluxVector(:,1) = dmodel.lb;
        fluxVector(:,2) = dmodel.ub;
        fluxVector(:,3) = fba.x;
 
        T = table(dmodel.rxns,...
                  dmodel.equations,...
                  string(dmodel.subSystems));
        T2 = array2table(fluxVector);
        T = [T T2];
        colheads = {'ID','EQUATION','SUBSYSTEM','LB','UB','FLUX'};
        T.Properties.VariableNames = colheads;
        writetable(T,'output/FBA/'+string(experiment)+'.xlsx','Sheet',char(string(diets(j))+'_'+models(i).modelName),'WriteVariableNames',true);
    end
end

