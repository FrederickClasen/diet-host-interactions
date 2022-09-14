function newModel=dietConstraint(model,diet)

tModel = model; % mModel is the merged model 
dietfile = 'data/models/'+diet+'.csv';
% block uptake of all metabolites and add the diet constraint
exchRxns = cellstr(tModel.rxns(contains(tModel.rxns,'EX_')));
tModel = setParam(tModel,'lb',exchRxns,0);
tModel = setParam(tModel,'ub',exchRxns,1000);
diet_constraints = importdata(dietfile);
tModel = addConstraints(tModel,diet_constraints,1);

sinks = {'EX_nmn(e)';'EX_ddca(e)';'EX_pnto_R(e)';'EX_4hbz(e)';'EX_gthox(e)';...
         'EX_12dgr180(e)';'EX_spmd(e)';'EX_sheme(e)';'EX_q8(e)';...
         'EX_pheme(e)';'EX_cgly(e)';'EX_2dmmq8(e)';'EX_26dap_M(e)';...
         'EX_zn2(e)';'EX_so4(e)';'EX_mg2(e)';'EX_mn2(e)';'EX_k(e)';'EX_cu2(e)';...
         'EX_cobalt2(e)';'EX_ca2(e)';'EX_cl(e)';'EX_thm(e)';'EX_fe2(e)';'EX_fe3(e)';...
         'EX_etoh(e)';'EX_fol(e)';'EX_pi(e)';'EX_ribflv(e)';...
         'EX_na1(e)';'EX_thm(e)';'EX_adocbl(e)';'EX_h2o(e)';'EX_zn2(e)';'EX_pydx(e)';...
         'EX_chol(e)';'EX_mqn7(e)';'EX_mqn8(e)';'EX_xan(e)';'EX_Lcystin(e)';'EX_4hpro(e)';...
         'EX_cys_L(e)';'EX_h(e)';'EX_cl(e)';'EX_adn(e)';'EX_amp(e)';'EX_dad_2(e)';...
         'EX_dcyt(e)';'EX_h(e)';'EX_hxan(e)';'EX_ins(e)';'EX_thymd(e)';'EX_ura(e)';...
         'EX_uri(e)';'EX_lac_L(e)';'EX_lac_D(e)';'EX_rib_D(e)';'EX_arab_L(e)';'EX_no3(e)';'EX_so4(e)'};
 
% 'EX_ocdcea(e)' - oleic acid
% 'EX_ttdca(e)' - myristic acid
%'EX_glyc3p(e)';'EX_but(e)';

if ismember('EX_arab_D(e)',tModel.rxns)
    tModel = setParam(tModel,'lb','EX_arab_D(e)',-0.5);
    tModel = setParam(tModel,'ub','EX_arab_D(e)',0);
end
  
tModel = setParam(tModel,'lb',sinks,-0.5);
tModel = setParam(tModel,'ub',sinks,0);

block = {'EX_co2(e)','EX_hdca(e)','EX_ocdca(e)','EX_h2(e)'};
tModel = setParam(tModel,'ub',block,0);

newModel = tModel;

end