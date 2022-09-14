function v = blockDietRxns(model,rxns,queryRxn)

v = zeros(length(rxns),1);

for k=1:length(rxns)
    kModel = setParam(model,'ub',rxns{k},0);
    fba = solveLP(kModel,1);
    flux = fba.x(getIndexes(kModel,queryRxn,'rxns'));
    v(k,1) = flux;
end