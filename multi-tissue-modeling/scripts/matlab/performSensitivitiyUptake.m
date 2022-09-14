function msens = performSensitivitiyUptake(model,rxn)

steps = 20;
increment = 1;
tModel = setParam(model,'eq',rxn,0);
msens = zeros(length(model.rxns),steps);
constraint = -increment;
for k=1:steps
    constraint = constraint + increment;
    sModel = setParam(tModel,'eq',rxn,constraint);
    fba = solveLP(sModel,1);
    if fba.stat == 1
        msens(:,k) = fba.x;
    else
        msens(:,k) = zeros(length(tModel.rxns,1));
    end
end

    
