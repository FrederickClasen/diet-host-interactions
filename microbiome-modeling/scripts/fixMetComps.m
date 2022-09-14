function newModel=fixMetComps(model)

for i=1:length(model.mets)
    met = model.mets(i);
    for j=1:length(model.comps)
        comp = char('['+string(model.comps(j))+']');
        if contains(met,comp)
            model.metComps(i) = j;
        end
    end
end

newModel = model;