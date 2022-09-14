function [outModel] = ratioConstrain(model,fc1,fc2)

% first Eflux constraint
constraints1 = importdata(fc1);
constraints1.textdata(1,:) = [];
constraints1.textdata(:,2) = [];
constraints1.textdata(:,2) = [];
constraints1.rxn = constraints1.textdata;
constraints1.lb = constraints1.data(:,[1]);
constraints1.ub = constraints1.data(:,[2]);

%second Eflux constraint
constraints2 = importdata(fc2);
constraints2.textdata(1,:) = [];
constraints2.textdata(:,2) = [];
constraints2.textdata(:,2) = [];
constraints2.rxn = constraints2.textdata;
constraints2.lb = constraints2.data(:,[1]);
constraints2.ub = constraints2.data(:,[2]);

ubRatio = constraints1.ub./constraints2.ub;
lbRatio = constraints1.lb./constraints2.lb;
lbRatio = -lbRatio;
lbRatio(isnan(lbRatio)) = 0;
ubRatio(isnan(ubRatio)) = 0;

outModel = setParam(model,'ub',constraints1.rxn,ubRatio);
outModel = setParam(outModel,'lb',constraints1.rxn,lbRatio);
















