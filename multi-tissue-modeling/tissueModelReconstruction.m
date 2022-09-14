%% IMPORT PATHS
addpath(genpath('/Users/clasenf/OneDrive - The Francis Crick Institute/Integration'));      % path to directory
addpath(genpath('/Users/clasenf/GEM/Raven/'));                                              % RAVEN toolbox
addpath(genpath('/Users/clasenf/GEM/COBRA/'));                                              % COBRA toolbox
changeCobraSolver('gurobi','all');
setRavenSolver('gurobi');

%% IMPORT MMRN

MMRN = importExcelModel('GEM/outModel.xlsx',false);

%% IMPORT TPM EXPRESSION

tpm = readtable('data/transcriptomics/max_tpm_all_tissues.csv');
expression_struct.tissues = tpm.Properties.VariableNames(2:end)';  
expression_struct.genes = tpm.Gene;  
expression_struct.levels = table2array(tpm(:, 2:end));  
expression_struct.threshold = 1;

%% tINIT

essentialTasks = parseTaskList('data/growthTasks.xlsx');
refModel = MMRN;  % the reference model from which the GEM will be extracted
tissue = 'Liver';  % must match the tissue name in data_struct.tissues
celltype = []; % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = expression_struct;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
%taskFile = 'data/growthTasks.xlsx';  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
paramsFT.TimeLimit = 99999999999999;
params.TimeLimit = 99999999999999;
%params = [];  % additional optimization parameters for the INIT algorithm
%paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm

liverGEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, ...
                         removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);


%%

clear
load('GEM/WAT/wat - generic.mat','outModel');
genericWAT = outModel;
save('GEM/WAT/genericWAT.mat','genericWAT');

load('GEM/Kidney/kidney - generic.mat','outModel');
genericKidney = outModel;
save('GEM/Kidney/genericKidney.mat','genericKidney');









