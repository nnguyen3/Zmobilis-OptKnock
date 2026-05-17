%% Optknock 05/16

%isobutanol pathway
N=load("/Users/nhinguyen/Desktop/Z.mobilis/Models/Zm_model_april_27_anaerobic_GF_2026.mat")
model=N.model_AN

%rxns = findRxnsFromMets(model,'isobald_c'); % KIVD
% fix metadata isobald_c
% isobald_c already exists in the model but had incorrect metadata,
% so its formula and charge are corrected here before using it in reactions

% add/fix isobutyraldehyde
if findMetIDs(model,'isobald_c') == 0
    model = addMetabolite(model,'isobald_c','Isobutyraldehyde','C4H8O','','','','',0,0);
else
    idx = find(strcmp(model.mets,'isobald_c'));
    model.metNames{idx} = 'Isobutyraldehyde';
    model.metFormulas{idx} = 'C4H8O';
    model.metCharges(idx) = 0;
end
% rxn1 : 2-ketoisovalerate → isobutyraldehyde + CO2
% rxn2 : isobutyraldehyde + NADH + H+ <=> isobutanol + NAD+ 
% rxn2 (note) :(a primary alcohol + NAD+ = an aldehyde + NADH + H+  )

% add rxn1/ Formula : (CH3)2CHC(O)CO2H → (CH3)2CHCHO + CO2 (wiki)
% this is the decarboxylation step before isobutanol formation
% 2-ketoisovalerate = 3-methyl-2-oxobutanoate in bigg = 3mob = C5H7O3
% => 3mob_c → isobutyraldehyde + CO2 
% add metabolie 3mob_c aka 3-methyl-2-oxobutanoate / 2-ketoisovalerate
findMetIDs(model,'3mob_c') % our model has this met = no need to add 
%model = addMetabolite(model,'3mob_c','2-ketoisovalerate','C5H7O3','','','','',-1,0);

%add rxn1 (kdcA / KIVD) , i made up a name kivd because bigg does not have
%this rxn
rxnName = {'KIVD','2-ketoisovalerate decarboxylase'};
metaboliteList = {'3mob_c','isobald_c','co2_c'}; % doi: 10.1186/s13068-025-02687-6
stoichCoeffList = [-1 1 1];

revFlag = false;
lowerBound = 0;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Valine, leucine and isoleucine metabolism';
grRule = '';
geneNameList = '';
systNameList = '';
checkDuplicate = true;
confScores = {1};

model.csense(length(model.mets),1)='E';
model = addReaction2(model,rxnName,metaboliteList,stoichCoeffList, ...
    revFlag,lowerBound,upperBound,objCoeff,subSystem, ...
    grRule,geneNameList,systNameList,checkDuplicate,confScores);


% Formula rxn2 : (CH3)2CHCHO+ NADH + H -> (CH3)2CHCH2OH + NAD (wikipedia)
% KEGG : a primary alcohol + NAD+ = an aldehyde + NADH + H+  
% there for (CH3)2CHCHO = isobutyraldehyde = C4H8O
% and (CH3)2CHCH2OH = Isobutanol = C4H10O

%add metabolite for isobutanol
if findMetIDs(model,'ibtol_c') == 0
    model = addMetabolite(model,'ibtol_c','Isobutanol','C4H10O','','','','',0,0);
else
    idx = find(strcmp(model.mets,'ibtol_c'));
    model.metNames{idx} = 'Isobutanol';
    model.metFormulas{idx} = 'C4H10O';
    model.metCharges(idx) = 0;
end
%add metabolite for isobutyraldehyde
if findMetIDs(model,'isobald_c') ~= 0
    idx = find(strcmp(model.mets,'isobald_c'));
    model.metNames{idx} = 'Isobutyraldehyde';
    model.metFormulas{idx} = 'C4H8O';
    model.metCharges(idx) = 0;
else
    model = addMetabolite(model,'isobald_c','Isobutyraldehyde','C4H8O','','','','',0,0);
end
%add reaction: isobutyraldehyde + NADH + H+ <=> isobutanol + NAD+
% E.C 1.1.1.1 alcohol dehydrogenase 
% BiGG doesn't have this exact reaction because ADH works on many aldehydes,
% so they only show it as a general reaction.
% This one is still real (based on EC 1.1.1.1), so I added it manually
% with a BiGG-like name to make the model produce isobutanol

rxnName = {'IBADH','Isobutyraldehyde reductase'};
metaboliteList = {'isobald_c','nadh_c','h_c','ibtol_c','nad_c'};
stoichCoeffList = [-1 -1 -1 1 1];
revFlag = true;
lowerBound = -1000;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Valine, leucine and isoleucine metabolism';
grRule = '';
geneNameList = '';
systNameList = '';
checkDuplicate = true;
confScores = {1};

model.csense(length(model.mets),1)='E';
model = addReaction2(model,rxnName,metaboliteList,stoichCoeffList, ...
    revFlag,lowerBound,upperBound,objCoeff,subSystem, ...
    grRule,geneNameList,systNameList,checkDuplicate,confScores);
% add rxn tranport c->p
model = addMetabolite(model,'ibtol_p','Isobutanol','C4H10O','','','','',0,0);
rxnName = {'IBTOLtcp','Isobutanol transport cytosol to periplasm'};
metaboliteList = {'ibtol_c','ibtol_p'};
stoichCoeffList = [-1 1];
revFlag = true;
lowerBound = -1000;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Transport';
grRule = '';
geneNameList = '';
systNameList = '';
checkDuplicate = true;
confScores = {1};

model.csense(length(model.mets),1)='E';
model = addReaction2(model,rxnName,metaboliteList,stoichCoeffList, ...
    revFlag,lowerBound,upperBound,objCoeff,subSystem, ...
    grRule,geneNameList,systNameList,checkDuplicate,confScores);
% add rxn tranport p->e
model = addMetabolite(model,'ibtol_e','Isobutanol','C4H10O','','','','',0,0);
rxnName = {'IBTOLtpe','Isobutanol transport'};
metaboliteList = {'ibtol_p','ibtol_e'};
stoichCoeffList = [-1 1];
revFlag = true;          
lowerBound = -1000;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Transport';
grRule = '';
geneNameList = '';
systNameList = '';
checkDuplicate = true;
confScores = {1};

model.csense(length(model.mets),1)='E';
model = addReaction2(model,rxnName,metaboliteList,stoichCoeffList, ...
    revFlag,lowerBound,upperBound,objCoeff,subSystem, ...
    grRule,geneNameList,systNameList,checkDuplicate,confScores);
% add exchange rxn for isobutanol
rxnName = {'EX_ibtol_e','Isobutanol exchange'};
metaboliteList = {'ibtol_e'};
stoichCoeffList = [-1];
revFlag = false;
lowerBound = 0;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Exchange';
grRule = '';
geneNameList = '';
systNameList = '';
checkDuplicate = true;
confScores = {1};

model.csense(length(model.mets),1)='E';
model = addReaction2(model,rxnName,metaboliteList,stoichCoeffList, ...
    revFlag,lowerBound,upperBound,objCoeff,subSystem, ...
    grRule,geneNameList,systNameList,checkDuplicate,confScores);

% check if they are added correctly
printRxnFormula(model,'KIVD') % rxn1 to produce isobutyraldehyde 
printRxnFormula(model,'IBADH') % rxn2 to produce isobutanol
printRxnFormula(model,'IBTOLtpe') % transport isobutanol p-e
printRxnFormula(model,'IBTOLtcp') % transport isobutanol c-p
printRxnFormula(model,'EX_ibtol_e') % exchange so isobutanol go out to dead end

model_AN = model;

save('/Users/nhinguyen/Desktop/Z.mobilis/Models/Zm_model_may_16_anaerobic_GF_2026_isobutanol.mat', ...
     'model_AN');
%% OptKnock ethanol benchmark for April 27 anaerobic model
% Goal:
% 1. Run OptKnock using ethanol exchange (EX_etoh_e) as the target reaction.
% 2. Check whether LDH_D exists and is included in the knockout candidate list.
% 3. Compare broad vs filtered candidate lists.

% 0. Clear workspace and setup COBRA solver
clear
clc

initCobraToolbox(false)
changeCobraSolver('gurobi','LP')

% 1. Load model
% Anaerobic model with isobutanol pathway added back.(using model April 27 with isobutanol added)

N = load('/Users/nhinguyen/Desktop/Z.mobilis/Models/Zm_model_may_16_anaerobic_GF_2026_isobutanol.mat');
model = N.model_AN;

% 2. Set anaerobic medium
% Glucose uptake = -10
% Ammonium uptake is opened
% Oxygen is blocked

model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'b');
model = changeRxnBounds(model, 'EX_nh4_e', -1000, 'l');
model = changeRxnBounds(model, 'EX_nh4_e', 1000, 'u');
model = changeRxnBounds(model, 'EX_o2_e', 0, 'b');


% 3. Define biomass and target reaction
% Biomass is used as the growth objective.
% EX_etoh_e is the ethanol exchange reaction.

biomassRxn = 'BIOMASS_core';
targetRxn  = 'EX_etoh_e';

% 4. Calculate wild-type growth
model = changeObjective(model, biomassRxn);
sol_WT_growth = optimizeCbModel(model);

model = changeObjective(model, targetRxn);
sol_WT_etoh = optimizeCbModel(model,'max');

fprintf('WT growth = %.6f\n', sol_WT_growth.f);
fprintf('WT max ethanol flux = %.6f\n', sol_WT_etoh.f);

% WT growth = 0.057545
% WT max ethanol flux = 20.102276

% 5. Build broad candidate knockout list
% Remove exchange, demand, and sink reactions.
% Also remove biomass and target reactions from knockout candidates.
% (biomass_core, EX_etoh_e) , we dont want to remove these 

isExchangeLike = startsWith(model.rxns, 'EX_') | ...
                 startsWith(model.rxns, 'DM_') | ...
                 startsWith(model.rxns, 'SK_');

selectedRxnList = model.rxns(~isExchangeLike);
selectedRxnList = setdiff(selectedRxnList, {biomassRxn, targetRxn});


% 6. Set OptKnock options for 2 knockouts
options = struct();
options.targetRxn = targetRxn;
% allow up to 2 knockouts
options.numDel = 2;
% less than or equal to 2 deletions
options.numDelSense = 'L';


% 7. Add growth constraint
% Mutant must maintain at least 50% of WT growth.

constrOpt = struct();
constrOpt.rxnList = {biomassRxn};
constrOpt.values = 0.5 * fbaWT.f;
constrOpt.sense = 'G';


% 8. Run OptKnock with broad candidate list
[optKnockSol_2KO_broad, bilevelMILPproblem_2KO_broad] = OptKnock( ...
    model, selectedRxnList, options, constrOpt);


% 9. Display broad 2-KO result
disp('Broad candidate list - recommended 2 knockouts:');
% Broad candidate list - recommended 2 knockouts:{'T2DECAI' }{'G6PDH2xr'}
disp(optKnockSol_2KO_broad.rxnList);

growthFlux = optKnockSol_2KO_broad.fluxes(strcmp(model.rxns, biomassRxn));

targetFlux = optKnockSol_2KO_broad.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Broad 2-KO predicted mutant growth = %.6f\n', growthFlux);

fprintf('Broad 2-KO predicted ethanol flux = %.6f\n', targetFlux);
% Broad 2-KO predicted mutant growth = 0.028772
% Broad 2-KO predicted ethanol flux = 19.660772

% 10. Run single knockout with broad candidate list
% Diagnostic test using only one knockout.

options.numDel = 1;
options.numDelSense = 'L';

[optKnockSol_1KO_broad, bilevelMILPproblem_1KO_broad] = OptKnock( ...
    model, selectedRxnList, options, constrOpt);

disp('Broad candidate list - recommended single knockout:');
disp(optKnockSol_1KO_broad.rxnList);
%Broad candidate list - recommended single knockout:{'Htex'}
growthFlux = optKnockSol_1KO_broad.fluxes(strcmp(model.rxns, biomassRxn));

targetFlux = optKnockSol_1KO_broad.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Broad 1-KO predicted mutant growth = %.6f\n', growthFlux);

fprintf('Broad 1-KO predicted ethanol flux = %.6f\n', targetFlux);
% Broad 1-KO predicted mutant growth = 0.073869
% Broad 1-KO predicted ethanol flux = 18.771976

% 11. Check whether LDH_D exists and is included
% Check:
% 1. Does LDH_D exist in the model?
% 2. What is the reaction formula?
% 3. Is LDH_D included in the candidate knockout list?

disp('Checking LDH_D:');
ldhID = findRxnIDs(model, 'LDH_D');
fprintf('LDH_D reaction ID = %d\n', ldhID);
printRxnFormula(model, 'LDH_D');
isLDHcandidate = ismember('LDH_D', selectedRxnList); %1 , Yes


% 12. Manually test LDH_D knockout
% Test whether LDH_D knockout is feasible and ethanol-related.

model_LDH = changeRxnBounds(model, 'LDH_D', 0, 'b');

model_LDH = changeObjective(model_LDH, biomassRxn);

sol_LDH_growth = optimizeCbModel(model_LDH);

model_LDH = changeObjective(model_LDH, targetRxn);

sol_LDH_etoh = optimizeCbModel(model_LDH);

fprintf('Manual LDH_D KO growth = %.6f\n', sol_LDH_growth.f);

fprintf('Manual LDH_D KO max ethanol flux = %.6f\n', sol_LDH_etoh.f);
% Manual LDH_D KO growth = 0.057545
% Manual LDH_D KO max ethanol flux = 20.102276

% 13. Build filtered candidate list
% Remove transport-like and ATP synthase-like reactions
% to reduce artifact knockout solutions.

isBadCandidate = startsWith(model.rxns,'EX_') | ...   % exchange reactions
                 startsWith(model.rxns,'DM_') | ...   % demand reactions
                 startsWith(model.rxns,'SK_') | ...   % sink reactions
                 contains(lower(model.rxns),'tex') | ... % extracellular transport reactions
                 contains(lower(model.rxns),'tpp') | ... % periplasm transport reactions
                 contains(lower(model.rxns),'t2pp') | ... % proton/periplasm transport reactions
                 contains(lower(model.rxns),'transport') | ... % transport-related reactions
                 contains(lower(model.rxns),'biomass') | ... % biomass reactions
                 contains(lower(model.rxns),'atps'); % ATP synthase / energy-related reactions
selectedRxnList_filtered = model.rxns(~isBadCandidate);

selectedRxnList_filtered = setdiff(selectedRxnList_filtered, {targetRxn});


% 14. Run single knockout with filtered candidate list
options.numDel = 1;
options.numDelSense = 'L';

[optKnockSol_1KO_filtered, bilevelMILPproblem_1KO_filtered] = OptKnock( ...
    model, selectedRxnList_filtered, options, constrOpt);

disp('Filtered candidate list - recommended single knockout:');
% Filtered candidate list - recommended single knockout:{'GTHS'}
disp(optKnockSol_1KO_filtered.rxnList);

growthFlux = optKnockSol_1KO_filtered.fluxes(strcmp(model.rxns, biomassRxn));

targetFlux = optKnockSol_1KO_filtered.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Filtered 1-KO predicted mutant growth = %.6f\n', growthFlux);

fprintf('Filtered 1-KO predicted ethanol flux = %.6f\n', targetFlux);
% Filtered 1-KO predicted mutant growth = 0.089444
% Filtered 1-KO predicted ethanol flux = 18.602561

% 15. Run 2 knockouts with filtered candidate list
options.numDel = 2;
options.numDelSense = 'L';

[optKnockSol_2KO_filtered, bilevelMILPproblem_2KO_filtered] = OptKnock( ...
    model, selectedRxnList_filtered, options, constrOpt);

disp('Filtered candidate list - recommended 2 knockouts:');
% Filtered candidate list - recommended 2 knockouts:{'GLUDy'}{'GLUSy'}
disp(optKnockSol_2KO_filtered.rxnList);

growthFlux = optKnockSol_2KO_filtered.fluxes(strcmp(model.rxns, biomassRxn));

targetFlux = optKnockSol_2KO_filtered.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Filtered 2-KO predicted mutant growth = %.6f\n', growthFlux);

fprintf('Filtered 2-KO predicted ethanol flux = %.6f\n', targetFlux);
% Filtered 2-KO predicted mutant growth = 0.046769
% Filtered 2-KO predicted ethanol flux = 18.845292

% 16. Summary notes
% The model was feasible under anaerobic conditions.
% Wild-type growth = 0.057545.

% OptKnock successfully ran using EX_etoh_e
% as the ethanol target reaction.

% Broad candidate list results:
% 2-KO result = T2DECAI + G6PDH2xr
% growth = 0.028772
% max ethanol flux = 19.660772

% 1-KO result = Htex
% growth = 0.073869
% max ethanol flux = 18.771976

% Some results were transport-like reactions,
% so a filtered candidate list was also tested.

% Filtered candidate list results:
%1-KO result = GTHS
%growth = 0.089444
%max ethanol flux = 18.602561

%2-KO result = GLUDy + GLUSy
%growth = 0.046769
%max ethanol flux = 18.845292

% LDH_D exists in the model.
% LDH_D is included in the OptKnock candidate list:
% ismember('LDH_D', selectedRxnList) = 1

% Manual LDH_D knockout was feasible:
% growth = 0.057545
% max ethanol flux = 20.102276

% Overall:
% The OptKnock workflow appears to be working correctly.
% LDH_D is present and can be used as a knockout candidate.
% However, under the current model and constraints,
% OptKnock found other knockout solutions instead of LDH_D.


% Why manual LDH_D ethanol > OptKnock ethanol
% Manual LDH_D knockout gave higher ethanol flux
% than the OptKnock-selected solutions.
% This suggests that LDH_D is a valid ethanol-related knockout,
% but OptKnock identified alternative feasible solutions
% under the current search space and constraints.

% Manual broad 2-KO: T2DECAI + G6PDH2xr
model_broad2 = model;
model_broad2 = changeRxnBounds(model_broad2, 'T2DECAI', 0, 'b');
model_broad2 = changeRxnBounds(model_broad2, 'G6PDH2xr', 0, 'b');

model_broad2 = changeObjective(model_broad2, biomassRxn);
sol_broad2_growth = optimizeCbModel(model_broad2);

model_broad2 = changeObjective(model_broad2, targetRxn);
sol_broad2_etoh = optimizeCbModel(model_broad2,'max');

fprintf('Manual T2DECAI + G6PDH2xr KO growth = %.6f\n', sol_broad2_growth.f);
fprintf('Manual T2DECAI + G6PDH2xr KO max ethanol flux = %.6f\n', sol_broad2_etoh.f);
% Manual T2DECAI + G6PDH2xr KO growth = 0.057545
% Manual T2DECAI + G6PDH2xr KO max ethanol flux = 20.102276

% Manual filtered 2-KO: GLUDy + GLUSy
model_filtered2 = model;
model_filtered2 = changeRxnBounds(model_filtered2, 'GLUDy', 0, 'b');
model_filtered2 = changeRxnBounds(model_filtered2, 'GLUSy', 0, 'b');

model_filtered2 = changeObjective(model_filtered2, biomassRxn);
sol_filtered2_growth = optimizeCbModel(model_filtered2);

model_filtered2 = changeObjective(model_filtered2, targetRxn);
sol_filtered2_etoh = optimizeCbModel(model_filtered2,'max');

fprintf('Manual GLUDy + GLUSy KO growth = %.6f\n', sol_filtered2_growth.f);
fprintf('Manual GLUDy + GLUSy KO max ethanol flux = %.6f\n', sol_filtered2_etoh.f);

% Manual GLUDy + GLUSy KO growth = 0.029418
% Manual GLUDy + GLUSy KO max ethanol flux = 20.102276
% 6. Manual broad 1-KO: Htex
model_Htex = changeRxnBounds(model, 'Htex', 0, 'b');

model_Htex = changeObjective(model_Htex, biomassRxn);
sol_Htex_growth = optimizeCbModel(model_Htex);

model_Htex = changeObjective(model_Htex, targetRxn);
sol_Htex_etoh = optimizeCbModel(model_Htex,'max');

fprintf('Manual Htex KO growth = %.6f\n', sol_Htex_growth.f);
fprintf('Manual Htex KO max ethanol flux = %.6f\n', sol_Htex_etoh.f);

% Manual Htex KO growth = 0.046464
% Manual Htex KO max ethanol flux = 20.061667

% 7. Manual filtered 1-KO: GTHS
model_GTHS = changeRxnBounds(model, 'GTHS', 0, 'b');

model_GTHS = changeObjective(model_GTHS, biomassRxn);
sol_GTHS_growth = optimizeCbModel(model_GTHS);

model_GTHS = changeObjective(model_GTHS, targetRxn);
sol_GTHS_etoh = optimizeCbModel(model_GTHS,'max');

fprintf('Manual GTHS KO growth = %.6f\n', sol_GTHS_growth.f);
fprintf('Manual GTHS KO max ethanol flux = %.6f\n', sol_GTHS_etoh.f);

% Manual GTHS KO growth = 0.000000
% Manual GTHS KO max ethanol flux = 20.102276


% WT growth = 0.057545
% WT max ethanol flux ≈ 20.102276

% T2DECAI + G6PDH2xr:
% growth = 0.057545
% ethanol = 20.102276
% => no improvement, but feasible

% GLUDy + GLUSy:
% growth = 0.029418
% ethanol = 20.102276
% => feasible, but growth lower

% Htex:
% growth = 0.046464
% ethanol = 20.061667
% => feasible, ethanol slightly lower

% GTHS:
% growth = 0.000000
% ethanol = 20.102276
% => not useful because no growth

% LDH_D KO growth = 0.057545
% LDH_D KO max ethanol = 20.102276
% LDH_D knockout was feasible and maintained biomass growth.
% The maximum ethanol flux was 20.102276, which was similar to the WT


%% Force ethanol secretion before OptKnock

model_force = model;

% Force small ethanol export
model_force = changeRxnBounds(model_force, 'EX_etoh_e', 1, 'l');

% Check growth with forced ethanol
model_force = changeObjective(model_force, biomassRxn);
sol_force_growth = optimizeCbModel(model_force);

% Check ethanol flux during growth
etohFlux_force = sol_force_growth.x(strcmp(model_force.rxns, targetRxn));

fprintf('Forced ethanol growth = %.6f\n', sol_force_growth.f);
fprintf('Forced ethanol flux during growth = %.6f\n', etohFlux_force);

% max ethanol under forced condition
model_force = changeObjective(model_force, targetRxn);
sol_force_etoh = optimizeCbModel(model_force, 'max');

fprintf('Forced condition max ethanol = %.6f\n', sol_force_etoh.f);

% Forced ethanol growth = 0.057545
% Forced ethanol flux during growth = 1.000000
% Forced condition max ethanol = 20.102276

for etohLB = [5 10 15 18 19 20]
    model_test = model;
    model_test = changeRxnBounds(model_test, 'EX_etoh_e', etohLB, 'l');

    model_test = changeObjective(model_test, biomassRxn);
    sol = optimizeCbModel(model_test);

    fprintf('Forced ethanol LB = %.1f | growth = %.6f\n', etohLB, sol.f);
end

% Forced ethanol LB = 5.0 | growth = 0.057545
% Forced ethanol LB = 10.0 | growth = 0.057545
% Forced ethanol LB = 15.0 | growth = 0.057545
% Forced ethanol LB = 18.0 | growth = 0.057545
% Forced ethanol LB = 19.0 | growth = 0.057545
% Ethanol LB = 19 was selected because it was the highest feasible value
% that maintained normal biomass growth.
% Forced ethanol LB = 20.0 | growth = 0.005866, decline

model_force = model;
model_force = changeRxnBounds(model_force, 'EX_etoh_e', 19, 'l');

model_force = changeObjective(model_force, biomassRxn);
sol_force_growth = optimizeCbModel(model_force);

constrOpt_force = struct();
constrOpt_force.rxnList = {biomassRxn};
constrOpt_force.values = 0.5 * sol_force_growth.f;
constrOpt_force.sense = 'G';

options.numDel = 1;
options.numDelSense = 'L';
options.targetRxn = 'EX_etoh_e';

[optKnockSol_1KO_force19, bilevelMILPproblem_1KO_force19] = OptKnock( ...
    model_force, selectedRxnList_filtered, options, constrOpt_force);

disp(optKnockSol_1KO_force19.rxnList); %{'GLUCYS'}

% Manual validation: GLUCYS KO under forced ethanol condition

model_GLUCYS = model_force;
model_GLUCYS = changeRxnBounds(model_GLUCYS, 'GLUCYS', 0, 'b');

model_GLUCYS = changeObjective(model_GLUCYS, biomassRxn);
sol_GLUCYS_growth = optimizeCbModel(model_GLUCYS);

etohFlux_GLUCYS_growth = sol_GLUCYS_growth.x(strcmp(model_GLUCYS.rxns, targetRxn));

model_GLUCYS = changeObjective(model_GLUCYS, targetRxn);
sol_GLUCYS_etoh = optimizeCbModel(model_GLUCYS, 'max');

fprintf('GLUCYS KO growth = %.6f\n', sol_GLUCYS_growth.f);
fprintf('GLUCYS KO ethanol flux during growth = %.6f\n', etohFlux_GLUCYS_growth);
fprintf('GLUCYS KO max ethanol flux = %.6f\n', sol_GLUCYS_etoh.f);

% GLUCYS KO growth = 0.000000
% GLUCYS KO ethanol flux during growth = 19.000000
% GLUCYS KO max ethanol flux = 20.102276
% GLUCYS was rejected because manual validation showed zero biomass growth,
% even though ethanol export was still feasible.


% Remove GLUCYS from the candidate list and re run 1 OK again
selectedRxnList_force_noGLUCYS = setdiff(selectedRxnList_filtered, {'GLUCYS'});

[optKnockSol_1KO_force19_noGLUCYS, bilevelMILPproblem_1KO_force19_noGLUCYS] = OptKnock( ...
    model_force, selectedRxnList_force_noGLUCYS, options, constrOpt_force);

disp(optKnockSol_1KO_force19_noGLUCYS.rxnList); % {'GTHS'}

% Manual validation: GTHS KO under forced ethanol condition

model_GTHS_force = model_force;
model_GTHS_force = changeRxnBounds(model_GTHS_force, 'GTHS', 0, 'b');

model_GTHS_force = changeObjective(model_GTHS_force, biomassRxn);
sol_GTHS_force_growth = optimizeCbModel(model_GTHS_force);

etohFlux_GTHS_force = sol_GTHS_force_growth.x(strcmp(model_GTHS_force.rxns, targetRxn));

fprintf('GTHS KO forced growth = %.6f\n', sol_GTHS_force_growth.f);
fprintf('GTHS KO forced ethanol flux = %.6f\n', etohFlux_GTHS_force);
% GTHS KO forced growth = 0.000000
% GTHS KO forced ethanol flux = 19.000000
% GTHS was rejected because manual validation showed zero biomass growth
% under the forced ethanol condition.

% Rerun OptKnock after removing  {'GLUCYS','GTHS'});candidates

selectedRxnList_force_noBad = setdiff( ...
    selectedRxnList_filtered, ...
    {'GLUCYS','GTHS'});

options.numDel = 1;
options.numDelSense = 'L';
options.targetRxn = 'EX_etoh_e';

[optKnockSol_1KO_force19_noBad, bilevelMILPproblem_1KO_force19_noBad] = OptKnock( ...
    model_force, selectedRxnList_force_noBad, options, constrOpt_force);

disp('Forced ethanol condition without GLUCYS and GTHS:');
disp(optKnockSol_1KO_force19_noBad.rxnList);

growthFlux = optKnockSol_1KO_force19_noBad.fluxes(strcmp(model_force.rxns, biomassRxn));
etohFlux = optKnockSol_1KO_force19_noBad.fluxes(strcmp(model_force.rxns, targetRxn));

fprintf('Predicted 1-KO growth = %.6f\n', growthFlux);
fprintf('Predicted 1-KO ethanol flux = %.6f\n', etohFlux);

% Forced ethanol condition without GLUCYS and GTHS:{'GLUDy'}
% Predicted 1-KO growth = 0.085104
% Predicted 1-KO ethanol flux = 18.591738

% Manual validation: GLUDy KO under forced ethanol condition

model_GLUDy_force = model_force;
model_GLUDy_force = changeRxnBounds(model_GLUDy_force, 'GLUDy', 0, 'b');

model_GLUDy_force = changeObjective(model_GLUDy_force, biomassRxn);
sol_GLUDy_force_growth = optimizeCbModel(model_GLUDy_force);

etohFlux_GLUDy_force = sol_GLUDy_force_growth.x(strcmp(model_GLUDy_force.rxns, targetRxn));

model_GLUDy_force = changeObjective(model_GLUDy_force, targetRxn);
sol_GLUDy_force_etoh = optimizeCbModel(model_GLUDy_force, 'max');

fprintf('GLUDy KO forced growth = %.6f\n', sol_GLUDy_force_growth.f);
fprintf('GLUDy KO forced ethanol flux during growth = %.6f\n', etohFlux_GLUDy_force);
fprintf('GLUDy KO max ethanol flux = %.6f\n', sol_GLUDy_force_etoh.f);

% GLUDy was kept because manual validation showed positive biomass growth
% while maintaining the forced ethanol flux of 19.

%% Run 2-KO OptKnock under forced ethanol condition

options.numDel = 2;
options.numDelSense = 'L';
options.targetRxn = 'EX_etoh_e';

[optKnockSol_2KO_force19, bilevelMILPproblem_2KO_force19] = OptKnock( ...
    model_force, selectedRxnList_force_noBad, options, constrOpt_force);

disp('Forced ethanol condition - recommended 2 knockouts:');
disp(optKnockSol_2KO_force19.rxnList);

growthFlux = optKnockSol_2KO_force19.fluxes(strcmp(model_force.rxns, biomassRxn));
etohFlux = optKnockSol_2KO_force19.fluxes(strcmp(model_force.rxns, targetRxn));

fprintf('Predicted 2-KO growth = %.6f\n', growthFlux);
fprintf('Predicted 2-KO ethanol flux = %.6f\n', etohFlux);

% OptKnock returned only one knockout candidate (NDPK4) options.numDelSense = 'L';
% even though up to 2 knockouts were allowed.
% Manual validation: NDPK4 KO under forced ethanol condition

model_NDPK4_force = model_force;
model_NDPK4_force = changeRxnBounds(model_NDPK4_force, 'NDPK4', 0, 'b');

model_NDPK4_force = changeObjective(model_NDPK4_force, biomassRxn);
sol_NDPK4_force_growth = optimizeCbModel(model_NDPK4_force);

etohFlux_NDPK4_force = sol_NDPK4_force_growth.x(strcmp(model_NDPK4_force.rxns, targetRxn));

model_NDPK4_force = changeObjective(model_NDPK4_force, targetRxn);
sol_NDPK4_force_etoh = optimizeCbModel(model_NDPK4_force, 'max');

fprintf('NDPK4 KO forced growth = %.6f\n', sol_NDPK4_force_growth.f);
fprintf('NDPK4 KO forced ethanol flux during growth = %.6f\n', etohFlux_NDPK4_force);
fprintf('NDPK4 KO max ethanol flux = %.6f\n', sol_NDPK4_force_etoh.f);

% NDPK4 KO forced growth = 0.057545
% NDPK4 KO forced ethanol flux during growth = 19.000000
% NDPK4 KO max ethanol flux = 20.102276

% With numDelSense = 'L', OptKnock returned only NDPK4 because it allowed
% less than or equal to 2 knockouts.
% NDPK4 was feasible, but it did not improve maximum ethanol production.

% With numDelSense = 'E', OptKnock was forced to return exactly 2 knockouts.
% OptKnock options.numDelSense = 'E'; to see what happen
options.numDel = 2;
options.numDelSense = 'E';
options.targetRxn = 'EX_etoh_e';

[optKnockSol_2KO_force19, bilevelMILPproblem_2KO_force19] = OptKnock( ...
    model_force, selectedRxnList_force_noBad, options, constrOpt_force);

disp('Forced ethanol condition - recommended 2 knockouts:');
disp(optKnockSol_2KO_force19.rxnList);

growthFlux = optKnockSol_2KO_force19.fluxes(strcmp(model_force.rxns, biomassRxn));
etohFlux = optKnockSol_2KO_force19.fluxes(strcmp(model_force.rxns, targetRxn));

fprintf('Predicted 2-KO growth = %.6f\n', growthFlux);
fprintf('Predicted 2-KO ethanol flux = %.6f\n', etohFlux);

% % The recommended 2-KO was GLUDy + GLUSy.

% Predicted 2-KO growth = 0.046769
% Predicted 2-KO ethanol flux = 18.843082

% Manual validation: GLUDy + GLUSy KO under forced ethanol condition

model_GLUDy_GLUSy_force = model_force;
model_GLUDy_GLUSy_force = changeRxnBounds(model_GLUDy_GLUSy_force, 'GLUDy', 0, 'b');
model_GLUDy_GLUSy_force = changeRxnBounds(model_GLUDy_GLUSy_force, 'GLUSy', 0, 'b');

model_GLUDy_GLUSy_force = changeObjective(model_GLUDy_GLUSy_force, biomassRxn);
sol_GLUDy_GLUSy_force_growth = optimizeCbModel(model_GLUDy_GLUSy_force);

etohFlux_GLUDy_GLUSy_force = sol_GLUDy_GLUSy_force_growth.x(strcmp(model_GLUDy_GLUSy_force.rxns, targetRxn));

model_GLUDy_GLUSy_force = changeObjective(model_GLUDy_GLUSy_force, targetRxn);
sol_GLUDy_GLUSy_force_etoh = optimizeCbModel(model_GLUDy_GLUSy_force, 'max');

fprintf('GLUDy + GLUSy KO forced growth = %.6f\n', sol_GLUDy_GLUSy_force_growth.f);
fprintf('GLUDy + GLUSy KO forced ethanol flux during growth = %.6f\n', etohFlux_GLUDy_GLUSy_force);
fprintf('GLUDy + GLUSy KO max ethanol flux = %.6f\n', sol_GLUDy_GLUSy_force_etoh.f);
% GLUDy + GLUSy was feasible under forced ethanol production.
% The mutant maintained ethanol flux = 19, but growth decreased to 0.029418.
% Maximum ethanol flux was unchanged compared to WT, so this KO did not improve ethanol production.

% Ethanol production in the WT model was already near the theoretical maximum.
% Several feasible knockout candidates were identified (GLUDy, NDPK4),
% but none increased the maximum ethanol flux beyond the WT value.

% WT:
% growth = 0.057545
% max ethanol flux = 20.102276

% GLUDy KO:
% growth = 0.053530
% max ethanol flux = 20.102276

% NDPK4 KO:
% growth = 0.057545
% max ethanol flux = 20.102276

% Although LDH_D KO was feasible and had similar maximum ethanol flux,
% OptKnock preferred other knockout solutions under the current optimization constraints.

