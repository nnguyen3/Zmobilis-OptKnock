%% OptKnock isobutanol workflow May 16 - cleaned full version
% Goal:
% 1. Run OptKnock using EX_ibtol_e as the isobutanol target.
% 2. Remove artifact knockout candidates.
% 3. Manually validate OptKnock candidates.
% 4. Check whether isobutanol production is growth-coupled.
% 5. Force minimum isobutanol production and rerun OptKnock.
% 6. Identify the strongest knockout candidate for growth-associated isobutanol production.

clear
clc

initCobraToolbox(false)
changeCobraSolver('gurobi','LP')
% 1. Load model

N = load('/Users/nhinguyen/Desktop/Z.mobilis/Models/Zm_model_may_06_anaerobic_GF_2026_isobutanol.mat');
model = N.model_AN;

% 2. Confirm isobutanol pathway reactions exist

findRxnIDs(model,'KIVD')
findRxnIDs(model,'IBADH')
findRxnIDs(model,'IBTOLtcp')
findRxnIDs(model,'IBTOLtpe')
findRxnIDs(model,'EX_ibtol_e')

% 3. Set anaerobic medium

model = changeRxnBounds(model, 'EX_glc__D_e', -10, 'b');
model = changeRxnBounds(model, 'EX_nh4_e', -1000, 'l');
model = changeRxnBounds(model, 'EX_nh4_e', 1000, 'u');
model = changeRxnBounds(model, 'EX_o2_e', 0, 'b');

biomassRxn = 'BIOMASS_core';
targetRxn  = 'EX_ibtol_e';

% 4. Check WT growth and max isobutanol
% Purpose:
% Check if the WT model can grow and produce isobutanol.

model = changeObjective(model, biomassRxn);
fbaWT = optimizeCbModel(model);
fprintf('WT growth = %.6f\n', fbaWT.f);

model = changeObjective(model, targetRxn);
fbaIso = optimizeCbModel(model);
fprintf('WT max isobutanol flux = %.6f\n', fbaIso.f);

% Result:
% WT growth = 0.059677
% WT max isobutanol flux = 10.030751


% Interpretation:
% The model can grow anaerobically
% and can theoretically produce isobutanol.


% 5. Broad candidate list
% Purpose:
% Create an initial list of knockout candidates.

isExchangeLike = startsWith(model.rxns,'EX_') | ...
                 startsWith(model.rxns,'DM_') | ...
                 startsWith(model.rxns,'SK_');

selectedRxnList_broad = model.rxns(~isExchangeLike);
selectedRxnList_broad = setdiff(selectedRxnList_broad, {biomassRxn, targetRxn});

% Interpretation:
% This list still contains transport and artifact reactions.

% 6. Growth constraint for initial OptKnock
% Purpose:
% Keep at least 50% WT growth during OptKnock.



constrOpt = struct();
constrOpt.rxnList = {biomassRxn};
constrOpt.values = 0.5 * fbaWT.f;
constrOpt.sense = 'G';

% Interpretation:
% This avoids knockout strategies that stop growth completely.

% 7. Broad 2-KO OptKnock
% Purpose:
% Search for possible 2-reaction knockout strategies.


options = struct();
options.targetRxn = targetRxn;
options.numDel = 2;
options.numDelSense = 'E';

[optKnockSol_2KO_broad, bilevelMILPproblem_2KO_broad] = OptKnock( ...
    model, selectedRxnList_broad, options, constrOpt);

disp('Broad candidate list - recommended 2 knockouts:');
disp(optKnockSol_2KO_broad.rxnList);

growthFlux = optKnockSol_2KO_broad.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = optKnockSol_2KO_broad.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Broad 2-KO predicted growth = %.6f\n', growthFlux);
fprintf('Broad 2-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% CYSItpp + PSURItpp
% growth = 0.069544
% isobutanol flux = 9.450471
% These are transport/periplasm

% Interpretation:
% The results included transport/periplasm reactions,
% so more filtering was needed.

% 8. Broad 1-KO OptKnock
% Purpose:
% Test whether one knockout can improve isobutanol production.

options.numDel = 1;
options.numDelSense = 'E';

[optKnockSol_1KO_broad, bilevelMILPproblem_1KO_broad] = OptKnock( ...
    model, selectedRxnList_broad, options, constrOpt);

disp('Broad candidate list - recommended 1 knockout:');
disp(optKnockSol_1KO_broad.rxnList);

growthFlux = optKnockSol_1KO_broad.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = optKnockSol_1KO_broad.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Broad 1-KO predicted growth = %.6f\n', growthFlux);
fprintf('Broad 1-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% sink_2hxmp_c
% growth = 0.069544
% isobutanol flux = 9.450471
% This is a sink reaction 

% Interpretation:
% The result suggested a sink reaction,
% which is not a good biological target.

% 9. Filter candidate list
% Purpose:
% Remove transport, sink, exchange,
% and other non-ideal reactions.


% =>This creates a cleaner candidate list.

isBadCandidate = startsWith(model.rxns,'EX_') | ...
                 startsWith(model.rxns,'DM_') | ...
                 startsWith(model.rxns,'SK_') | ...
                 startsWith(lower(model.rxns),'sink_') | ...
                 contains(lower(model.rxns),'tex') | ...
                 contains(lower(model.rxns),'tpp') | ...
                 contains(lower(model.rxns),'t2pp') | ...
                 contains(lower(model.rxns),'transport') | ...
                 contains(lower(model.rxns),'biomass') | ...
                 contains(lower(model.rxns),'atps') | ...
                 contains(lower(model.rxns),'abc') | ...
                 contains(lower(model.rxns),'exchange');

selectedRxnList_filtered = model.rxns(~isBadCandidate);
selectedRxnList_filtered = setdiff(selectedRxnList_filtered, {biomassRxn, targetRxn});

% 10. Filtered 1-KO OptKnock
% Purpose:
% Rerun OptKnock using the filtered list.

options.numDel = 1;
options.numDelSense = 'E';

[optKnockSol_1KO_filtered, bilevelMILPproblem_1KO_filtered] = OptKnock( ...
    model, selectedRxnList_filtered, options, constrOpt);

disp('Filtered candidate list - recommended 1 knockout:');
disp(optKnockSol_1KO_filtered.rxnList);

growthFlux = optKnockSol_1KO_filtered.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = optKnockSol_1KO_filtered.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Filtered 1-KO predicted growth = %.6f\n', growthFlux);
fprintf('Filtered 1-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% HPN6
% growth = 0.069544
% isobutanol flux = 9.450471

% Interpretation:
% HPN6 was identified as a possible candidate.


% 11. Filtered 2-KO OptKnock
% Purpose:
% Test whether two filtered knockouts improve production.


options.numDel = 2;
options.numDelSense = 'E';

[optKnockSol_2KO_filtered, bilevelMILPproblem_2KO_filtered] = OptKnock( ...
    model, selectedRxnList_filtered, options, constrOpt);

disp('Filtered candidate list - recommended 2 knockouts:');
disp(optKnockSol_2KO_filtered.rxnList);

growthFlux = optKnockSol_2KO_filtered.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = optKnockSol_2KO_filtered.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Filtered 2-KO predicted growth = %.6f\n', growthFlux);
fprintf('Filtered 2-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% NTP3pp + ECAP3pp
% growth = 0.069544
% isobutanol flux = 9.450471
% These are still periplasm-associated.

% Interpretation:
% Periplasm-related reactions still appeared,
% so more filtering was needed.

% 12. Add pp filter and rerun filtered 2-KO
% Purpose:
% Remove more periplasm-associated reactions.


% => This gave cleaner metabolic candidates.

isBadCandidate = isBadCandidate | contains(lower(model.rxns),'pp');

selectedRxnList_filtered_pp = model.rxns(~isBadCandidate);
selectedRxnList_filtered_pp = setdiff(selectedRxnList_filtered_pp, {biomassRxn, targetRxn});

options.numDel = 2;
options.numDelSense = 'E';

[optKnockSol_2KO_filtered_pp, bilevelMILPproblem_2KO_filtered_pp] = OptKnock( ...
    model, selectedRxnList_filtered_pp, options, constrOpt);

disp('Filtered candidate list with pp removed - recommended 2 knockouts:');
disp(optKnockSol_2KO_filtered_pp.rxnList);

growthFlux = optKnockSol_2KO_filtered_pp.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = optKnockSol_2KO_filtered_pp.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Filtered pp 2-KO predicted growth = %.6f\n', growthFlux);
fprintf('Filtered pp 2-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% G6PDA + MDH3
% growth = 0.069544
% isobutanol flux = 9.450471

% 13. Manual validation of initial candidates
% Purpose:
% Check whether HPN6 and G6PDA + MDH3
% truly improve isobutanol production.

% HPN6
model_HPN6 = changeRxnBounds(model, 'HPN6', 0, 'b');
model_HPN6 = changeObjective(model_HPN6, biomassRxn);
sol_HPN6_growth = optimizeCbModel(model_HPN6);

model_HPN6 = changeObjective(model_HPN6, targetRxn);
sol_HPN6_iso = optimizeCbModel(model_HPN6);

fprintf('Manual HPN6 KO growth = %.6f\n', sol_HPN6_growth.f);
fprintf('Manual HPN6 KO max isobutanol flux = %.6f\n', sol_HPN6_iso.f);

% G6PDA + MDH3
model_G6PDA_MDH3 = model;
model_G6PDA_MDH3 = changeRxnBounds(model_G6PDA_MDH3, 'G6PDA', 0, 'b');
model_G6PDA_MDH3 = changeRxnBounds(model_G6PDA_MDH3, 'MDH3', 0, 'b');

model_G6PDA_MDH3 = changeObjective(model_G6PDA_MDH3, biomassRxn);
sol_G6PDA_MDH3_growth = optimizeCbModel(model_G6PDA_MDH3);

model_G6PDA_MDH3 = changeObjective(model_G6PDA_MDH3, targetRxn);
sol_G6PDA_MDH3_iso = optimizeCbModel(model_G6PDA_MDH3);

fprintf('Manual G6PDA + MDH3 KO growth = %.6f\n', sol_G6PDA_MDH3_growth.f);
fprintf('Manual G6PDA + MDH3 KO max isobutanol flux = %.6f\n', sol_G6PDA_MDH3_iso.f);

% Result:
% HPN6 KO growth = 0.139088
% HPN6 KO max isobutanol flux = 10.000000

% G6PDA + MDH3 KO growth = 0.139088
% G6PDA + MDH3 KO max isobutanol flux = 10.000000

% These did not improve maximum isobutanol production.
% because the wild-type model already reached the same maximum
% isobutanol flux.

% WT max isobutanol = 10.000000
% HPN6 KO max isobutanol = 10.000000
% G6PDA + MDH3 KO max isobutanol = 10.000000

% Therefore, these knockouts did not increase the maximum production
% beyond the WT capacity.




% 14. Check isobutanol flux during normal growth
% Purpose:
% Check whether WT naturally produces isobutanol during growth.



model_WT = changeObjective(model, biomassRxn);
sol_WT = optimizeCbModel(model_WT);

id_iso = findRxnIDs(model_WT, targetRxn);

fprintf('WT growth = %.6f\n', sol_WT.f);
fprintf('WT isobutanol flux during growth = %.6f\n', sol_WT.v(id_iso));

% Result:
% WT growth = 0.059677
% WT isobutanol flux during growth = 0.000000

% The model grows but does not naturally produce isobutanol during growth.

% Interpretation:
% WT growth produced zero isobutanol flux,
% so production was not naturally growth-coupled.

% 15. Validate HPN6 and G6PDA + MDH3 during growth
% Purpose:
% Check whether these knockouts produce isobutanol during growth.

% HPN6 during growth
model_HPN6_growth = changeRxnBounds(model, 'HPN6', 0, 'b');
model_HPN6_growth = changeObjective(model_HPN6_growth, biomassRxn);
sol_HPN6_growth_only = optimizeCbModel(model_HPN6_growth);

id_iso = findRxnIDs(model_HPN6_growth, targetRxn);

fprintf('HPN6 KO growth = %.6f\n', sol_HPN6_growth_only.f);
fprintf('HPN6 KO isobutanol flux during growth = %.6f\n', sol_HPN6_growth_only.v(id_iso));

% G6PDA + MDH3 during growth
model_G6PDA_MDH3_growth = model;
model_G6PDA_MDH3_growth = changeRxnBounds(model_G6PDA_MDH3_growth, 'G6PDA', 0, 'b');
model_G6PDA_MDH3_growth = changeRxnBounds(model_G6PDA_MDH3_growth, 'MDH3', 0, 'b');

model_G6PDA_MDH3_growth = changeObjective(model_G6PDA_MDH3_growth, biomassRxn);
sol_G6PDA_MDH3_growth_only = optimizeCbModel(model_G6PDA_MDH3_growth);

id_iso = findRxnIDs(model_G6PDA_MDH3_growth, targetRxn);

fprintf('G6PDA + MDH3 KO growth = %.6f\n', sol_G6PDA_MDH3_growth_only.f);
fprintf('G6PDA + MDH3 KO isobutanol flux during growth = %.6f\n', sol_G6PDA_MDH3_growth_only.v(id_iso));

% Result:
% HPN6 KO isobutanol flux during growth = 0.000000
% G6PDA + MDH3 KO isobutanol flux during growth = 0.000000

% These knockouts do not couple isobutanol production to growth.
% G6PDA + MDH3, was manually validated under the normal growth condition.
% It did not produce isobutanol during growth, so it was not considered
% a true growth-coupled isobutanol knockout.

% Interpretation:
% => Both knockouts still produced zero isobutanol during growth.

% 16. Test if growth and isobutanol production can occur together
% Purpose:
% Test whether the model can grow
% while producing isobutanol.


for lb = [0.1 0.5 1 2 5 8 9]
    model_test = changeRxnBounds(model, targetRxn, lb, 'l');
    model_test = changeObjective(model_test, biomassRxn);
    sol = optimizeCbModel(model_test);

    id_iso = findRxnIDs(model_test, targetRxn);

    fprintf('EX_ibtol_e LB = %.2f | growth = %.6f | iso flux = %.6f | stat = %d\n', ...
        lb, sol.f, sol.v(id_iso), sol.stat);
end

% Result:
% EX_ibtol_e LB = 0.10 | growth = 0.059439 | iso flux = 0.100000 | stat = 1
% EX_ibtol_e LB = 0.50 | growth = 0.058491 | iso flux = 0.500000 | stat = 1
% EX_ibtol_e LB = 1.00 | growth = 0.057305 | iso flux = 1.000000 | stat = 1
% EX_ibtol_e LB = 2.00 | growth = 0.054933 | iso flux = 2.000000 | stat = 1
% EX_ibtol_e LB = 5.00 | growth = 0.047817 | iso flux = 5.000000 | stat = 1
% EX_ibtol_e LB = 8.00 | growth = 0.040702 | iso flux = 8.000000 | stat = 1
% EX_ibtol_e LB = 9.00 | growth = 0.038330 | iso flux = 9.000000 | stat = 1

% Isobutanol production is feasible during growth,
% but it is not naturally growth-coupled.

% Interpretation:
% Growth was still possible while producing isobutanol.

% 17. Force minimum isobutanol production and rerun OptKnock

% Purpose:
% Force the model to produce isobutanol during growth
% and search for better knockouts.


model_OK = changeRxnBounds(model, targetRxn, 1, 'l');
model_OK = changeObjective(model_OK, biomassRxn);

sol_force = optimizeCbModel(model_OK);
fprintf('Forced iso growth = %.6f\n', sol_force.f);
%Forced iso growth = 0.057305
% Forced WT isobutanol flux = 1.000000

% Interpretation:
% This helps identify growth-associated production strategies.


% Build a cleaner candidate list for the forced isobutanol condition.
% I kept only better knockout candidates and removed reactions that are
% not good biological targets, such as exchange, sink, transport,
% periplasm, biomass, ATP-related, and reactions without gene rules.

% Build clean candidate list with gene rules only
isBadCandidate_force = startsWith(model.rxns,'EX_') | ...
                       startsWith(model.rxns,'DM_') | ...
                       startsWith(model.rxns,'SK_') | ...
                       startsWith(lower(model.rxns),'sink_') | ...
                       contains(lower(model.rxns),'tex') | ...
                       contains(lower(model.rxns),'tpp') | ...
                       contains(lower(model.rxns),'t2pp') | ...
                       contains(lower(model.rxns),'transport') | ...
                       contains(lower(model.rxns),'biomass') | ...
                       contains(lower(model.rxns),'atps') | ...
                       contains(lower(model.rxns),'abc') | ...
                       contains(lower(model.rxns),'exchange') | ...
                       contains(lower(model.rxns),'pp');

selectedRxnList_force = model.rxns(~isBadCandidate_force);
selectedRxnList_force = setdiff(selectedRxnList_force, {biomassRxn, targetRxn});


% Keep only reactions that have gene rules.
% This removes reactions without associated genes,
% because reactions without genes are not useful
% biological knockout targets.
hasGene = ~cellfun(@isempty, model.grRules);
selectedRxnList_force = selectedRxnList_force( ...
    hasGene(ismember(model.rxns, selectedRxnList_force)));

% Rerun OptKnock under the forced isobutanol condition.
% Here, the model is required to produce at least 1 unit of isobutanol,
% and OptKnock searches for one knockout that can improve production
% while still keeping at least 90% of the forced-growth value.
options = struct();
options.targetRxn = targetRxn;
options.numDel = 1;
options.numDelSense = 'E';

constrOpt_force = struct();
constrOpt_force.rxnList = {biomassRxn};
constrOpt_force.values = 0.9 * sol_force.f;
constrOpt_force.sense = 'G';

[result_force, bilevelMILPproblem_force] = OptKnock( ...
    model_OK, selectedRxnList_force, options, constrOpt_force);

disp('Forced iso condition - recommended 1 knockout:');
disp(result_force.rxnList);

% Result:
% Forced iso growth = 0.057305
% Recommended 1-KO = ALCD2x

% Interpretation:
% Under the forced isobutanol condition, OptKnock recommended ALCD2x
% as the best single knockout candidate.
% This means ALCD2x may help redirect flux toward isobutanol production
% while still allowing the model to grow.

% 18. Manual validation of ALCD2x
% Purpose:
% Manually knock out ALCD2x to check if the OptKnock result is real.

model_ALCD2x = model_OK;
model_ALCD2x = changeRxnBounds(model_ALCD2x, 'ALCD2x', 0, 'b');

model_ALCD2x = changeObjective(model_ALCD2x, biomassRxn);
sol_ALCD2x = optimizeCbModel(model_ALCD2x);

id_iso = findRxnIDs(model_ALCD2x, targetRxn);

fprintf('ALCD2x KO growth with forced iso = %.6f\n', sol_ALCD2x.f);
fprintf('ALCD2x KO iso flux during growth = %.6f\n', sol_ALCD2x.v(id_iso));

% Result:
% ALCD2x KO growth = 0.036737
% ALCD2x KO isobutanol flux during growth = 9.671502

% ALCD2x strongly increased isobutanol production during growth.

% Interpretation:
% This step confirms whether ALCD2x actually increases isobutanol
% production during growth, instead of only being an OptKnock prediction.

% 19. Check ALCD2x information

printRxnFormula(model, 'ALCD2x')
model.rxnNames(findRxnIDs(model,'ALCD2x'))
model.grRules(findRxnIDs(model,'ALCD2x'))

% Result:
% ALCD2x KO growth = 0.036737
% ALCD2x KO isobutanol flux during growth = 9.671502

% Interpretation:
% This is a strong result. The model still grows, and isobutanol
% production increases from the forced minimum of 1.000000 to 9.671502.
% Therefore, ALCD2x is a strong knockout candidate for growth-associated
% isobutanol production.

% 20. Try 2-KO under forced isobutanol condition
% Purpose:
% Test whether two knockouts can improve isobutanol production
% more than the ALCD2x single knockout.

options.numDel = 2;
options.numDelSense = 'E';

[result_force_2KO, bilevelMILPproblem_force_2KO] = OptKnock( ...
    model_OK, selectedRxnList_force, options, constrOpt_force);

disp('Forced iso condition - recommended 2 knockouts:');
disp(result_force_2KO.rxnList);

growthFlux = result_force_2KO.fluxes(strcmp(model.rxns, biomassRxn));
targetFlux = result_force_2KO.fluxes(strcmp(model.rxns, targetRxn));

fprintf('Forced iso 2-KO predicted growth = %.6f\n', growthFlux);
fprintf('Forced iso 2-KO predicted isobutanol flux = %.6f\n', targetFlux);

% Result:
% Even with exact 2-KO requested, OptKnock only returned ALCD2x.
% This suggests ALCD2x is the main beneficial knockout.

% 21. Remove ALCD2x and test backup candidate
% Purpose:
% Remove ALCD2x from the candidate list to see if OptKnock
% can find another possible knockout.

% Interpretation:
% PGCD was found as a backup candidate.
% However, manual validation showed that PGCD only kept isobutanol
% at the forced minimum level of 1.000000.
% Therefore, PGCD is feasible but not a strong improvement.
selectedRxnList_noALCD2x = setdiff(selectedRxnList_force, {'ALCD2x'});

[result_force_noALCD2x, bilevelMILPproblem_force_noALCD2x] = OptKnock( ...
    model_OK, selectedRxnList_noALCD2x, options, constrOpt_force);

disp('Forced iso condition without ALCD2x:');
disp(result_force_noALCD2x.rxnList);

% Result:
% PGCD

model_PGCD = model_OK;
model_PGCD = changeRxnBounds(model_PGCD, 'PGCD', 0, 'b');

model_PGCD = changeObjective(model_PGCD, biomassRxn);
sol_PGCD = optimizeCbModel(model_PGCD);

id_iso = findRxnIDs(model_PGCD, targetRxn);

fprintf('PGCD KO growth with forced iso = %.6f\n', sol_PGCD.f);
fprintf('PGCD KO iso flux during growth = %.6f\n', sol_PGCD.v(id_iso));

% Result:
% PGCD KO growth = 0.055732
% PGCD KO isobutanol flux during growth = 1.000000

% PGCD is feasible, but it does not increase isobutanol above the forced minimum.

% Interpretation:
% PGCD was found as a backup candidate.
% However, manual validation showed that PGCD only kept isobutanol
% at the forced minimum level of 1.000000.
% Therefore, PGCD is feasible but not a strong improvement.

% 22. Final summary
% Final interpretation:
% The WT model can grow anaerobically and can produce isobutanol
% when isobutanol is directly optimized.

% However, during normal growth, WT produces 0 isobutanol.
% This means isobutanol production is possible, but not naturally
% growth-coupled.

% Initial OptKnock candidates, including HPN6 and G6PDA + MDH3,
% did not produce isobutanol during normal growth after manual validation.

% After forcing EX_ibtol_e lower bound to 1, OptKnock suggested ALCD2x.
% Manual validation showed that ALCD2x increased isobutanol flux during
% growth from 1.000000 to 9.671502.

% Therefore, ALCD2x is the strongest knockout candidate found for
% growth-associated isobutanol production.