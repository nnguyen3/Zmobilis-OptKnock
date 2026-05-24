%% FluxRETAP
%the core idea is to force increasing levels of product flux, 
% gradually increase isobutanol production from 0 to the maximum level
% run Flux Variability Analysis (FVA) at each production level
% compare reaction flux ranges between low and high-isobutanol conditions
% identify reactions whose flux ranges change the most
% FluxRETAP ranks reactions based on how different the flux distributions 
% are between the two conditions
% Small overlap between flux ranges suggests the reaction is more important
% for isobutanol production


%% FluxRETAP-inspired analysis for isobutanol in MATLAB
% Goal:
% Force isobutanol production from low to high levels,
% run FVA at each level, and rank reactions by flux change.

clear; clc;

% 1. Load model
N = load('/Users/nhinguyen/Desktop/Z.mobilis/Models/Zm_model_may_06_anaerobic_GF_2026_isobutanol.mat');
model = N.model_AN;

% If your model variable is not named "model", rename it here
% Example:
% model = Zm_model_may_16_anaerobic_GF_2026_isobutanol;

biomassRxn = 'BIOMASS_core';
targetRxn  = 'EX_ibtol_e';

% 2. Set solver
changeCobraSolver('gurobi','LP');

% 3. Set anaerobic medium
model = changeRxnBounds(model,'EX_glc__D_e',-10,'l');
model = changeRxnBounds(model,'EX_o2_e',0,'b');
model = changeRxnBounds(model,'EX_nh4_e',-1000,'l');
model = changeRxnBounds(model,'EX_pi_e',-1000,'l');
model = changeRxnBounds(model,'EX_so4_e',-1000,'l');

% 4. Check normal growth
model_growth = changeObjective(model, biomassRxn);
sol_growth = optimizeCbModel(model_growth,'max');

fprintf('\nWT growth = %.6f\n', sol_growth.f);
% WT growth = 0.059677

% 5. Find maximum isobutanol production
model_iso = changeObjective(model, targetRxn);
sol_iso = optimizeCbModel(model_iso,'max');

maxIso = sol_iso.f;
fprintf('Max isobutanol = %.6f\n', maxIso);
%Max isobutanol = 10.030751

%% 6. Run FVA at increasing isobutanol levels
fractions = [0 0.1 0.25 0.5 0.75 0.9 1.0];

numRxns = length(model.rxns);
numLevels = length(fractions);

minMat = zeros(numRxns,numLevels);
maxMat = zeros(numRxns,numLevels);

for i = 1:numLevels

    model_temp = model;

    forcedIso = fractions(i) * maxIso;

    % force minimum isobutanol production
    model_temp = changeRxnBounds(model_temp,targetRxn,forcedIso,'l');

    fprintf('\nRunning FVA at %.0f%% max isobutanol, forced iso = %.6f\n', ...
        fractions(i)*100, forcedIso);

    [minFlux,maxFlux] = fluxVariability(model_temp,100,'max',model.rxns);

    minMat(:,i) = minFlux;
    maxMat(:,i) = maxFlux;

    fprintf('Finished %.0f%% max isobutanol\n', fractions(i)*100);
end
% Running FVA at 0% max isobutanol, forced iso = 0.000000
% Finished 0% max isobutanol
% 
% Running FVA at 10% max isobutanol, forced iso = 1.003075
% Finished 10% max isobutanol
% 
% Running FVA at 25% max isobutanol, forced iso = 2.507688
% Finished 25% max isobutanol
% 
% Running FVA at 50% max isobutanol, forced iso = 5.015375
% Finished 50% max isobutanol
% 
% Running FVA at 75% max isobutanol, forced iso = 7.523063
% Finished 75% max isobutanol
% 
% Running FVA at 90% max isobutanol, forced iso = 9.027675
% Finished 90% max isobutanol
% 
% Running FVA at 100% max isobutanol, forced iso = 10.030751
% Finished 100% max isobutanol


%% 7. Compare low vs high isobutanol states
lowIdx = [1 2];                      % 0% and 10%
highIdx = [numLevels-1 numLevels];   % 90% and 100%

lowMean = mean((minMat(:,lowIdx) + maxMat(:,lowIdx))/2,2);
highMean = mean((minMat(:,highIdx) + maxMat(:,highIdx))/2,2);

lowRange = mean(maxMat(:,lowIdx) - minMat(:,lowIdx),2);
highRange = mean(maxMat(:,highIdx) - minMat(:,highIdx),2);

deltaFlux = highMean - lowMean;
absDeltaFlux = abs(deltaFlux);
rangeChange = highRange - lowRange;

% 8. Create result table
FluxRETAP_table = table( ...
    model.rxns, ...
    model.rxnNames, ...
    model.grRules, ...
    lowMean, ...
    highMean, ...
    deltaFlux, ...
    absDeltaFlux, ...
    lowRange, ...
    highRange, ...
    rangeChange, ...
    'VariableNames', {'rxn','rxnName','GPR','lowMean','highMean','deltaFlux','absDeltaFlux','lowRange','highRange','rangeChange'});

FluxRETAP_table = sortrows(FluxRETAP_table,'absDeltaFlux','descend');
disp(FluxRETAP_table(1:30,:));
% ALCD2x
% lowMean  = -17.951
% highMean = -0.63016
% delta    = +17.321

%% 9. Save output
writetable(FluxRETAP_table,'FluxRETAP_inspired_isobutanol_results.csv');
fprintf('\nSaved result to FluxRETAP_inspired_isobutanol_results.csv\n');

% 10. Display top 30 reactions
disp(FluxRETAP_table(1:30,:));

%% 11. Check important OptKnock candidates
candidateRxns = {'ALCD2x','HPN6','G6PDA','PGCD','MDH3','EX_ibtol_e'};

fprintf('\nImportant candidate reactions:\n');

for j = 1:length(candidateRxns)
    idx = strcmp(FluxRETAP_table.rxn,candidateRxns{j});

    if any(idx)
        disp(FluxRETAP_table(idx,:));
    else
        fprintf('%s not found in model/table\n', candidateRxns{j});
    end
end

%% Important candidate reaction interpretation:
%
% ALCD2x showed one of the largest flux changes between low- and
% high-isobutanol conditions. The reaction flux decreased strongly
% as isobutanol production increased, suggesting that ethanol
% synthesis competes with isobutanol production.
%
% This result is consistent with the previous OptKnock and
% manual knockout validation results, supporting ALCD2x as
% a promising engineering target.
%
% HPN6 showed zero flux under both low- and high-isobutanol
% conditions, indicating that the reaction was inactive in
% the model solution. This explains why HPN6 knockout had
% little effect during manual validation.
%
% G6PDA showed only a very small flux change, suggesting
% that this reaction is not strongly associated with
% isobutanol production under the tested condition.
%
% PGCD showed a moderate increase in flux as isobutanol
% production increased. However, the change was much smaller
% than the changes observed for ALCD2x and the core
% isobutanol synthesis pathway reactions.
%
% EX_ibtol_e, KIVD, IBADH (core isobutanol pathway), ACLS,KARA1, and DHAD1 (were already present in the model.)
% increased strongly
% as isobutanol production increased, consistent with their
% roles in the isobutanol synthesis pathway.
% These reactions are upstream reactions in valine/leucine/isoleucine metabolism and provide
% the precursor 3-methyl-2-oxobutanoate (3mob_c) for isobutanol production.
% Therefore, their increase at high isobutanol production is biologically
% reasonable and supports that the model is routing flux through the expected
% precursor pathway.
%
% Overall, the FluxRETAP-inspired analysis, together with OptKnock and manual
% validation results, consistently identified ALCD2x as the strongest
% engineering target for improving isobutanol production in the anaerobic Z. mobilis model.
