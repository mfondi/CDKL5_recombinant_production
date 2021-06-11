function [tableFVA, analysed_table, important_reactions,important_new_nutrient_out, important_recomb_GG_out,   increasing_flux_reactions , activating_reactions] = CompareFluxes_WT_vs_Rec(model, new_nutrient_rxn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    wt_model = model;
%    new_nutrient_rxn =     'EX_cpd01262_e';

    %    help fluxVariability
    %table(CDKL5model_constrained_exp_FluxDiff.rxns,CDKL5model_constrained_exp_FluxDiff.rxnNames ,minFlux_biomass, maxFlux_biomass)
    sol_new_nutrient = optimizeCbModel(wt_model, 'max');
 
    recomb_model_GG = changeRxnBounds(wt_model, 'RXNbiomass', sol_new_nutrient.f*.74  , 'b');
    recomb_model_GG = changeObjective(recomb_model_GG,'EX_cdkl5[c]');
    sol_recomb_GG = optimizeCbModel(recomb_model_GG,'max');
    sol_recomb_GG.f
    %sol_CDKL5 = optimizeCbModel(CDKL5model_constrained_exp_FluxDiff, 'max')
    [minFlux_recomb, maxFlux_recomb] = fluxVariability(recomb_model_GG)
    
    %% recomb model with new nutrient
    new_nutrient_model = changeRxnBounds(recomb_model_GG, new_nutrient_rxn, -0.5, 'l')
	sol_new_nutrient = optimizeCbModel(new_nutrient_model, 'max');
    sol_new_nutrient.f
    fprintf('\nThe log of ratio is  %d\n\n', log10(sol_new_nutrient.f/sol_recomb_GG.f));

    [minFlux_new_nutrient, maxFlux_wt_new_nutrient] = fluxVariability(new_nutrient_model);

    tableFVA =table(wt_model.rxns,wt_model.rxnNames ,minFlux_recomb, maxFlux_recomb, minFlux_new_nutrient, maxFlux_wt_new_nutrient, (abs(minFlux_recomb)-abs(maxFlux_recomb)), (abs(minFlux_new_nutrient)-abs(maxFlux_wt_new_nutrient)) )
    important_new_nutrient = [];
    important_recomb = [];

    for i = 1:1:length(sol_new_nutrient.v)
        
    if sol_new_nutrient.v(i) >= 0
         
            if ((minFlux_new_nutrient(i) >= 0.7*sol_new_nutrient.v(i)) & (maxFlux_wt_new_nutrient(i) <= 1.3*sol_new_nutrient.v(i)) & (minFlux_new_nutrient(i) ~= 0) & (maxFlux_wt_new_nutrient(i) ~= 0))
               important_new_nutrient(i) = 1;
            else important_new_nutrient(i) = 0;
            end
            
    elseif sol_new_nutrient.v(i) < 0
        
             if ((minFlux_new_nutrient(i) >= 1.3 * sol_new_nutrient.v(i)) & (maxFlux_wt_new_nutrient(i) <= 0.7*sol_new_nutrient.v(i)) & (minFlux_new_nutrient(i) ~= 0) & (maxFlux_wt_new_nutrient(i) ~= 0))
                  important_new_nutrient(i) =1;
             else important_new_nutrient(i) = 0;
             end
    end
    
    end
    
    
    for i = 1:1:length(sol_recomb_GG.v)
        
    if sol_recomb_GG.v(i) >= 0
         
            if ((minFlux_recomb(i) >= 0.7*sol_recomb_GG.v(i)) & (maxFlux_recomb(i) <= 1.3*sol_recomb_GG.v(i)) & (minFlux_recomb(i) ~= 0) & (maxFlux_recomb(i) ~= 0) )
               important_recomb(i) = 1;
            else important_recomb(i) = 0;
            end
            
    elseif sol_recomb_GG.v(i) < 0
        
             if ((minFlux_recomb(i) >= 1.3 * sol_recomb_GG.v(i)) & (maxFlux_recomb(i) <= 0.7*sol_recomb_GG.v(i)) & (minFlux_recomb(i) ~= 0) & (maxFlux_recomb(i) ~= 0))
                  important_recomb(i) =1;
             else important_recomb(i) = 0;
             end
    end
    
    
    end
    
    
   important_recomb_GG_rxn = recomb_model_GG.rxns(important_recomb'>0)
   important_new_nutrient_rxn = wt_model.rxns(important_new_nutrient'>0)
   length(important_recomb_GG_rxn)
   length(important_new_nutrient_rxn)

          
   important_reactions = intersect(important_new_nutrient_rxn, important_recomb_GG_rxn)
   length(important_reactions)
    
   important_new_nutrient_out = setdiff(important_new_nutrient_rxn,important_recomb_GG_rxn)
   important_recomb_GG_out = setdiff(important_recomb_GG_rxn,important_new_nutrient_rxn)

   analysed_table = table(tableFVA.Var1(tableFVA.Var7==0 & tableFVA.Var8==0), tableFVA.Var4(tableFVA.Var7==0 & tableFVA.Var8==0), tableFVA.Var5(tableFVA.Var7==0 & tableFVA.Var8==0))

   idx = ismember(tableFVA.Var1, important_reactions);
   tableFVA_important_reactions =  tableFVA(idx,:);
   tableFVA_important_reactions_out = table(tableFVA_important_reactions.Var1(abs(tableFVA_important_reactions.Var6) > abs(tableFVA_important_reactions.Var4)), tableFVA_important_reactions.Var2(abs(tableFVA_important_reactions.Var6) > abs(tableFVA_important_reactions.Var4)), tableFVA_important_reactions(abs(tableFVA_important_reactions.Var6) > abs(tableFVA_important_reactions.Var4),4), tableFVA_important_reactions(abs(tableFVA_important_reactions.Var6) > abs(tableFVA_important_reactions.Var4),6))
   
   increasing_flux_reactions = tableFVA_important_reactions_out.Var1;
   
   idx_2 = ismember(tableFVA.Var1, important_new_nutrient_out);
   activating_reactions =  tableFVA(idx_2,1);
   

end

