function [tableFVA, analysed_table, important_reactions,important_wt_out, important_recomb_out] = CompareFluxes_WT_vs_Rec(model)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    wt_model = model;
    %    help fluxVariability
    [minFlux_wt, maxFlux_wt] = fluxVariability(wt_model);
    %table(CDKL5model_constrained_exp_FluxDiff.rxns,CDKL5model_constrained_exp_FluxDiff.rxnNames ,minFlux_biomass, maxFlux_biomass)
    sol_wt = optimizeCbModel(wt_model, 'max');
 
    recomb_model = changeRxnBounds(wt_model, 'RXNbiomass', sol_wt.f*.734  , 'b');
    recomb_model = changeObjective(recomb_model,'EX_cdkl5[c]');
    sol_recomb = optimizeCbModel(recomb_model,'max');
    
    a = sol_wt.x ~= 0;
    sum(a);
	
    fprintf('\nThe number of reactions with flux (FBA) diffrent from zero in  the wt model is %d  \n\n', sum(a)); 

    
    b = sol_recomb.x ~= 0;
    sum(b);

    fprintf('\nThe number of reactions with flux (FBA) diffrent from zero in  the recomb model is %d  \n\n', sum(b)); 
    
   
    %sol_CDKL5 = optimizeCbModel(CDKL5model_constrained_exp_FluxDiff, 'max')
    [minFlux_recomb, maxFlux_recomb] = fluxVariability(recomb_model);
    tableFVA =table(wt_model.rxns,wt_model.rxnNames ,minFlux_recomb, maxFlux_recomb, minFlux_wt, maxFlux_wt, (abs(minFlux_recomb)-abs(maxFlux_recomb)), (abs(minFlux_wt)-abs(maxFlux_wt)) );
    important_wt = [];
    important_recomb = [];

    for i = 1:1:length(sol_wt.v)
        
    if sol_wt.v(i) >= 0
         
            if ((minFlux_wt(i) >= 0.7*sol_wt.v(i)) && (maxFlux_wt(i) <= 1.3*sol_wt.v(i)))
                important_wt(i) = 1;
            else important_wt(i) = 0;
                
            end
    
            
    elseif sol_wt.v(i) < 0
        
             if ((minFlux_wt(i) >= 1.3 * sol_wt.v(i)) && (maxFlux_wt(i) <= 0.7*sol_wt.v(i)))
                  important_wt(i) =1;
             else important_wt(i) = 0;

             end
    end
    
    end
    
    
    for i = 1:1:length(sol_recomb.v)
        
    if sol_recomb.v(i) >=    0
         
            if ((minFlux_recomb(i) >= 0.7*sol_recomb.v(i)) && (maxFlux_recomb(i) <= 1.3*sol_recomb.v(i)))
               important_recomb(i) = 1;
            else important_recomb(i) = 0;
                
            end
    end
    
            
    if sol_recomb.v(i) < 0
        
             if ((minFlux_recomb(i) >= 1.3 * sol_recomb.v(i)) && (maxFlux_recomb(i) <= 0.7*sol_recomb.v(i)))
                  important_recomb(i) =1;
             else important_recomb(i) = 0;
                 
             end
    end
    
    
    end
    
   important_recomb_rxn = recomb_model.rxns(important_recomb'>0);
   important_wt_rxn = wt_model.rxns(important_wt'>0);
   fprintf('\nThe number of active reactions for the recomb model is %d  \n\n', length(important_recomb_rxn)); 
   fprintf('\nThe number of active reactions for the wt model is %d  \n\n', length(important_wt_rxn)); 

          
   important_reactions = intersect(important_wt_rxn, important_recomb_rxn);
   length(important_reactions);
   fprintf('\nThe number of active reactions shared by recomb and wt models  is %d  \n\n', length(important_reactions)); 
    
   important_wt_out = setdiff(important_wt_rxn,important_recomb_rxn);
   important_recomb_out = setdiff(important_recomb_rxn,important_wt_rxn);
   fprintf('\nThere are %d reactions that are active in the wt model and are not active in the recomb model  \n\n', length(important_wt_out)); 
   fprintf('\nThere are %d reactions that are active in the recomb model and are not active in the wt model  \n\n', length(important_recomb_out)); 

   analysed_table = table(tableFVA.Var1(tableFVA.Var7==0 & tableFVA.Var8==0), tableFVA.Var4(tableFVA.Var7==0 & tableFVA.Var8==0), tableFVA.Var5(tableFVA.Var7==0 & tableFVA.Var8==0));


end

