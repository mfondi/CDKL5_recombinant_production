function [constrained_model, Rect_constrained_model, sol_wt_constrained_mu, sol_rec_constrained_mu, sol_rec_constrained_cdkl5] = Constraint_model_with_experiments(model, RxnExchangeGlcGln, Glutamate_UR,Gluconate_UR, biomass_ratio_wt_cdkl5)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

constrained_model = changeRxnBounds(model, RxnExchangeGlcGln, [-Glutamate_UR,-Gluconate_UR], 'l');
constrained_model = changeObjective(constrained_model,'RXNbiomass');
sol_constrained_model_WT = optimizeCbModel(constrained_model, 'max', 'one')

fprintf('\nThe production of biomass after constraining the model with experimental data (WT model) is %d \n', sol_constrained_model_WT.f);

reduced_biomass_solution = sol_constrained_model_WT.f*biomass_ratio_wt_cdkl5;

% set the biomass constraint to 99% value and predict CDKL5 flux
constrained_model = changeRxnBounds(constrained_model, 'RXNbiomass', reduced_biomass_solution  , 'b' );
constrained_model = changeObjective(constrained_model,'EX_cdkl5[c]');
sol_constrained_model_CDKL5 = optimizeCbModel(constrained_model, 'max');
fprintf('\nThe production of CDKL5 after constraining the model with experimental data and  biomass production flux is %d \n', sol_constrained_model_CDKL5.f);
constrained_model = changeRxnBounds(constrained_model, 'RXNbiomass', 0 , 'l' );
constrained_model = changeRxnBounds(constrained_model, 'RXNbiomass', 1000 , 'u' );
constrained_model = changeObjective(constrained_model,'EX_cdkl5[c]');
sol_constrained_model_CDKL5_no_biomass_constraint = optimizeCbModel(constrained_model, 'max');
fprintf('\nThe production of CDKL5 after constraining the model with experimental data BUT NOT with biomass production flux (i.e. CDKL theoretycal yield) is %d \n', sol_constrained_model_CDKL5_no_biomass_constraint.f);


%experimentally determined CDKL5 production flux
CDKL5_production_flux_exp = sol_constrained_model_CDKL5.f;

%set the CDKL5 production flux in the model to the 
constrained_model_biomass_reduction = changeRxnBounds(constrained_model, 'RXNbiomass', 0 , 'l' );
constrained_model_biomass_reduction = changeRxnBounds(constrained_model_biomass_reduction, 'RXNbiomass', 1000 , 'u' );
constrained_model_biomass_reduction = changeRxnBounds(constrained_model_biomass_reduction, 'EX_cdkl5[c]', CDKL5_production_flux_exp , 'b' );
constrained_model_biomass_reduction = changeObjective(constrained_model_biomass_reduction,'RXNbiomass');
sol_constrained_model_biomass_reduction = optimizeCbModel(constrained_model_biomass_reduction, 'max', 'one');
fprintf('\nThe production of biomass after constraining the model with experimental data of CDKL produtction is %d \n', sol_constrained_model_biomass_reduction.f);

fprintf('\nThe ratio of biomass production between WT and recombinant strain is %d \n', sol_constrained_model_biomass_reduction.f/sol_constrained_model_WT.f);

constrained_model = changeObjective(constrained_model,'RXNbiomass');
Rect_constrained_model = changeObjective(constrained_model,'EX_cdkl5[c]');
Rect_constrained_model = changeRxnBounds(Rect_constrained_model, 'RXNbiomass', reduced_biomass_solution  , 'b' );

sol_wt_constrained_mu= sol_constrained_model_WT.f;
sol_rec_constrained_mu= sol_constrained_model_biomass_reduction.f;
sol_rec_constrained_cdkl5= sol_constrained_model_CDKL5.f;

end

