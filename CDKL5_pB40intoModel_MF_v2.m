clear all 
close all
initCobraToolbox;
TAC125Model = readCbModel('iMF721_SMBL3FBC2.xml')

% close all ex reactions
EX_reactions= TAC125Model.rxns(~cellfun(@isempty, regexp(TAC125Model.rxns,'^EX_')));
TAC125Model = changeRxnBounds(TAC125Model, EX_reactions, 0, 'l');
FBAsolutionNoEX =optimizeCbModel(TAC125Model,'max' )

% define GG medium
RxnExchangeSchatz = { 'EX_cpd00254_e' , 'EX_cpd00012_e' , 'EX_cpd00971_e' , 'EX_cpd00067_e' , 'EX_cpd00063_e' , 'EX_cpd00048_e' , 'EX_cpd00205_e' , 'EX_cpd00099_e' , 'EX_cpd00009_e' , 'EX_cpd00013_e' , 'EX_cpd10515_e',  'EX_cpd00209_e' };
TAC125Model = changeRxnBounds(TAC125Model, RxnExchangeSchatz, -1000, 'l');
RxnExchangeGlcGln = {'EX_cpd00222_e' , 'EX_cpd00023_e'};
% definisco i LB dei GG (no sperimentali)
TAC125ModelGG = changeRxnBounds(TAC125Model, RxnExchangeGlcGln, [-10 -10], 'l');
optimizeCbModel(TAC125ModelGG)

%biomassa calcolata con BOF piÃ¹ il plasmide a se 
BOF = addReaction(TAC125ModelGG, 'pB40', '21 cpd00001[c] + 21 cpd00002[c] +  57 cpd00115[c] + 43 cpd00241[c] + 43 cpd00356[c] + 57 cpd00357[c] 	->	plasmid[c] + 21 cpd00008[c] + 21 cpd00009[c] + 21 cpd00067[c]');
BOF = addExchangeRxn(BOF, 'plasmid[c]' , 0,  1000);
BOF = addReaction(BOF, 'RXNbiomass', '0.000341 plasmid[c] + 0.0114 LOS[c] + 40 cpd00001[c] + 40.18733392998440582 cpd00002[c] + 0.181174243687705 cpd00062[c] + 0.13236500641593288 cpd00052[c] + 0.13749453170020254 cpd00038[c] + 0.014829810919920318 cpd00115[c] + 0.010737901824505394 cpd00356[c] + 0.009432412454815721 cpd00241[c] +0.015269223874917404 cpd00357[c] + 0.00215 cpd00003[c] + 5e-06 cpd00004[c] + 0.0004 cpd00005[c] + 0.00013 cpd00006[c] + 5e-06 cpd00010[c] + 5e-06 cpd00015[c] + 0.001 cpd00018[c] + 5e-05 cpd00022[c] + 0.003 cpd00026[c] + 3e-06 cpd00078[c] + 0.035 cpd00118[c] + 0.47 cpd00155[c] + 0.007 cpd00264[c] + 0.05 cpd00345[c] + 0.010152 cpd00908[c] + 6.9e-05 cpd11422[c] + 0.069 cpd11456[c] + 0.000487 cpd11463[c] + 0.0200076 cpd11652[c] + 0.0277812 cpd16661[e]  -> 40 cpd00008[c] + 40 cpd00009[c] + 40 cpd00067[c] + cpd11416[c] +  0.688637060862405 cpd00012[c]');
%BOF = addReaction(BOF, 'CDKL5.c', ' 0.001 plasmid[c] + 59 cpd00035[c] + 6 cpd00084[c]  + 63 cpd00041[c] +  76 cpd00023[c] + 32 cpd00066[c]  + 72 cpd00033[c] + 49 cpd00119[c] + 39 cpd00322[c] + 84 cpd00039[c] + 101 cpd00107[c] + 20 cpd00060[c] + 59 cpd00132[c] + 80 cpd00129[c] + 53 cpd00053[c] + 76 cpd00051[c] + 140 cpd00054[c] + 56 cpd00161[c] + 41 cpd00156[c] + 6 cpd00065[c] + 32 cpd00069[c] + 2288 cpd00002[c] + 2286 cpd00038[c]  ->   0.01 cdkl5[c] + 2288 cpd00018[c] + 2286 cpd00031[c] + 4574 cpd00009[c]');
%new version of CDKL5 
BOF = addReaction(BOF, 'CDKL5.c', ' 0.426 cpd00052[c] + 0.574 cpd00062[c] + 59 cpd00035[c] + 6 cpd00084[c]  + 63 cpd00041[c] +  76 cpd00023[c] + 32 cpd00066[c]  + 72 cpd00033[c] + 49 cpd00119[c] + 39 cpd00322[c] + 84 cpd00039[c] + 101 cpd00107[c] + 20 cpd00060[c] + 59 cpd00132[c] + 80 cpd00129[c] + 53 cpd00053[c] + 76 cpd00051[c] + 140 cpd00054[c] + 56 cpd00161[c] + 41 cpd00156[c] + 6 cpd00065[c] + 32 cpd00069[c] + 2288.574 cpd00002[c] + 2286.426 cpd00038[c]  ->  0.01 cdkl5[c] + 2288 cpd00018[c] + 2286 cpd00031[c] + 4572 cpd00009[c]');
BOF = addExchangeRxn(BOF, 'cdkl5[c]' , 0,  1000);

iMF721_v2= addReaction(TAC125Model, 'RXNbiomass', '0.0114 LOS[c] + 40 cpd00001[c] + 40.18733392998440582 cpd00002[c] + 0.181174243687705 cpd00062[c] + 0.13236500641593288 cpd00052[c] + 0.13749453170020254 cpd00038[c] + 0.014829810919920318 cpd00115[c] + 0.010737901824505394 cpd00356[c] + 0.009432412454815721 cpd00241[c] +0.015269223874917404 cpd00357[c] + 0.00215 cpd00003[c] + 5e-06 cpd00004[c] + 0.0004 cpd00005[c] + 0.00013 cpd00006[c] + 5e-06 cpd00010[c] + 5e-06 cpd00015[c] + 0.001 cpd00018[c] + 5e-05 cpd00022[c] + 0.003 cpd00026[c] + 3e-06 cpd00078[c] + 0.035 cpd00118[c] + 0.47 cpd00155[c] + 0.007 cpd00264[c] + 0.05 cpd00345[c] + 0.010152 cpd00908[c] + 6.9e-05 cpd11422[c] + 0.069 cpd11456[c] + 0.000487 cpd11463[c] + 0.0200076 cpd11652[c] + 0.0277812 cpd16661[e]  -> 40 cpd00008[c] + 40 cpd00009[c] + 40 cpd00067[c] + cpd11416[c] +  0.688637060862405 cpd00012[c]');
iMF721_v2 = changeRxnBounds(iMF721_v2, EX_reactions, -1, 'l');
optimizeCbModel(iMF721_v2)
save('iMF217_v2')


% verify CDKL5 flux 
optCDKL5_TACpB40_CDKL5 = changeObjective(BOF, 'EX_cdkl5[c]');
sol_optCDKL5_TACpB40_CDKL5 = optimizeCbModel(optCDKL5_TACpB40_CDKL5)
% verify Biomass flux
optBIOMASS_TACpB40_CDKL5 = changeObjective(optCDKL5_TACpB40_CDKL5, 'RXNbiomass');
sol_optBIOMASS_TACpB40_CDKL5 = optimizeCbModel(optBIOMASS_TACpB40_CDKL5)


%% costrained Biomass rxn with GG uptake rates from experimental results and 0.99 of biomass production
% load experimental growth data from wt and cdkl4 strain to get biomass
% production ratio in the two strains after 8 hours (4th time point)
Growth_curves = importdata('Growth_curves.txt')

%% Function to compute uptake rates from experimental data
[Glutamate_UR,Gluconate_UR] =  Compute_GG_consumption()

%% Compute cdkl5 production rates based on experimental data
final_cdkl5_mg = 5.2; % mg/l;
hours = 8;
biomass_ratio_wt_cdkl5 = Growth_curves.data(4,3)/Growth_curves.data(4,2);

% conditions are : 5 mM IPTG, 8hrs post_induction: 2.03 µg in 1 OD  ? 5.2 mg in 1 L culture
molecular_weight_cdkl5 = 128082.77 % mg/mml-1
OD_to_gL_scaling_factor = .74;
final_biomass_wt = Growth_curves.data(4,2)*OD_to_gL_scaling_factor;
final_biomass_cdkl5 = Growth_curves.data(4,3)*OD_to_gL_scaling_factor;
cdkl5_production_rate = final_cdkl5_mg/8/molecular_weight_cdkl5/final_biomass_cdkl5

%% Constrain the model with experimental data
[constrained_model, Rec_constrained_model,sol_wt_constrained_mu, sol_rec_constrained_mu, sol_rec_constrained_cdkl5 ]= Constraint_model_with_experiments(optBIOMASS_TACpB40_CDKL5, RxnExchangeGlcGln, Glutamate_UR,Gluconate_UR, biomass_ratio_wt_cdkl5);
%constrained_model = changeRxnBounds(constrained_model, 'rxn09240', 0, 'l'); 

% experimental growth rate of wt strain
mu_sim_wt = sol_wt_constrained_mu;
% experimental growth rate of rec strain
mu_sim_rec = sol_rec_constrained_mu;
% experimental growth CDKL5 prod rate 
cdkl5_sim_rec = sol_rec_constrained_cdkl5;

sol_rec_constrained = optimizeCbModel(Rec_constrained_model);
sol_wt_constrained = optimizeCbModel(constrained_model);

% experimental growth rate of wt strain
mu_exp_wt = (log(Growth_curves.data(4,2)) - log(Growth_curves.data(2,2)))/hours;
% experimental growth CDKL5 prod rate 
mu_exp_cdkl5 = (log(Growth_curves.data(4,3)) - log(Growth_curves.data(2,3)))/hours;


%% Plot the comparison between model predictions and experimental data
figure('Name','Figure1')
subplot(1,4,1)
plot(Growth_curves.data(:,1), Growth_curves.data(:,2), 'Marker', 'O', 'Color', 'k', 'LineWidth', 2)
hold on
plot(Growth_curves.data(:,1), Growth_curves.data(:,3), 'Marker', 'O', 'Color', '[0.4940 0.1840 0.5560]', 'LineWidth', 2)
legend('WT','hCDKL5', 'Location','northwest','Box','off')
set(gca,'FontName','Arial','fontsize',16)
title('A', 'FontSize', 15);
ylabel('OD (600nm)')
xlabel('Time since induction (hours)')


ModelGR = [mu_sim_wt, mu_sim_rec]
ExperimentalGR = [mu_exp_wt, mu_exp_cdkl5]
Combined = [ModelGR(:), ExperimentalGR(:)]
subplot(1,4,2)
bar1 = bar(Combined, 'LineWidth',1.5)
bar1(1).FaceColor = [0.4940 0.1840 0.5560];
bar1(2).FaceColor = [0.4660 0.6740 0.1880];
ylabel('\mu (h^-^1)')
xticklabels({'WT', 'Recombinant'})
set(gca,'FontName','Arial','fontsize',16)
legend('Model','Experimental', 'Location','northeast','Box','off')
title('B', 'FontSize', 15);
hold on


model_cdkl5 = cdkl5_sim_rec;
exp_cdkl5= cdkl5_production_rate;
Combined_cdkl5 = [model_cdkl5  cdkl5_production_rate];
subplot(1,4,3);
bar2 = bar(Combined_cdkl5', 'FaceColor','flat','LineWidth',1.5);
bar2(1).CData(2,:)=[0.4660 0.6740 0.1880];
bar2(1).CData(1,:)=[0.4940 0.1840 0.5560];

ylabel('mmmol/g_{CDW}*h^{-1}');
xticklabels({'hCDKL5 flux'});
set(gca,'FontName','Arial','fontsize',16);
title('C ', 'FontSize', 15);

subplot(1,4,4);
productionEnvelope(constrained_model, '', 'k', 'EX_cdkl5[c]',  'RXNbiomass', 'geneDelFlag', 10);
ylabel('hCDKL5 flux (mmmol/g_{CDW}*h^{-1})');
xlabel('\mu (h^{-1})');
set(gca,'FontName','Arial','fontsize',16);
title('D', 'FontSize', 15);


%% compare the flux distributions between the WT and the recombinant strains

[tableFVA, analysed_table, Intersect_important_reactions,WT_only, Rec_only]= CompareFluxes_WT_vs_Rec(constrained_model)    

sol_rec = optimizeCbModel(Rec_constrained_model)
sol_rec.x(findRxnIDs(Rec_constrained_model, Intersect_important_reactions))

sol_WT = optimizeCbModel(constrained_model)
sol_WT.x(findRxnIDs(constrained_model, Intersect_important_reactions))
constrained_model.rxnNames(findRxnIDs(constrained_model, Intersect_important_reactions))

%construct a table with fluxes for important reactions shared by the two
%models
table_flux_comparison = table(constrained_model.rxns(findRxnIDs(constrained_model, Intersect_important_reactions)), constrained_model.rxnNames(findRxnIDs(constrained_model, Intersect_important_reactions)), sol_rec.x(findRxnIDs(Rec_constrained_model, Intersect_important_reactions)), sol_WT.x(findRxnIDs(constrained_model, Intersect_important_reactions)))

table_flux_ratio_WT_rec = table(constrained_model.rxnNames(findRxnIDs(constrained_model, Intersect_important_reactions)), Intersect_important_reactions, ((sol_WT.x(findRxnIDs(constrained_model, Intersect_important_reactions))+0.000000001)./(sol_rec.x(findRxnIDs(Rec_constrained_model, Intersect_important_reactions))+0.000000001)))
increased_flux_in_recomb = table(table_flux_ratio_WT_rec.Intersect_important_reactions(table_flux_ratio_WT_rec.Var3<1), table_flux_ratio_WT_rec.Var1(table_flux_ratio_WT_rec.Var3<1))

%genes_in_increased_flux_in_recomb = findGenesFromRxns(constrained_model, table_flux_ratio_WT_rec.Intersect_important_reactions(table_flux_ratio_WT_rec.Var3<1));

increased_flux_in_recomb = table(table_flux_ratio_WT_rec.Intersect_important_reactions(table_flux_ratio_WT_rec.Var3<1), table_flux_ratio_WT_rec.Var1(table_flux_ratio_WT_rec.Var3<1))

fprintf('\nThe wt and recomb model share a set of  %d  important reactions \n\n', length(Intersect_important_reactions));

fprintf('\nThe recomb model shows an increased flux in the following %d reactions, compared to the  wt model \n\n', length(increased_flux_in_recomb.Var1));

writetable(increased_flux_in_recomb,'increased_flux_in_recomb.csv','Delimiter','	','QuoteStrings',false)

char(constrained_model.rxns(findRxnIDs(constrained_model, Rec_only)))


constrained_model.rxns(findRxnIDs(constrained_model, Rec_only))
constrained_model.rxnNames(findRxnIDs(constrained_model, Rec_only))
%writetable(char(constrained_model.rxns(findRxnIDs(constrained_model, Rec_only))),'active_only_in_recomb.csv','Delimiter','	','QuoteStrings',false)

%writecell(WT_only,'WT.csv','Delimiter','	','QuoteStrings',false)


% table(table_flux_ratio_WT_rec.Var1(table_flux_ratio_WT_rec.Var3>2 | table_flux_ratio_WT_rec.Var3<-2), table_flux_ratio_WT_rec.Intersect_important_reactions(table_flux_ratio_WT_rec.Var3>2 | table_flux_ratio_WT_rec.Var3<-2), table_flux_ratio_WT_rec.Var3(table_flux_ratio_WT_rec.Var3>2 | table_flux_ratio_WT_rec.Var3<-2))
% 
% %construct a table with important reactions for the WT
% important_only_in_WT = table( constrained_model.rxnNames(findRxnIDs(constrained_model, WT_only)), WT_only)
% important_only_in_rec = table( constrained_model.rxnNames(findRxnIDs(constrained_model, Rec_only)), Rec_only)
% 

%% search for overproduction strategies with try alternative carbon sources
    
[ratio_between_media, tested_mets, new_sol, model_EX_reactions_descr] = BetterMedia(constrained_model)

    % plot results on a log10  scale
    new_medium_res_table =table(ratio_between_media', new_sol', model_EX_reactions_descr,tested_mets, log10(ratio_between_media') );
    new_medium_res_table_sorted = sortrows(new_medium_res_table, 'Var1', 'descend')
    new_medium_res_table_sorted_log_GT1 = new_medium_res_table_sorted(new_medium_res_table_sorted.Var5 > 0,:);
    new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr =  strrep(new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr, 'EX_', '')
    new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr =  strrep(new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr, '_e', '')
    new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr =  strrep(new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr, '_', '')
    new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr =  strrep(new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr, 'EX', '')

    valueset = new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr(1:30)

    X = categorical(new_medium_res_table_sorted_log_GT1.model_EX_reactions_descr(1:30), valueset);
    figure('Name','Figure2')
    bar_media = bar(X,  new_medium_res_table_sorted_log_GT1.Var5(1:30))
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',18)
    xtickangle(45)
    ylabel('Log FC of CDKL5 production flux')
    bar_media(1).FaceColor = [0.4940 0.1840 0.5560];
    hold on
    

%% compare fluxes between GG medium and amended medium
%amylotriose
new_nutrient =     'EX_cpd01262_e';
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_amylotriose, Rec_GG_only_BetterMedium, increasing_flux_reactions_amylotriose , activating_reactions_amylotriose] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient);

%maltose
new_nutrient =     'EX_cpd00179_e';
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_maltose, Rec_GG_only_BetterMedium, increasing_flux_reactions_maltose , activating_reactions_maltose] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient)    

table_increasing_amylotriose_maltose = table(intersect(increasing_flux_reactions_amylotriose, increasing_flux_reactions_maltose), constrained_model.rxnNames(findRxnIDs (constrained_model, intersect(increasing_flux_reactions_amylotriose, increasing_flux_reactions_maltose))))

writecell(table_increasing_amylotriose_maltose.Var1,'increased_flux_amylotriose_maltose.csv','Delimiter','	','QuoteStrings',false)


%mannitol
new_nutrient =     'EX_cpd00314_e' ;
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_mannitole, Rec_GG_only_BetterMedium, increasing_flux_reactions_mannitol , activating_reactions_mannitol] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient)

writecell(increasing_flux_reactions_mannitol,'increased_flux_mannitol.csv','Delimiter','	','QuoteStrings',false)

increasing_flux_amylotriose_maltose_mannitol = intersect(intersect(increasing_flux_reactions_amylotriose, increasing_flux_reactions_maltose), increasing_flux_reactions_mannitol)

writecell(increasing_flux_amylotriose_maltose_mannitol,'increased_flux_amylotriose_maltose_mannitol.csv','Delimiter','	','QuoteStrings',false)

%thymidine
new_nutrient =     'EX_cpd00184_e' ;
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_mannitole, Rec_GG_only_BetterMedium, increasing_flux_reactions_thymidine , activating_reactions_thymidine] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient)

writecell(increasing_flux_reactions_thymidine,'increased_flux_thymidine.csv','Delimiter','	','QuoteStrings',false)

%galactose
new_nutrient =     'EX_cpd00108_e' ;
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_mannitole, Rec_GG_only_BetterMedium, increasing_flux_reactions_galactose , activating_reactions_galactose] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient)

writecell(increasing_flux_reactions_galactose,'increased_flux_galactose.csv','Delimiter','	','QuoteStrings',false)

%leucine
new_nutrient =     'EX_cpd00107_e' ;
[tableFVA_BetterMedium, analysed_table_BetterMedium, Intersect_important_reactions_BetterMedium, WT_new_nutrient_only_BetterMedium_mannitole, Rec_GG_only_BetterMedium, increasing_flux_reactions_leucine , activating_reactions_leucine] = ...
    CompareFluxes_WT_vs_Rec_better_medium(constrained_model, new_nutrient)

writecell(increasing_flux_reactions_leucine,'increased_flux_leucine.csv','Delimiter','	','QuoteStrings',false)


%% identification of putative targets for optimization of CDKL5 production - overexpression -

% the following code is run to include the .rev field in the constrained
% model 
raven_model = importModel('iMF721_SMBL3FBC2_v2.xml');
printRxnFormula(constrained_model, constrained_model.rxns);
rev_vector = [raven_model.rev(1:1325);0;0;0;0];
length(raven_model.rev);
length(constrained_model.rxns);
constrained_model.rev = rev_vector;
 
%raven_model = importModel('iMF721_SMBL3FBC2_v2.xml')
%table(constrained_model.rxns(1:1329),  raven_model.rxns(1:1329), raven_model.rev(1:1329), printRxnFormula(constrained_model, constrained_model.rxns(1:1329)))
targets=FSEOF(constrained_model,'RXNbiomass','EX_cdkl5[c]',100, 0.9, 'FSEOF.txt');
FSEOF_table = readtable('FSEOF.txt');
FSEOF_table_sorted = sortrows(FSEOF_table, 'Slope', 'descend');


constrained_model_no_ATPsulfate = changeRxnBounds(constrained_model, 'rxn00379', 0, 'b'); 
constrained_model_no_ATPsulfate = changeRxnBounds(constrained_model_no_ATPsulfate, 'rxn09240', 0, 'l'); 
optimizeCbModel(constrained_model_no_ATPsulfate)
optimizeCbModel(constrained_model)

%%%%%%%

raven_model = importModel('iMF721_SMBL3FBC2_v2.xml');
printRxnFormula(constrained_model_no_ATPsulfate, constrained_model_no_ATPsulfate.rxns);
rev_vector = [raven_model.rev(1:1325);0;0;0;0];
length(raven_model.rev);
length(constrained_model_no_ATPsulfate.rxns);
constrained_model_no_ATPsulfate.rev = rev_vector;
 
%raven_model = importModel('iMF721_SMBL3FBC2_v2.xml')
%table(constrained_model.rxns(1:1329),  raven_model.rxns(1:1329), raven_model.rev(1:1329), printRxnFormula(constrained_model, constrained_model.rxns(1:1329)))
targets=FSEOF(constrained_model_no_ATPsulfate,'RXNbiomass','EX_cdkl5[c]',100, 0.9, 'FSEOF.txt');
FSEOF_table = readtable('FSEOF.txt');
FSEOF_table_sorted = sortrows(FSEOF_table, 'Slope', 'descend');


%find the reactions that are common to FSEOF and flux change analysis.
intersect(FSEOF_table_sorted.EnzymeID, increased_flux_in_recomb.Var1)

FSEOF_subsystems = readtable('FSEOF_SubSystems/CountsOfInvolvedSubSys-FSEOF.txt');
FSEOF_subsystems_sorted = sortrows(FSEOF_subsystems, 'Var2', 'ascend')


%get only pathways with 2 or more reactions and remove no info categories
MoreThan2 = FSEOF_subsystems_sorted(FSEOF_subsystems_sorted.Var2>1,:);
MoreThan2 = MoreThan2(~(strcmp('No Info',MoreThan2.Var1)), :);
MoreThan2 = MoreThan2(~(strcmp('No Subsystem found',MoreThan2.Var1)), :);

figure('Name','Figure3');
bar_FSEOF= barh( MoreThan2.Var2);
set(gca,'XTick',1:length(MoreThan2.Var2));
set(gca,'YTickLabel',MoreThan2.Var1);

xlabel('N. of reactions');
set(gca,'FontName','Arial','fontsize',16);
bar_FSEOF(1).FaceColor = [0.4940 0.1840 0.5560];



