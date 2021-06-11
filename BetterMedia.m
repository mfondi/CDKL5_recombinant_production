function [ratio, model_EX_reactions, new_sol, model_EX_reactions_descr] = BetterMedia(model, nc)
%BetterMedia searched for suitable better media for CDKL5 overproduction by
%scanning all the model metabolites in order to find better combination of
%nutrients to improve the biosynthesis of the compound of interest
%   "model" is the input model and "nc" is the number of compouds that have to
%   be added to the medium in order to generate a new medium

model_better_media = changeObjective(model,'RXNbiomass');
biomass_constrained_sol_new_media = optimizeCbModel(model_better_media, 'max')
biomass_constraint_for_CDKL5_production_new_media = biomass_constrained_sol_new_media.f*0.74;
CDKL5model_constrained_exp_new_media = changeRxnBounds(model_better_media, 'RXNbiomass', biomass_constraint_for_CDKL5_production_new_media  , 'b' );
model_better_media = changeObjective(CDKL5model_constrained_exp_new_media,'EX_cdkl5[c]');
optimizeCbModel(model_better_media, 'max')
sol_original_medium = optimizeCbModel(model_better_media);
sol_original_medium_out =  sol_original_medium.f
ratio = [];
tested_mets = [];
new_sol = [];
model_EX_reactions = model_better_media.rxns(~cellfun(@isempty, regexp(model_better_media.rxns,'^EX_')));
model_EX_reactions_descr = model_better_media.rxnNames(~cellfun(@isempty, regexp(model_better_media.rxns,'^EX_')));
has_plasmid = cellfun(@isempty, regexp(model_EX_reactions, 'plasmid'));
model_EX_reactions = model_EX_reactions(has_plasmid);
model_EX_reactions_descr = model_EX_reactions_descr(has_plasmid);

    for i=1:1:length(model_EX_reactions)
    %TF = contains(model_EX_reactions{i}, 'plasmid')
    %if (TF ==0)
    model_in_new_medium = changeRxnBounds(model_better_media, model_EX_reactions{i}, -0.5 , 'l');
    sol_in_new_medium = optimizeCbModel(model_in_new_medium);
    sol_in_new_medium_out = sol_in_new_medium.f;
    ratio(end+1) = sol_in_new_medium_out/sol_original_medium_out;
    new_sol(end+1) = sol_in_new_medium_out;
    %tested_mets(end+1) = findRxnIDs(model_better_media,   model_EX_reactions{i});
    %else
    end


end

