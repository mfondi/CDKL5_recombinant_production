
    function [Glutamate_UR,Gluconate_UR] = Compute_GG_consumption()
%Copmute uptake rates and other growth parameters from experimental (NMR)
%data. Outputs are the utilization rates to be used to constrain the model
%(lowe bounds) when replicating the experimental conditions (growth in GG
%plus Shatz salts


    GG_consumption= importdata('GG_consumption.txt')

    % Plot Glutamate and Gluconate consumption

    figure('Name','SX')
    subplot(1,2,2)
    errorbar(GG_consumption.data(:,1), GG_consumption.data(:,2),GG_consumption.data(:,3) ,'Marker', 'O', 'Color', 'k', 'LineWidth', 2)
    hold on
    errorbar(GG_consumption.data(:,1), GG_consumption.data(:,4),GG_consumption.data(:,5), 'Marker', 'O', 'Color', 'r', 'LineWidth', 2) 
    %ylim([.14,.2]);
    %xlim([0,10]);
    legend('Glutamate','Gluconate', 'Location','southwest','Box','off')
    set(gca,'FontName','Arial','fontsize',16)
    ylabel('Concentration (mM)')
    xlabel('Time (hours)')
    title('B')
    hold on
    % Plot Biomass
    subplot(1,2,1)
    errorbar(GG_consumption.data(:,1), GG_consumption.data(:,6),GG_consumption.data(:,7) ,'Marker', 'O', 'Color', 'k', 'LineWidth', 2)
    set(gca,'FontName','Arial','fontsize',16)
    ylabel('OD (600 nM)')
    xlabel('Time (hours)')
	title('A')

    %% Compute Uptake rates
    OD_to_gL_scaling_factor = .74;
    hours = 8;
    volume = 1.6 %l
    %consumed mmol for glutamate (MW = 147.13 g/mol - initial concentration = 5g/l - volume = 1.6 l)
    Glutamate_consumed_moles =GG_consumption.data(1,2) -  GG_consumption.data(end,2);
    %consumed mmol for gluconate (MW = 196.16 g/mol - initial concentration = 5g/l - volume = 1.6 l)
    Gluconate_consumed_moles =GG_consumption.data(1,4) -  GG_consumption.data(end,4);
    %overall biomasss (final OD * scaling factor * total volume (1.6 L)
    biomass = GG_consumption.data(end,6)*OD_to_gL_scaling_factor*volume;
    %Growth rate (ln(OD_final) - ln(OD_init)/hours)
    mu = (log(GG_consumption.data(end,6)) - log(GG_consumption.data(1,6)))/hours;
    %Glutamate yield
    Glutamate_yield = biomass/Glutamate_consumed_moles;
    %Gluconate yield
    Gluconate_yield = biomass/Gluconate_consumed_moles;
    %Glutamate uptake rate (Growht rate / Yield)
    Glutamate_UR = mu/Glutamate_yield;
    %Gluconate uptake rate
    Gluconate_UR = mu/Gluconate_yield;

    
    
    
    end



