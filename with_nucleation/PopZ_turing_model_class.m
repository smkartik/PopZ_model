classdef PopZ_turing_model
    
    % A class implimenting the odes for the PopZ turing model
    
    % EXAMPLE RUN
    % -----------
    %
    % >> m = PopZ_turing_model()
    % >> tspan = [0 150]
    % >> [t y] = ode15s(@m.odes, tspan, m.get_initial_vales());
    %
    % Retrieve the observables
    %
    % >> y_obs = m.get_observables(y)
    %
    properties
        observables
        parameters
    end
    
    methods
        function self = PopZ_turing_model()
            % Assign default parameter values
            self.parameters = struct(...
                'PopZ_polymer_init', polymer_init, ...
                'PopZ_monomer_init', 0, ...
                'PopZ_mRNA_init', 0, ...
                'kpopz_syn', 8, ... % units => per min
                'kpopz_mdeg', 0.05, ...
                'kpopz_pdeg', 0.05, ...
                'kdnv', 35, ...
                'kaut', 2.5, ...
                'kdepol', 0.1, ...
                'ksyn_sw_mrna', 3, ...
                'kdeg_mrna', 0.5, ...
                'kaut_prA', 1, ...
                'D_popZm', 100, ... % units => micrometer^2/min
                'D_popZp', 0.0005, ...
                'D_popZ_mr', 0.05);
            
            % Define species indices (first row) and coefficients
                
                
            