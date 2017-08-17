function parameters(cond)

% Global variables
global p
 
% Growth conditions
if cond == 0              % cell size is fixed
    p.mu = 0;             % growth rate constant    
elseif cond == 1          % cell is growing 
    p.mu = 0.0055;        % growth rate constant = 0.0055 (units => 1/min)    
end

% Diffusion parameters (units => um^2/min)
p.D_popZm = 100;         % PopZ monomer dffusion  
p.D_popZp = 0.0005;%5;        % PopZ polymer diffusion
p.D_popZ_mr = 1;%10;%0.05;%10;%0.05;         % popZ mRNA diffusion  %0.5 

% Reaction rates (units => 1/min) 
p.popz_syn = 5;%10;       % popZ monomer synthesis
p.popz_mdeg = 0.05;     % popZ monomer degradation
p.popz_pdeg = 0.05;     % popZ polymer degradation
p.dnv = 0;%35;            % denovo polymerization
p.aut = 2.5;%3;             % autocatalytic polymerization
p.depol = 0.1;%0.05;         % deplymerization 
p.ksyn_sw_mrna = 3;    % synthesis of mRNA from gene in swarmer half
p.kdeg_mrna = 0.5;


% % Diffusion parameters (units => um^2/min)
% p.D_popZm = 100;         % PopZ monomer dffusion  
% p.D_popZp = 0.005;        % PopZ polymer diffusion
% p.D_popZ_mr = 0.05;         % popZ mRNA diffusion  %0.5 
% 
% % Reaction rates (units => 1/min) 
% p.popz_syn = 30;       % popZ monomer synthesis
% p.popz_mdeg = 0.5;     % popZ monomer degradation
% p.popz_pdeg = 0.5;     % popZ polymer degradation
% p.dnv = 30;            % denovo polymerization
% p.aut = 3;             % autocatalytic polymerization
% p.depol = 0.1;         % deplymerization 
% p.ksyn_sw_mrna = 3;    % synthesis of mRNA from gene in swarmer half
% p.kdeg_mrna = 0.5;
