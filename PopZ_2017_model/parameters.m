function parameters(cond, range_value)

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
p.D_popZp = 0.0005;        % PopZ polymer diffusion
p.D_popZ_mr = 1;         % popZ mRNA diffusion  %0.5 

% Reaction rates (units => 1/min) 
p.popz_syn = range_value;       % popZ monomer synthesis
p.popz_mdeg = 0.1; %0.1     % popZ monomer degradation
p.popz_pdeg = 0.05;     % popZ polymer degradation
p.dnv = 60;            % denovo polymerization
p.aut = 2.5;          % autocatalytic polymerization
p.depol = 0.1;         % deplymerization 
p.ksyn_sw_mrna =3;    % synthesis of mRNA from gene in swarmer half
p.kdeg_mrna = 0.5;
p.prA = 0;
p.kf_pa_tn = 1e-1;
p.kr_pa_tn = 1e-2;
p.kf_pa_pz = 1e-5;%1e-5;
p.kr_pa_pz = 1e-2;

p.spmx_syn = 0;
p.spmx_deg = 0.001;
p.sxp = 10;
p.kf_sx_pz = 1e-5;
p.kr_sx_pz = 1e-2;
p.D_spmx = 1;

p.prAdnv = 0;%10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Diffusion parameters (units => um^2/min)
% p.D_popZm = 100;         % PopZ monomer dffusion  
% p.D_popZp = 0.0005;        % PopZ polymer diffusion
% p.D_popZ_mr = 1;         % popZ mRNA diffusion  %0.5 
% 
% % Reaction rates (units => 1/min) 
% p.popz_syn = 5;       % popZ monomer synthesis
% p.popz_mdeg = 0.05; %0.1     % popZ monomer degradation
% p.popz_pdeg = 0.05;     % popZ polymer degradation
% p.dnv = 60;%35;            % denovo polymerization
% p.aut = 2.5;          % autocatalytic polymerization
% p.depol = 0.1;         % deplymerization 
% p.ksyn_sw_mrna =3;    % synthesis of mRNA from gene in swarmer half
% p.kdeg_mrna = 0.5;
% p.prA = 1;


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
