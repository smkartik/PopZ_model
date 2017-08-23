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
p.D_popZ_mr = 100;         % popZ mRNA diffusion  %0.5 
p.D_parA = 100;
p.D_spmx = 1;

% Reaction rates (units => 1/min) 
p.popz_syn = range_value;   % popZ monomer synthesis
p.popz_mdeg = 0.1;          % popZ monomer degradation
p.popz_pdeg = 0.05;         % popZ polymer degradation
p.dnv = 60;                 % denovo polymerization
p.aut1 = 0.5;               % autocatalytic polymerization
p.aut2 = 2.5;               % autocatalytic polymerization
p.depol = 0.1;              % deplymerization 
p.ksyn_sw_mrna =3;          % synthesis of mRNA from gene in swarmer half
p.kdeg_mrna = 0.5;          % mRNA degradation
p.prA = 1;                  % 1 if ParA-dependent autocatalysis
p.spmx_syn = 0;             % SpmX synthesis
p.spmx_deg = 0.1;           % SpmX degradation
p.sxp = 10;                 % 10 if SpmX-dependent autocatalysis
p.kf_sx_pz = 1e-2;          % PopZ dependent localization of SpmX at pole
p.kr_sx_pz = 1e-2;          % SpmX dissocaites from polar scaffold

% ecoli_params
p.D_x= 100;
p.kfb = 1;
p.krb = 1e-5;

