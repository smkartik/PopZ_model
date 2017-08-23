% odes for PopZ model

function dydt = caulobacter_model_equations(t,y)

global p;
   
%% Equations

dydt=zeros(901,1);

% equation for cell growth
dydt(901) = p.mu*y(901);

%% Equations for PopZ monomer
% Bin 1
dydt(1) = p.popz_syn*y(201) - p.popz_mdeg*y(1)... % synthesis & degradation 
    - p.dnv*y(1)...% denovo polymerization
    - p.aut1*y(101)*y(101)*y(1)...% autocatalytic polymerization
    - p.prA*y(301)*(p.aut2*y(101)*y(101))*y(1)...% ParA-dependent autocatalytic polymerization
    - p.sxp*y(501)*(p.aut2*y(101)*y(101))*y(1)...% SpmX-dependent autocatalytic polymerization
    + p.depol*y(101)...% PopZ depolymerization
    + p.D_popZm*(y(2)-y(1))/(y(901)^2) - p.mu*y(1); % monomer diffusion % dilution
% Bin 2 to 99
for m = 2:99
dydt(m) = p.popz_syn*y(m+200) - p.popz_mdeg*y(m)...% synthesis & degradation
    - p.dnv*y(m)...% denovo polymerization
    - p.aut1*y(m+100)*y(m+100)*y(m)...% autocatalytic polymerization
    - p.prA*y(m+300)* (p.aut2*y(m+100)*y(m+100))*y(m)...% ParA-dependent autocatalytic polymerization
    - p.sxp*y(m+500)*(p.aut2*y(m+100)*y(m+100))*y(m)...% SpmX-dependent autocatalytic polymerization
    + p.depol*y(m+100)...% PopZ depolymerization
    + p.D_popZm*(y(m-1)-2*y(m)+y(m+1))/(y(901)^2) - p.mu*y(m);  % monomer diffusion % dilution
end
% Bin 100
dydt(100) = p.popz_syn*y(300) - p.popz_mdeg*y(100)...
    - p.dnv*y(100)...% denovo polymerization
    - p.aut1*y(200)*y(200)*y(100)...% autocatalytic polymerization
    - p.prA*y(400)*(p.aut2*y(200)*y(200))*y(100)...% ParA-dependent autocatalytic polymerization
    - p.sxp*y(600)*(p.aut2*y(200)*y(200))*y(100)...% SpmX-dependent autocatalytic polymerization
    + p.depol*y(200)...% PopZ depolymerization
    + p.D_popZm*(y(99)-y(100))/(y(901)^2) - p.mu*y(100); % monomer diffusion % dilution
 
%% Equations for PopZ polymer
% Bin 1
dydt(101) = -p.popz_pdeg*y(101)...
    + p.dnv*y(1)...% denovo polymerization
    + p.aut1*y(101)*y(101)*y(1)...% autocatalytic polymerization
    + p.prA*y(301)*(p.aut2*y(101)*y(101))*y(1)...% ParA-dependent autocatalytic polymerization
    + p.sxp*y(501)*(p.aut2*y(101)*y(101))*y(1)...% SpmX-dependent autocatalytic polymerization
    - p.depol*y(101)...% PopZ depolymerization
    + p.D_popZp*(y(102)-y(101))/(y(901)^2) - p.mu*y(101); % polymer diffusion % dilution
% Bin 2 to 99
for pl = 102:199
dydt(pl) = -p.popz_pdeg*y(pl)...
    + p.dnv*y(pl-100)...% denovo polymerization
    + p.aut1*y(pl)*y(pl)*y(pl-100)...% autocatalytic polymerization
    + p.prA*y(pl+200)*(p.aut2*y(pl)*y(pl))*y(pl-100)...% ParA-dependent autocatalytic polymerization
    + p.sxp*y(pl+400)*(p.aut2*y(pl)*y(pl))*y(pl-100)...% SpmX-dependent autocatalytic polymerization
    - p.depol*y(pl)...% PopZ depolymerization
    + p.D_popZp*(y(pl+1)-2*y(pl)+y(pl-1))/(y(901)^2) - p.mu*y(pl); % polymer diffusion % dilution
end
% Bin 100
dydt(200) = -p.popz_pdeg*y(200)...
    + p.dnv*y(100)...% denovo polymerization
    + p.aut1*y(200)*y(200)*y(100)...% autocatalytic polymerization
    + p.prA*y(400)*(p.aut2*y(200)*y(200))*y(100)...% ParA-dependent autocatalytic polymerization
    + p.sxp*y(600)*(p.aut2*y(200)*y(200))*y(100)...% SpmX-dependent autocatalytic polymerization
    - p.depol*y(200)...% PopZ depolymerization
    + p.D_popZp*(y(199)-y(200))/(y(901)^2) - p.mu*y(200); % polymer diffusion % dilution

%% Equations for popZ mRNA
% Compute position of popZ genes in cytoplasm
gene_st_pos = ceil(56*1.3/(100*y(901)));
gene_sw_pos = 100 - gene_st_pos;

% Bin 1
dydt(201) =  - p.kdeg_mrna*y(201)...% PopZ mRNA degradation
    + p.D_popZ_mr*(y(202)-y(201))/(y(901)^2) - p.mu*y(201);% mRNA diffusion & dilution
% Bin 2 to 99
for mr = 202:299
    if mr == (gene_sw_pos + 200) % synthesize mRNA if popZ gene in bin, else only degrade and diffuse
        dydt(mr) = p.ksyn_sw_mrna...% synthesis of PopZ mRNA from gene in old chromosome
            - p.kdeg_mrna*y(mr)...% PopZ mRNA degradation
            + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(901)^2) - p.mu*y(mr); % mRNA diffusion & dilution
      elseif mr == (gene_st_pos + 200)
          dydt(mr) = p.ksyn_st_mrna ...% synthesis of PopZ mRNA from gene in new chromosome
              - p.kdeg_mrna*y(mr)...% PopZ mRNA degradation
              + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(901)^2) - p.mu*y(mr);% mRNA diffusion & dilution
     else dydt(mr) =  - p.kdeg_mrna*y(mr)...% mRNA degradation 
             + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(901)^2) - p.mu*y(mr);% mRNA diffusion & dilution
    end
end
% Bin 100
dydt(300) =  - p.kdeg_mrna*y(300)...% PopZ mRNA degradation
    + p.D_popZ_mr*(y(299)-y(300))/(y(901)^2) - p.mu*y(300);% mRNA diffusion & dilution

%% Equation for ParA
dydt(301:400) = 0;

%% Equation for SpmX free
dydt(401) = p.spmx_syn - p.spmx_deg*y(401)...
    - p.kf_sx_pz*y(401)*y(101)*y(601) + p.kr_sx_pz*y(501)...
    + p.D_spmx*(y(402) - y(401))/ (y(901)^2) - p.mu*y(401);
for spf_idx=402:499
    dydt(spf_idx) = p.spmx_syn - p.spmx_deg*y(spf_idx)...
    - p.kf_sx_pz*y(spf_idx)*y(spf_idx-300)*y(spf_idx+200) + p.kr_sx_pz*y(spf_idx+100)...
    + p.D_spmx*(y(spf_idx-1) - 2*y(spf_idx) + y(spf_idx+1))/ (y(901)^2) - p.mu*y(spf_idx);
end
dydt(500) = p.spmx_syn - p.spmx_deg*y(500)...
    - p.kf_sx_pz*y(500)*y(200)*y(700) + p.kr_sx_pz*y(600)...
    + p.D_spmx*(y(499) - y(500))/ (y(901)^2) - p.mu*y(500);

%% Equation for Spmx bound to poles
for spc_idx=501:600
    dydt(spc_idx) = p.kf_sx_pz*y(spc_idx-100)*y(spc_idx-400)*y(spc_idx+100)...
        - p.kr_sx_pz*y(spc_idx);
end

%% Polar regions
dydt(601:700) = 0;

%% dummy variables for extension
dydt(701:900) = 0;

%% cell growth equation
dydt(901) = p.mu*y(901);