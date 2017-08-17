% odes for PopZ model

function dydt = mrna_equations(t,y)

global p;

%% Equations

dydt=zeros(701,1);

% equation for cell growth
dydt(701) = p.mu*y(701);

% Equations for PopZ monomer

dydt(1) = p.popz_syn*y(201) - p.popz_mdeg*y(1) - (p.dnv + (1+p.prA*(y(301)+y(401)))*(p.aut*y(101)*y(101)))*y(1) + p.depol*y(101) + p.D_popZm*(y(2)-y(1))/(y(701)^2) - p.mu*y(1); 
for m = 2:99
dydt(m) = p.popz_syn*y(m+200) - p.popz_mdeg*y(m) - (p.dnv + (1+p.prA*(y(m+300)+y(m+400)))*(p.aut*y(m+100)*y(m+100)))*y(m) + p.depol*y(m+100) + p.D_popZm*(y(m-1)-2*y(m)+y(m+1))/(y(701)^2) - p.mu*y(m);
end
dydt(100) = p.popz_syn*y(300) - p.popz_mdeg*y(100) -(p.dnv + (1+p.prA*(y(400)+y(500)))*(p.aut*y(200)*y(200)))*y(100) + p.depol*y(200) + p.D_popZm*(y(99)-y(100))/(y(701)^2) - p.mu*y(100);
 
% Equations for PopZ polymer

dydt(101) = -p.popz_pdeg*y(101) + (p.dnv + (1+p.prA*(y(301)+y(401)))*(p.aut*y(101)*y(101)))*y(1) - p.depol*y(101) + p.D_popZp*(y(102)-y(101))/(y(701)^2) - p.mu*y(101) ;
for pl = 102:199
dydt(pl) = -p.popz_pdeg*y(pl) + (p.dnv + (1+p.prA*(y(pl+200)+y(pl+300)))*(p.aut*y(pl)*y(pl)))*y(pl-100) - p.depol*y(pl) + p.D_popZp*(y(pl+1)-2*y(pl)+y(pl-1))/(y(701)^2) - p.mu*y(pl);
end
dydt(200) = -p.popz_pdeg*y(200) + (p.dnv + (1+p.prA*(y(400)+y(500)))*(p.aut*y(200)*y(200)))*y(100) - p.depol*y(200) + p.D_popZp*(y(199)-y(200))/(y(701)^2) - p.mu*y(200) ;

% Equations for popZ mRNA

% Compute position of popZ genes in cytoplasm
%gene_st_pos = ceil(56*1.3/(100*y(701)));
%gene_sw_pos = 100 - gene_st_pos;
 gene_st_pos = 10;
 gene_sw_pos = 90;

% gene_pos = [270,299];%[202,gene_st_pos+200,gene_sw_pos+200, 275 299];
dydt(201) =  - p.kdeg_mrna*y(201) + p.D_popZ_mr*(y(202)-y(201))/(y(701)^2) - p.mu*y(201);
for mr = 202:299
    %if (~isempty(find(gene_pos == mr))) 
    if mr == (gene_sw_pos + 200)
        dydt(mr) = p.ksyn_sw_mrna - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(701)^2) - p.mu*y(mr);
      elseif mr == (gene_st_pos + 200)
          dydt(mr) = p.ksyn_st_mrna - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(701)^2) - p.mu*y(mr);
     else dydt(mr) =  - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(701)^2) - p.mu*y(mr);
    end
end
dydt(300) =  - p.kdeg_mrna*y(300) + p.D_popZ_mr*(y(299)-y(300))/(y(701)^2) - p.mu*y(300);

% ParA:PopZ
for prab_idx = 301:400
    dydt(prab_idx) =  p.kf_pa_pz*y(prab_idx-200)*y(prab_idx+200)-  p.kr_pa_pz*y(prab_idx);
end

% ParA:TipN
for praz_idx = 401:500
    dydt(praz_idx) = p.kf_pa_tn*y(praz_idx+200)*y(praz_idx+100) -  p.kr_pa_tn*y(praz_idx);
end

% Free ParA
for praf_idx = 501:600
    dydt(praf_idx) = - p.kf_pa_pz*y(praf_idx-400)*y(praf_idx) +  p.kr_pa_pz*y(praf_idx-200) - p.kf_pa_tn*y(praf_idx+100)*y(praf_idx) +  p.kr_pa_pz*y(praf_idx-100);
end

% TipN proxy
dydt(601:700) = 0;
