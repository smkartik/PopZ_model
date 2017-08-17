% odes for PopZ model

function dydt = mrna_equations(t,y)

global p;

%% Equations

dydt=zeros(301,1);

% equation for cell growth
dydt(301) = p.mu*y(301);

% Equations for PopZ monomer

dydt(1) = p.popz_syn*y(201) - p.popz_mdeg*y(1) - (p.dnv + (p.aut*y(101)*y(101)))*y(1) + p.depol*y(101) + p.D_popZm*(y(2)-y(1))/(y(301)^2) - p.mu*y(1); 
for m = 2:99
dydt(m) = p.popz_syn*y(m+200) - p.popz_mdeg*y(m) - (p.dnv + (p.aut*y(m+100)*y(m+100)))*y(m) + p.depol*y(m+100) + p.D_popZm*(y(m-1)-2*y(m)+y(m+1))/(y(301)^2) - p.mu*y(m);
end
dydt(100) = p.popz_syn*y(300) - p.popz_mdeg*y(100) -(p.dnv + (p.aut*y(200)*y(200)))*y(100) + p.depol*y(200) + p.D_popZm*(y(99)-y(100))/(y(301)^2) - p.mu*y(100);
 
% Equations for PopZ polymer

dydt(101) = -p.popz_pdeg*y(101) + (p.dnv + (p.aut*y(101)*y(101)))*y(1) - p.depol*y(101) + p.D_popZp*(y(102)-y(101))/(y(301)^2) - p.mu*y(101) ;
for pl = 102:199
dydt(pl) = -p.popz_pdeg*y(pl) + (p.dnv + (p.aut*y(pl)*y(pl)))*y(pl-100) - p.depol*y(pl) + p.D_popZp*(y(pl+1)-2*y(pl)+y(pl-1))/(y(301)^2) - p.mu*y(pl);
end
dydt(200) = -p.popz_pdeg*y(200) + (p.dnv + (p.aut*y(200)*y(200)))*y(100) - p.depol*y(200) + p.D_popZp*(y(199)-y(200))/(y(301)^2) - p.mu*y(200) ;

% Equations for popZ mRNA

% Compute position of popZ genes in cytoplasm
gene_st_pos = ceil(56*1.3/(100*y(301)));
gene_sw_pos = 100 - gene_st_pos;
 %gene_st_pos = 34;
 %gene_sw_pos = 68;

dydt(201) =  - p.kdeg_mrna*y(201) + p.D_popZ_mr*(y(202)-y(201))/(y(301)^2) - p.mu*y(201);
for mr = 202:299
    if mr == (gene_sw_pos + 200)
        dydt(mr) = p.ksyn_sw_mrna - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(301)^2) - p.mu*y(mr);
    elseif mr == (gene_st_pos + 200)
        dydt(mr) = p.ksyn_st_mrna - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(301)^2) - p.mu*y(mr);
    else dydt(mr) =  - p.kdeg_mrna*y(mr) + p.D_popZ_mr*(y(mr-1)-2*y(mr)+y(mr+1))/(y(301)^2) - p.mu*y(mr);
    end
end
dydt(300) =  - p.kdeg_mrna*y(300) + p.D_popZ_mr*(y(299)-y(300))/(y(301)^2) - p.mu*y(300);

