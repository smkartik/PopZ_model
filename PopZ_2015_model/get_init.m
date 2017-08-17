% Main file for PopZ 
function popz_init = get_init() 
%clear M PopZ mRNA gene
load('polymer_init.mat'); 
global p

parameters(0);

% Initial conditions
y0 = zeros(301,1);


y0(1:100) = 0;                  % PopZ monomer
y0(101:200) = polymer_init.';      % PopZ polymer 
y0(301) = 0.013;                % initial size of bin => 0.01*size of swarmer cell (1.3 um) 
y0(101:200)
p.ksyn_st_mrna = 0; % synthesis from second gene set to zero for first 50 minutes

% Integration parameters

t0=0;
tf=100; % average time for complete cell cycel from new born swarmer to predivisonal cell is 150 mins

[t,y]=ode15s(@mrna_equations,[t0 tf],y0); % the ode solver 15s takes file containing equations, start and end time, inital conditions as arguments

popz_init.init = y(end,:);

end