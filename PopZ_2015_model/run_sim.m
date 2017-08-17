% Main file for PopZ 
clear M PopZ mRNA gene
global p
load('PopZ_polymer_init.mat');
load('polymer_init_round.mat');
load('ovex_init_full.mat');
parameters(1);
p.ksyn_st_mrna = 0;

% Initial conditions
y0 = zeros(301,1);
%y0 = popz_init.init.';

y0(101:200) = polymer_init_round.';
y0(1:100) = 0;
y0(201:300) = 0;
y0(301) = 0.013;

% Integration parameters
t0  =  0;       % Start time
tf  = 150;     % End time

options  =  odeset('Events',@popz_event,'RelTol',1e-4,'AbsTol',1e-6);
tout = t0;
y0 = y0.';
yout = y0;
teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
    
    [t,y,te,ye,ie] = ode15s(@mrna_equations,[t0 tf],y0,options);
    
    nt = length(t);
    
    tout = [tout;t(2:nt)];
    yout = [yout;y(2:nt,:)];
    teout  =  [teout;te];
    yeout  =  [yeout;ye];
    ieout  =  [ieout;ie];
    
    y0  =  y(nt,:);
    
    if isscalar(ie)  ==  0
        ie  =  0;
    end
    
    if ie  ==  1
       p.ksyn_st_mrna = p.ksyn_sw_mrna;
       %p.ksyn_sw_mrna = 0;
    end
    
    t0 = t(nt);
    
    if t0 >= tf
        break;
    end
    ie
end
%% Generating a space-time plot

% In the plot, the X-axis refers to space and Y-axis is time. Since the
% cell grows in time, the X-axis co-odinates in terms of um  for each of 
% the 100 frid points changes at every time step. We therefore define the
% X-axis by a matrix of 100 columns, where each row contains the X-axis
% co-ordinates in um for a given time step

for n = 51:100
    M(:,n) = yout(:,301)*(n/100)-0.5*(yout(:,301)*(.5) + yout(:,301)*.51); % Each element from column 51 to 100 is assigned the 'n'th fraction of total cell length at a given time step followed by subtracting the mid point value so that the centre of the grid takes the value 0
end

M(:,1:50) = fliplr(M(:,51:100));    % The 1 st 50 colums is the flipped version of the next 50 eg: 3 2 1 1 2 3
M(:,1:50) = -M(:,1:50);             % now its -3 -2 -1 1 2 3
M = 100*M;

PopZ(:,1:100) = yout(:,101:200) + yout(:,1:100);
mRNA(:,1:100) = yout(:,201:300);

PopZ = fliplr(PopZ);
mRNA = fliplr(mRNA);


PopZ = PopZ.';
mRNA = mRNA.';
M = M.';

% Record gene position
gene = zeros(size(mRNA));

for time_idx = 1:length(tout)
    gene_st_pos = ceil(100*0.56*0.013/yout(time_idx,301));
    gene_sw_pos = 100 - gene_st_pos;
    gene(gene_sw_pos, time_idx) = 1;
    if tout(time_idx) >= 50
    gene(gene_st_pos, time_idx) = 1;
    end
end

gene = flipud(gene);

% syntax to genearte space-timeplot
figure(1)

ax1 = subplot(3,1,1);
pcolor(tout,M, gene)
shading interp
colormap(jet)
caxis([0  1])
ylabel('popZ gene')

ax2 = subplot(3,1,2);
pcolor(tout,M, mRNA)
shading interp
colormap(jet)
caxis([0  0.25])
ylabel('popZ mRNA')

ax3 = subplot(3,1,3);
pcolor(tout,M, PopZ)
shading interp
colormap(jet)
%caxis([0  20])
ylabel('PopZ polymer')
xlabel('time (min)')

