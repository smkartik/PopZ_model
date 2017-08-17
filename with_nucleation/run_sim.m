
function output = run_sim(range_value)
% Main file for PopZ 
clear M PopZ mRNA gene ParA
global p
%load('PopZ_polymer_init.mat');
load('polymer_init_round.mat');
load('thick_init.mat');
parameters(1, range_value);
p.ksyn_st_mrna = 0;

% Initial conditions
y0 = zeros(701,1);


y0(101:200) = 2*polymer_init_round.';
%y0(191:200) = 1.5*y0(191:200);
y0(1:100) = 0;
y0(201:300) = 0;
y0(601:700) = 0.1;
y0(501:600) = 1; %ParA = 1/bin
y0(701) = 0.013;

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
        %y0(301:310) = 0.9;
        y0(611:700) = 0.0111;
        y0(601:610) = 0.9;
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
    M(:,n) = yout(:,701)*(n/100)-0.5*(yout(:,701)*(.5) + yout(:,701)*.51); % Each element from column 51 to 100 is assigned the 'n'th fraction of total cell length at a given time step followed by subtracting the mid point value so that the centre of the grid takes the value 0
end

M(:,1:50) = fliplr(M(:,51:100));    % The 1 st 50 colums is the flipped version of the next 50 eg: 3 2 1 1 2 3
M(:,1:50) = -M(:,1:50);             % now its -3 -2 -1 1 2 3
M = 100*M;

PopZ(:,1:100) = yout(:,101:200) + yout(:,1:100);
mRNA(:,1:100) = yout(:,201:300);
ParA(:,1:100) = yout(:,301:400)+yout(:,401:500);

PopZ = fliplr(PopZ);
mRNA = fliplr(mRNA);
ParA = fliplr(ParA);


PopZ = PopZ.';
mRNA = mRNA.';
ParA = ParA.';
M = M.';

output.PopZ = PopZ;
output.time = tout;
output.ParA = ParA;
output.mRNA = mRNA;
output.grid = M;
%Record gene position
gene = zeros(size(mRNA));

for time_idx = 1:length(tout)
    gene_st_pos = ceil(100*0.56*0.013/yout(time_idx,701));
    % Define cell boundaries
    gene(1, time_idx) = 1;
    gene(100, time_idx) = 1;
    
    % Deine gene position
    gene_sw_pos = 100 - gene_st_pos;
    gene(gene_sw_pos, time_idx) = 1;
    if tout(time_idx) >= 30
    gene(gene_st_pos, time_idx) = 1;
    gene(gene_sw_pos, time_idx) = 1;
    end
end

gene = flipud(gene);

%syntax to genearte space-timeplot
figure(1)
hFig = figure(1);
xwidth = 200;
ywidth = 400;

set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 xwidth ywidth])

ax1 = subplot(3,1,1);
pcolor(tout,M, gene)
shading interp
colormap(ax1, flipud(gray))
colorbar()
caxis([0  1])
title('popZ gene')

ax2 = subplot(3,1,2);
pcolor(tout,M, mRNA)
shading interp
colormap(ax2,'jet')
colorbar()
caxis([0  0.25])
title('popZ mRNA')
ylabel('cell size')

ax3 = subplot(3,1,3);
pcolor(tout,M, PopZ)
shading interp
colormap(ax3,'jet')
colorbar()
%caxis([0  20])
title('PopZ polymer')
xlabel('time (min)')
% % print('ovex_prna_null_mrna_fst', '-dpng', '-r600')

output.gene = gene
save ./springer_ouput/popZ 
