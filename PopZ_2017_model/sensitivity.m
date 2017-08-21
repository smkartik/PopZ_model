% Sensitivity range for parameters

synthesis_range = (0:0.5:10);
dnv_range = (0:10:300);
aut_range = (0:0.5:10);
depol_range = (0:0.05:0.5);

input_range = synthesis_range;


PopZ_end_dist = zeros(100, length(input_range));
T_bipolar = zeros(1, length(input_range));

parfor idx = 1:length(input_range)
    output = run_sim(input_range(idx));
    PopZ_end_dist(:,idx) = output.PopZ(:,end);
    time_index = find(sum(output.PopZ(81:100,:),1) > 0.2*sum(output.PopZ,1),1, 'first');
    if isempty(time_index)
        T_bipolar(1, idx) = 0;
    else
    T_bipolar(1, idx) = output.time(time_index,1);
    end
end

y_range = 1:1:100;
input_label = 'autocatalytic rc';

figure(1)
hFig = figure(1);
xwidth = 200;
ywidth = 250;

set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 xwidth ywidth])
ax1 = subplot(2,1,1);
pcolor(input_range,y_range,PopZ_end_dist)
shading interp
colormap('jet')
%caxis([0 60])
%colorbar
%ylabel('distribution at t = 150 min')
xlabel(input_label)

ax2 = subplot(2,1,2);
plot(input_range, T_bipolar, 'ro') 
meanline = line([2.5, 2.5], ylim)
%ylabel('time when bipolar')
xlabel(input_label)
print('aut_sens_prna_null_mrna_slow', '-dpng', '-r600')

