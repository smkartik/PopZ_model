function [value,isterminal,direction] = popz_event(t,y)

S_phase_start = 30;
Spmx_phase_end = 50;

value = [sign(t - S_phase_start); sign(t-Spmx_phase_end)];
isterminal = [1; 1];
direction = [+1; +1];
end