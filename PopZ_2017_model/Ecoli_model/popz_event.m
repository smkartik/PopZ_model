function [value,isterminal,direction] = popz_event(t,y)

S_phase_start = 30;


value = [sign(t - S_phase_start)];
isterminal = [1];
direction = [+1];
end