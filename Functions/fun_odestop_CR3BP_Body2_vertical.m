% poincare map where x = 1 - mu
function [value,isterminal,direction] = fun_odestop_CR3BP_Body2_vertical(t,x,mu)

  value = x(1) - (1-mu);
  isterminal = 0; % if isterminal == 1, the calcuration will stop when the event occurs.
  direction = 0; % if direction ==  0, all zeros are to be computed (the default)
                 % if direction == +1, only the zeros where the event function increases
                 % if direction == -1, only the zeros where the event function decreases.
end