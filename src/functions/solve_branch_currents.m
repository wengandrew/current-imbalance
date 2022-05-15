function [Ia, Ib] = solve_branch_currents(I, dz, alpha, Ra, Rb)

    Ia = (+alpha*dz + I*Rb) / (Ra + Rb);
    Ib = (-alpha*dz + I*Ra) / (Ra + Rb);

end