function parallel_ocv_r_symbolic()
    % Parallel OCV-R models
    % 03/16/2022

    syms I Ia Ib Ic Id       real positive
    syms Q Qa Qb Qc Qd       real positive
    syms R Ra Rb Rc Rd       real positive
    syms za zb zc zd         real positive
    syms x1 x2 x3 x4         real
    syms Ua Ub Uc Ud         real positive
    syms la lb lc ld lbd     real positive
    syms Ua0 Ub0 Uc0 Ud0 U0  real positive
    syms Vt real positive
    
    %% General functions
    % OCV functions

    % Assume identical OCV functions for now
    Ua = U0 + lbd * za;
    Ub = U0 + lbd * zb;
    Uc = U0 + lbd * zc;
    Ud = U0 + lbd * zd;
       
    % Kirchoff's voltage law within each branch
    eqnVta = Vt == Ua - Ia * Ra;
    eqnVtb = Vt == Ub - Ib * Rb;
    eqnVtc = Vt == Uc - Ic * Rc;
    eqnVtd = Vt == Ud - Id * Rd;

    
    %% Solve the system of equations for 3 cells
    
    % Eliminate Vt using Kirchoff's law
    eqnIa_temp = isolate(0 == eliminate([eqnVta, eqnVtb], Vt), Ia);
    eqnIb_temp = isolate(0 == eliminate([eqnVtb, eqnVtc], Vt), Ib);
    eqnIc_temp = isolate(0 == eliminate([eqnVtc, eqnVta], Vt), Ic);
    eqnI = I == Ia + Ib + Ic;

    eqnIa = subs(eqnIa_temp, Ib, rhs(isolate(eqnI, Ib)))
    eqnIa = subs(eqnIa, Ic, rhs(eqnIc_temp))
    eqnIa = isolate(eqnIa, Ia)

    eqnIb = subs(eqnIb_temp, Ic, rhs(isolate(eqnI, Ic)))
    eqnIb = subs(eqnIb, Ia, rhs(eqnIa_temp))
    eqnIb = isolate(eqnIb, Ib)
  
    eqnIc = subs(eqnIc_temp, Ia, rhs(isolate(eqnI, Ia)))
    eqnIc = subs(eqnIc, Ib, rhs(eqnIb_temp))
    eqnIc = isolate(eqnIc, Ic)

    % Sanity check that the total current still adds up
    simplify(rhs(eqnIa) + rhs(eqnIb) + rhs(eqnIc))
    
    % Solve for zdot:
    eqnx1 = x1 == - rhs(eqnIa)/Qa + rhs(eqnIb)/Qb;
    eqnx1 = isolate(collect(simplify(eqnx1), [za-zb zb-zc zc-za I]), x1)

    % Same cells
    isolate(simplify(subs(eqnx1, {Qa Qb Qc Ra Rb Rc}, {Q Q Q R R R})), x1)

    eqnx2 = x2 == - rhs(eqnIb)/Qb + rhs(eqnIc)/Qc;
    eqnx2 = isolate(collect(simplify(eqnx2), [za-zb zb-zc zc-za I]), x2)

    % Same cells
    isolate(simplify(subs(eqnx2, {Qa Qb Qc Ra Rb Rc}, {Q Q Q R R R})), x2)

    eqnx3 = x3 == - rhs(eqnIc)/Qc + rhs(eqnIa)/Qa;
    eqnx3 = isolate(collect(simplify(eqnx3), [za-zb zb-zc zc-za I]), x3)

    % Same cells
    isolate(simplify(subs(eqnx3, {Qa Qb Qc Ra Rb Rc}, {Q Q Q R R R})), x3)


    %% Solve the system of equations for 4 cells
    eqnIa_temp = isolate(0 == eliminate([eqnVta, eqnVtb], Vt), Ia);
    eqnIb_temp = isolate(0 == eliminate([eqnVtb, eqnVtc], Vt), Ib);
    eqnIc_temp = isolate(0 == eliminate([eqnVtc, eqnVtd], Vt), Ic);
    eqnId_temp = isolate(0 == eliminate([eqnVtd, eqnVta], Vt), Id);

    eqnI = I == Ia + Ib + Ic + Id;

    eqnIa = subs(eqnIa_temp, Ib, rhs(isolate(eqnI, Ib)))
    eqnIa = subs(eqnIa, Ic, rhs(eqnIc_temp))
    eqnIa = subs(eqnIa, Id, rhs(eqnId_temp))
    eqnIa = isolate(eqnIa, Ia)

    eqnIb = subs(eqnIb_temp, Ic, rhs(isolate(eqnI, Ic)))
    eqnIb = subs(eqnIb, Id, rhs(eqnId_temp))
    eqnIb = subs(eqnIb, Ia, rhs(eqnIa_temp))
    eqnIb = isolate(eqnIb, Ib)

    eqnIc = subs(eqnIc_temp, Id, rhs(isolate(eqnI, Id)))
    eqnIc = subs(eqnIc, Ia, rhs(eqnIa_temp))
    eqnIc = subs(eqnIc, Ib, rhs(eqnIb_temp))
    eqnIc = isolate(eqnIc, Ic)

    eqnId = subs(eqnId_temp, Ia, rhs(isolate(eqnI, Ia)))
    eqnId = subs(eqnId, Ib, rhs(eqnIb_temp))
    eqnId = subs(eqnId, Ic, rhs(eqnIc_temp))
    eqnId = isolate(eqnId, Id)

    simplify(rhs(eqnIa) + rhs(eqnIb) + rhs(eqnIc) + rhs(eqnId))

    % Solve for zdot:
    eqnx1 = x1 == - rhs(eqnIa)/Qa + rhs(eqnIb)/Qb;
    eqnx1 = isolate(collect(simplify(eqnx1), [za-zb zb-zc zc-zd zd-za I]), x1)

    isolate(simplify(subs(eqnx1, {Qa Qb Qc Qd Ra Rb Rc Rd}, {Q Q Q Q R R R R})), x1)

    eqnx2 = x2 == - rhs(eqnIb)/Qb + rhs(eqnIc)/Qc;
    eqnx2 = isolate(collect(simplify(eqnx2), [za-zb zb-zc zc-zd zd-za I]), x2)
    isolate(simplify(subs(eqnx2, {Qa Qb Qc Qd Ra Rb Rc Rd}, {Q Q Q Q R R R R})), x2)

    eqnx3 = x3 == - rhs(eqnIc)/Qc + rhs(eqnId)/Qd;
    eqnx3 = isolate(collect(simplify(eqnx3), [za-zb zb-zc zc-zd zd-za I]), x3)
    isolate(simplify(subs(eqnx3, {Qa Qb Qc Qd Ra Rb Rc Rd}, {Q Q Q Q R R R R})), x3)

    eqnx4 = x4 == - rhs(eqnId)/Qd + rhs(eqnIa)/Qa;
    eqnx4 = isolate(collect(simplify(eqnx4), [za-zb zb-zc zc-zd zd-za I]), x4)
    isolate(simplify(subs(eqnx4, {Qa Qb Qc Qd Ra Rb Rc Rd}, {Q Q Q Q R R R R})), x4)

    keyboard

end