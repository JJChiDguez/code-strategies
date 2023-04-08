/*
Strategy for computing/evaluating Richelot (2^n,2^n)-isogeny starting at a predefined abelian surface.
This abelian surface is the Jacobian of a hyperelliptic curve defined by a Type-2 equation.
 */

load "22_richelot_isogeny.m";

getCurve := function(input, x)
    A, B, C, E := Explode(input);
    f := (x^2 - 1) * (x^2 - A) * (E*x^2 - B*x + C);
    return HyperellipticCurve(f);
end function;

RichelotChain22_strategy := function(type2_invariants, kernel, n: P_list:=[], S)
    /*
     * INPUT:
     *     - type2_invariants = [A,B,C,E] defining hyperelliptic curve
     *       y^2 = (x^2-1)(x^2-A)(E*x^2-B*x+C)
     *     - kernel = [P1,P2]: Mumford coordinates of P1,P2 in defining a
     *       (2^n,2^n)-group of the Jacobian of the hyperelliptic curve
     *     - n
     *     - P_list: list of points
     *     - S: strategy as in SIDH
     * OUTPUT:
     *     - type2_invariants defining the codomain curve of the (2^n,2^n)-isogeny
     *     - Q_list: list of the image points in P_list under the isogeny
     */

    multiply := 0;  // Counter
    evaluate := 0;  // Counter

    K4 :=[];
    k := 1;
    H := type2_invariants;
    Append(~K4, kernel);
    R := kernel;
    moves := [0];
    assert #S eq (n - 1);

    for i := 0 to (#S - 1) do
        prev := &+(moves);
        while prev ne (n - 1 - i) do
            Append(~moves, S[k]);
            R4 := R;
            for j := prev to (prev + S[k] - 2) do
                R4 := DoubleType2(R4, type2_invariants);
                multiply +:= #R4;
            end for;
            Append(~K4, R4);
            R := DoubleType2(R4, type2_invariants);
            multiply +:= #R4;
            prev +:= S[k];
            k +:= 1;
        end while;

        trafo, type1_invariants := FindTransformation(type2_invariants, R[1], R[2] : P:=K4[#K4][1]);
        Prune(~moves);
        Prune(~K4);
        K4_lists := MumfordTransformationga1(trafo, K4);
        _, K4 := RichelotType1(type1_invariants, K4_lists);
        evaluate +:= #K4_lists * 2; // <--- Increase counter
        Q_lists := MumfordTransformationga1(trafo, [P_list]);
        type2_invariants, Q_lists := RichelotType1(type1_invariants, Q_lists);
        evaluate +:= #Q_lists[1]; // <--- Increase counter
        P_list := Explode(Q_lists);

        // Next branch catches special cases determined by the public strategy
        if #K4 gt 1 then
            R := DoubleType2(K4[#K4], type2_invariants);
            multiply +:= #K4[#K4];
        else
            R := K4[#K4];
        end if;
    end for;

    // last (2,2)-isogeny
    assert #K4 eq 1;
    trafo, type1_invariants := FindTransformation(type2_invariants, R[1], R[2]);
    Q_lists := MumfordTransformationga1(trafo, [P_list]);
    type2_invariants, Q_lists := RichelotType1(type1_invariants, Q_lists);
    evaluate +:= #Q_lists[1]; // <--- Increase counter

    P_list := Explode(Q_lists);

    printf "\nNumber of point multiplications [3]:\t %o\n", multiply;
    printf "Number of (3,3)-isogeny evaluations:\t %o\n", evaluate;
    return type2_invariants, P_list;
end function;