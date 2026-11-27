# R21, rigorous scan (after the critic pass).  Only PROVED necessary conditions for a unimodular G-lattice with (BS):
#  (R1) every class of prime order q:  tau = 0 (regular) or tau = -n/(q-1) (free)                 [Round 5 rigidity]
#  (R2) class g of order m, prime l with l || m (exactly once) and g^(m/l) regular  ==>  tau(g) = 0   [fibrewise lemma]
#  (R3) class g of prime-power order l^i, i >= 2, with g^(m/l) regular            ==>  tau(g) = 0   [prime-power lemma]
# 'unres' counts classes with nonzero trace that are neither totally free nor covered by (R2)/(R3)
# (every regular prime divides the order at least twice and the order is not a prime power): not tested.
# 'strong' is the old heuristic flag of scan3.g (Claim E for all classes), kept for comparison.
LoadPackage("ctbllib", false);
Analyze := function(t, name, out, mingain)
  local ords, size, irr, ncl, done, i, chi, vals, fam, c, x, pos, tau, n, gain, adm, j, q, freeprimes, regprimes, isfree, isreg,
        rig, strong, unres, m, l, pm, tf, cov, want, cnt, anyreg;
  ords := OrdersClassRepresentatives(t); size := Size(t); irr := Irr(t); ncl := Length(ords); done := []; cnt := 0;
  for i in [1..Length(irr)] do
    if i in done then continue; fi;
    chi := irr[i]; vals := ValuesOfClassFunction(chi);
    fam := Set(List(PrimeResidues(Conductor(vals)), k -> List(vals, x -> GaloisCyc(x, k))));
    c := Length(fam);
    for x in fam do pos := Position(List(irr, ValuesOfClassFunction), x); if pos <> fail then AddSet(done, pos); fi; od;
    if Length(ClassPositionsOfKernel(chi)) > 1 then continue; fi;
    tau := Sum(fam); n := tau[1];
    gain := size / (2 * c * chi[1]^2);
    if gain < mingain then continue; fi;
    adm := true; freeprimes := []; regprimes := [];
    isfree := List([1..ncl], j -> false); isreg := List([1..ncl], j -> false);
    for j in [2..ncl] do
      q := ords[j];
      if IsPrimeInt(q) then
        if tau[j] = 0 then AddSet(regprimes, q); isreg[j] := true;
        elif tau[j] * (q - 1) = -n then AddSet(freeprimes, q); isfree[j] := true;
        else adm := false; break; fi;
      fi;
    od;
    if not adm then continue; fi;
    rig := true; strong := true; unres := 0;
    for j in [2..ncl] do
      m := ords[j];
      if not IsPrimeInt(m) then
        tf := true; cov := false; anyreg := false;
        for l in Set(Factors(m)) do
          pm := PowerMap(t, m / l);
          if isreg[pm[j]] then
            tf := false; anyreg := true;
            if Number(Factors(m), y -> y = l) = 1 then cov := true; fi;       # (R2)
            if IsPrimePowerInt(m) then cov := true; fi;                         # (R3)
          fi;
        od;
        if tf then want := n * MoebiusMu(m) / Phi(m); else want := 0; fi;
        if tau[j] <> want then strong := false; fi;
        if tf and tau[j] <> want then Print("INTERNAL: totally free class with unexpected trace in ", name, "\n"); fi;
        if cov and tau[j] <> 0 then rig := false; fi;
        if anyreg and not cov and tau[j] <> 0 then unres := unres + 1; fi;
      fi;
    od;
    cnt := cnt + 1;
    AppendTo(out, name, " | size=", size, " | deg=", chi[1], " | c=", c, " | ind=", Indicator(t, [chi], 2)[1],
       " | gain=", gain, " ~ ", Float(gain), " | free=", freeprimes, " | reg=", regprimes,
       " | rigorous=", rig, " | strong=", strong, " | unres=", unres, "\n");
  od;
  return cnt;
end;
out := OutputTextFile("scan4_lib.txt", false); SetPrintFormattingStatus(out, false);
names := AllCharacterTableNames(IsDuplicateTable, false);
tot := 0; ntab := 0;
for name in names do
  t := CharacterTable(name);
  if t = fail or NrConjugacyClasses(t) > 1200 then continue; fi;
  ntab := ntab + 1;
  tot := tot + Analyze(t, name, out, 1);
od;
AppendTo(out, "# tables scanned: ", ntab, "; pairs passing (R1) with gain >= 1: ", tot, "\n");
CloseStream(out);
QUIT;
