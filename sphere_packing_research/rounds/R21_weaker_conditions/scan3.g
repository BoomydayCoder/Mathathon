# R21, second scan: stronger cyclic test (C1'') and the family SL(2,q) computed from the groups.
# (C1'') for every class g of order m:  if all prime-order powers of g are free then tau(g) = n*mu(m)/phi(m), else tau(g) = 0.
LoadPackage("ctbllib", false);
Analyze := function(t, name, out, mingain)
  local ords, size, irr, ncl, done, i, chi, vals, fam, c, x, pos, tau, n, gain, adm, j, q, freeprimes, regprimes, isfreeclass,
        strong, m, tf, l, pm, want, cnt, sylc;
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
    adm := true; freeprimes := []; regprimes := []; isfreeclass := List([1..ncl], j -> false);
    for j in [2..ncl] do
      q := ords[j];
      if IsPrimeInt(q) then
        if tau[j] = 0 then AddSet(regprimes, q);
        elif tau[j] * (q - 1) = -n then AddSet(freeprimes, q); isfreeclass[j] := true;
        else adm := false; break; fi;
      fi;
    od;
    if not adm then continue; fi;
    strong := true; sylc := true;
    for j in [2..ncl] do if isfreeclass[j] and (size / SizesCentralizers(t)[j]) mod ords[j] = 0 then sylc := false; fi; od;
    for j in [2..ncl] do
      m := ords[j];
      if not IsPrimeInt(m) then
        tf := true;
        for l in Set(Factors(m)) do
          pm := PowerMap(t, m / l);
          if not isfreeclass[pm[j]] then tf := false; fi;
        od;
        if tf then want := n * MoebiusMu(m) / Phi(m); else want := 0; fi;
        if tau[j] <> want then strong := false; fi;
      fi;
    od;
    cnt := cnt + 1;
    AppendTo(out, name, " | size=", size, " | deg=", chi[1], " | c=", c, " | ind=", Indicator(t, [chi], 2)[1],
       " | gain=", gain, " ~ ", Float(gain), " | free=", freeprimes, " | reg=", regprimes, " | strong=", strong, " | freeSylowCentral=", sylc, "\n");
  od;
  return cnt;
end;

out := OutputTextFile("scan3_lib.txt", false); SetPrintFormattingStatus(out, false);
names := AllCharacterTableNames(IsDuplicateTable, false);
tot := 0;
for name in names do
  t := CharacterTable(name);
  if t = fail or NrConjugacyClasses(t) > 1200 then continue; fi;
  tot := tot + Analyze(t, name, out, 1);
od;
AppendTo(out, "# total admissible pairs: ", tot, "\n");
CloseStream(out);
QUIT;
