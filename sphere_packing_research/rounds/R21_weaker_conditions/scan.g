# R21: scan the GAP character table library for "admissible" pairs (G, chi).
# W = rational representation whose character is (Schur index) x tau, tau = sum of the Galois conjugates of chi.
# Necessary conditions for tensor stability (BS) used here (see notes.md):
#  (C1) every class of prime order q:  tau(g) = 0 (regular)  or  tau(g) = -tau(1)/(q-1) (free);
#  (C1') every class of composite order m: fixed dimension a = (1/m) sum_{d|m} phi(d) tau(g^(m/d)) is 0 or >= tau(1)/m.
# gain index  g = |G| / (2 dim_Q A) = |G| / (2 c chi(1)^2),  c = number of Galois conjugates.
LoadPackage("ctbllib", false);
LIMIT := 1200;
out := OutputTextFile("scan_out.txt", false);
SetPrintFormattingStatus(out, false);
names := AllCharacterTableNames(IsDuplicateTable, false);
PrintTo(out, "# tables: ", Length(names), "\n");
ntab := 0; npairs := 0;
for name in names do
  t := CharacterTable(name);
  if t = fail then continue; fi;
  if NrConjugacyClasses(t) > LIMIT then AppendTo(out, "# SKIP (too many classes) ", name, "\n"); continue; fi;
  ntab := ntab + 1;
  ords := OrdersClassRepresentatives(t);
  size := Size(t);
  irr := Irr(t);
  ncl := Length(ords);
  done := [];
  ok_pm := true;
  for i in [1..Length(irr)] do
    if i in done then continue; fi;
    chi := irr[i];
    fam := Set(List(PrimeResidues(Conductor(ValuesOfClassFunction(chi))), k -> List(ValuesOfClassFunction(chi), x -> GaloisCyc(x, k))));
    c := Length(fam);
    for x in fam do
      pos := Position(List(irr, ValuesOfClassFunction), x);
      if pos <> fail then AddSet(done, pos); fi;
    od;
    if Length(ClassPositionsOfKernel(chi)) > 1 then continue; fi;    # faithful only
    tau := Sum(fam);
    n := tau[1];
    gain := size / (2 * c * chi[1]^2);
    if gain < 1 then continue; fi;
    adm := true; nfree := 0; nreg := 0; freeprimes := []; regprimes := [];
    for j in [2..ncl] do
      q := ords[j];
      if IsPrimeInt(q) then
        if tau[j] = 0 then nreg := nreg + 1; AddSet(regprimes, q);
        elif tau[j] * (q - 1) = -n then nfree := nfree + 1; AddSet(freeprimes, q);
        else adm := false; break; fi;
      fi;
    od;
    if not adm then continue; fi;
    # composite orders
    cyc := true; nz := 0;
    for j in [2..ncl] do
      m := ords[j];
      if not IsPrimeInt(m) then
        s := 0;
        for d in DivisorsInt(m) do
          pm := PowerMap(t, m / d);
          if pm = fail then ok_pm := false; s := fail; break; fi;
          s := s + Phi(d) * tau[pm[j]];
        od;
        if s = fail then cyc := fail; break; fi;
        a := s / m;
        if a <> 0 and a * m < n then cyc := false; fi;
      fi;
    od;
    npairs := npairs + 1;
    AppendTo(out, name, " | size=", size, " | deg=", chi[1], " | c=", c, " | ind=", Indicator(t, [chi], 2)[1],
       " | gain=", gain, " ~ ", Float(gain), " | free=", freeprimes, " | reg=", regprimes, " | cyclic_ok=", cyc, "\n");
  od;
od;
AppendTo(out, "# scanned tables: ", ntab, ", admissible pairs with gain>=1: ", npairs, "\n");
CloseStream(out);
QUIT;
