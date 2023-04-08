strategy := function(n: mu := 1, eta := 1)
	// This is a dynamic programming algorithm for computing optimal strategies
	// inputs: mu (cost of scalar multiplication) and eta (cost of isogeny evaluation)
	// output: the optimal strategy and its cost

	// strategy per exponent
	S := AssociativeArray();
	S[1] := [];
	// strategy cost per exponent
	C := AssociativeArray();
	C[1] := 0;
	for i:=2 to (n+1) do
		cost, b := Minimum([ C[i-b] + C[b] + b*mu + (i-b)*eta : b in [1 .. (i-1)]]);
		S[i] := [b] cat S[i-b] cat S[b];
        C[i] := cost;
	end for;
	return S[n], C[n];
end function;

balanced_strategy := function(n)
	// This is a dynamic programming algorithm for computing optimal strategies
	// strategy per exponent
	S := AssociativeArray();
	S[1] := [];
	for i:=2 to (n+1) do
		b := Integers()!(Ceiling(i / 2));
		S[i] := [b] cat S[i-b] cat S[b];
	end for;
	return S[n];
end function;