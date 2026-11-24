Minimize
 obj: 9 qA + 7 qB + 8 qC + 120 useB + 0.2 largest
Subject To
 demand: qA + qB + qC >= 100
 b_off: useB = 0 -> qB <= 0
 b_min: useB = 1 -> qB >= 25
Bounds
 0 <= qA <= 60
 0 <= qB <= 70
 0 <= qC <= 50
 0 <= largest <= 70
Binaries
 useB
SOS
 a_or_c: S1 :: qA : 1 qC : 2
General Constraints
 g_largest: largest = MAX ( qA , qB , qC )
End
