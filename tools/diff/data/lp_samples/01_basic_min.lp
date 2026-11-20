\ LP format for model basic_min
\ A minimal minimization with <= rows and default bounds.
Minimize
 cost: 2 x1 + 3 x2
Subject To
 cap: x1 + x2 <= 10
 mix: x1 - x2 <= 2
End
