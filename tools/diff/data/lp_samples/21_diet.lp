\ LP format for model diet
\ Mixed senses, ranged variable bounds, maximize.
Maximize
 profit: 2 x1 + 3 x2 + 1 x3
Subject To
 fat: 2 x1 + 4 x2 + 0.5 x3 <= 10
 prot: x1 + x2 + 2 x3 >= 4
 exact: x1 - x2 = 0.5
Bounds
 x2 <= 3
 0.5 <= x3 <= 2.5
 x1 free
End
