Minimize
 obj: 2 x0 + 1 i0 + 90 y0 + 2 x1 + 1 i1 + 90 y1 + 2 x2 + 1 i2 + 90 y2 + 2 x3 + 1 i3 + 90 y3 + 2 x4 + 1 i4 + 90 y4 + 2 x5 + 1 i5 + 90 y5 + 2 x6 + 1 i6 + 90 y6 + 2 x7 + 1 i7 + 90 y7
Subject To
 bal_0: x0 - i0 = 20
 link_0: x0 - 5000 y0 <= 0
 bal_1: i0 + x1 - i1 = 40
 link_1: x1 - 5000 y1 <= 0
 bal_2: i1 + x2 - i2 = 0
 link_2: x2 - 5000 y2 <= 0
 bal_3: i2 + x3 - i3 = 35
 link_3: x3 - 5000 y3 <= 0
 bal_4: i3 + x4 - i4 = 50
 link_4: x4 - 5000 y4 <= 0
 bal_5: i4 + x5 - i5 = 15
 link_5: x5 - 5000 y5 <= 0
 bal_6: i5 + x6 - i6 = 45
 link_6: x6 - 5000 y6 <= 0
 bal_7: i6 + x7 - i7 = 30
 link_7: x7 - 5000 y7 <= 0
Bounds
 0 <= x0 <= 80
 0 <= x1 <= 80
 0 <= x2 <= 80
 0 <= x3 <= 80
 0 <= x4 <= 80
 0 <= x5 <= 80
 0 <= x6 <= 80
 0 <= x7 <= 80
 0 <= i0 <= 200
 0 <= i1 <= 200
 0 <= i2 <= 200
 0 <= i3 <= 200
 0 <= i4 <= 200
 0 <= i5 <= 200
 0 <= i6 <= 200
 0 <= i7 <= 200
Binaries
 y0 y1 y2 y3 y4 y5 y6 y7
End
