Minimize
 obj: 8 ship_0_0 + 6 ship_0_1 + 10 ship_0_2 + 9 ship_0_3 + 9 ship_1_0 + 12 ship_1_1 + 13 ship_1_2 + 7 ship_1_3 + 14 ship_2_0 + 9 ship_2_1 + 16 ship_2_2 + 5 ship_2_3
Subject To
 supply_0: ship_0_0 + ship_0_1 + ship_0_2 + ship_0_3 = 30
 supply_1: ship_1_0 + ship_1_1 + ship_1_2 + ship_1_3 = 40
 supply_2: ship_2_0 + ship_2_1 + ship_2_2 + ship_2_3 = 30
 demand_0: ship_0_0 + ship_1_0 + ship_2_0 = 20
 demand_1: ship_0_1 + ship_1_1 + ship_2_1 = 25
 demand_2: ship_0_2 + ship_1_2 + ship_2_2 = 30
 demand_3: ship_0_3 + ship_1_3 + ship_2_3 = 25
End
