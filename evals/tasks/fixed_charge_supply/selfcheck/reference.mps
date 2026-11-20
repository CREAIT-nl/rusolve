NAME FAC
ROWS
 N obj
 E demand_1
 E demand_2
 E demand_3
 E demand_4
 L cap_1
 L cap_2
 L cap_3
COLUMNS
    MARKER  'MARKER'  'INTORG'
    open_1 obj 100
    open_1 cap_1 -60
    open_2 obj 150
    open_2 cap_2 -90
    open_3 obj 120
    open_3 cap_3 -80
    MARKER  'MARKER'  'INTEND'
    ship_1_1 obj 4
    ship_1_1 demand_1 1
    ship_1_1 cap_1 1
    ship_1_2 obj 6
    ship_1_2 demand_2 1
    ship_1_2 cap_1 1
    ship_1_3 obj 9
    ship_1_3 demand_3 1
    ship_1_3 cap_1 1
    ship_1_4 obj 5
    ship_1_4 demand_4 1
    ship_1_4 cap_1 1
    ship_2_1 obj 5
    ship_2_1 demand_1 1
    ship_2_1 cap_2 1
    ship_2_2 obj 3
    ship_2_2 demand_2 1
    ship_2_2 cap_2 1
    ship_2_3 obj 4
    ship_2_3 demand_3 1
    ship_2_3 cap_2 1
    ship_2_4 obj 7
    ship_2_4 demand_4 1
    ship_2_4 cap_2 1
    ship_3_1 obj 8
    ship_3_1 demand_1 1
    ship_3_1 cap_3 1
    ship_3_2 obj 5
    ship_3_2 demand_2 1
    ship_3_2 cap_3 1
    ship_3_3 obj 3
    ship_3_3 demand_3 1
    ship_3_3 cap_3 1
    ship_3_4 obj 4
    ship_3_4 demand_4 1
    ship_3_4 cap_3 1
RHS
    rhs demand_1 30
    rhs demand_2 25
    rhs demand_3 40
    rhs demand_4 35
BOUNDS
  BV BND open_1
  BV BND open_2
  BV BND open_3
ENDATA
