* Source:     Pyomo MPS Writer
* Format:     Free MPS
*
NAME unknown
OBJSENSE
 MAX
ROWS
 N  o
 L  c_u_c1_
 G  r_l_c2_
 L  r_u_c2_
 E  c_e_c3_
 G  c_l_c4_
 G  c_l_c5_
 E  c_e_ONE_VAR_CONSTANT
COLUMNS
     x o 3
     x c_u_c1_ 1
     x r_l_c2_ 1
     x r_u_c2_ 1
     x c_e_c3_ 1
     y o 2
     y c_u_c1_ 1
     y r_l_c2_ -1
     y r_u_c2_ -1
     y c_l_c4_ 1
     y c_l_c5_ -1
     z o -1
     z c_u_c1_ 1
     z c_l_c5_ 1
     b o 4
     b c_e_c3_ 1
     w o 0.5
     v o 1
     v c_l_c4_ 1
     ONE_VAR_CONSTANT o 7
     ONE_VAR_CONSTANT c_e_ONE_VAR_CONSTANT 1
RHS
     RHS c_u_c1_ 10
     RHS r_l_c2_ 2
     RHS r_u_c2_ 8
     RHS c_e_c3_ 3
     RHS c_l_c4_ -5
     RHS c_l_c5_ -4
     RHS c_e_ONE_VAR_CONSTANT 1
BOUNDS
 LI BOUND x 0
 UI BOUND x 6
 LO BOUND y -2
 UP BOUND y 5
 FR BOUND z
 BV BOUND b
 LO BOUND w 2
 UP BOUND w 2
 MI BOUND v
 UP BOUND v 3
ENDATA
