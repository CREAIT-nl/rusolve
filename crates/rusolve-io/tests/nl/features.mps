* Source:     Pyomo MPS Writer
* Format:     Free MPS
*
NAME unknown
OBJSENSE
 MIN
ROWS
 N  obj
 E  c_e_eq_
 G  r_l_rng_
 L  r_u_rng_
 L  c_u_le_
 G  c_l_ge_
 G  c_l_cst_
 L  c_u_dup_
 G  r_l_box_
 L  r_u_box_
 E  c_e_ONE_VAR_CONSTANT
COLUMNS
     i obj 2
     i c_e_eq_ 1
     i c_l_ge_ -1
     i c_u_dup_ 2
     b obj -3
     b c_e_eq_ 1
     b c_u_dup_ -1
     f obj 1
     f r_l_rng_ 1
     f r_u_rng_ 1
     f c_u_le_ 1
     f r_l_box_ 1
     f r_u_box_ 1
     u obj -1
     u r_l_rng_ 1
     u r_u_rng_ 1
     u c_l_ge_ 1
     u r_l_box_ 1
     u r_u_box_ 1
     w obj 0.25
     w c_l_cst_ 1
     p obj 1
     p r_l_rng_ -1
     p r_u_rng_ -1
     p c_l_cst_ 1
     p r_l_box_ 1
     p r_u_box_ 1
     q obj 4
     q c_e_eq_ 1
     q c_u_le_ 2
     ONE_VAR_CONSTANT obj -3.5
     ONE_VAR_CONSTANT c_e_ONE_VAR_CONSTANT 1
RHS
     RHS c_e_eq_ 2
     RHS r_l_rng_ -1
     RHS r_u_rng_ 6
     RHS c_u_le_ 3
     RHS c_l_ge_ -7
     RHS c_l_cst_ 2
     RHS c_u_dup_ 9
     RHS r_l_box_ -20
     RHS r_u_box_ 20
     RHS c_e_ONE_VAR_CONSTANT 1
BOUNDS
 LI BOUND i -5
 UI BOUND i 5
 BV BOUND b
 FR BOUND f
 MI BOUND u
 UP BOUND u 4
 LO BOUND w 1.5
 UP BOUND w 1.5
 LO BOUND p 0
 LO BOUND q -3
 UP BOUND q 2
ENDATA
