NAME PROJ
OBJSENSE
    MAX
ROWS
 N obj
 L budget
 L E_needs_A
 L F_needs_B
 L F_needs_C
 L C_or_D
 G at_least_two
COLUMNS
    MARKER  'MARKER'  'INTORG'
    pA obj 40
    pA budget 30
    pA at_least_two 1
    pA E_needs_A -1
    pB obj 35
    pB budget 25
    pB at_least_two 1
    pB F_needs_B -1
    pC obj 55
    pC budget 40
    pC at_least_two 1
    pC F_needs_C -1
    pC C_or_D 1
    pD obj 50
    pD budget 35
    pD at_least_two 1
    pD C_or_D 1
    pE obj 70
    pE budget 45
    pE at_least_two 1
    pE E_needs_A 1
    pF obj 90
    pF budget 50
    pF at_least_two 1
    pF F_needs_B 1
    pF F_needs_C 1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs budget 120
    rhs C_or_D 1
    rhs at_least_two 2
BOUNDS
  BV BND pA
  BV BND pB
  BV BND pC
  BV BND pD
  BV BND pE
  BV BND pF
ENDATA
