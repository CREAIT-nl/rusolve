NAME INFPLAN
ROWS
 N obj
 G dem1
 G dem2
 L cap1
 L cap2
COLUMNS
    a obj 2
    a dem1 1
    a cap1 1
    b obj 3
    b dem1 1
    b cap2 1
    c obj 1
    c dem2 1
    c cap1 1
    d obj 4
    d dem2 1
    d cap2 1
RHS
    rhs dem1 30
    rhs dem2 25
    rhs cap1 40
    rhs cap2 40
BOUNDS
 UP bnd a 15
 UP bnd b 10
 UP bnd c 20
 UP bnd d 20
ENDATA
