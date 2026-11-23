NAME DIET
ROWS
 N obj
 G cal
 G protein
 G calcium
 L beans_share
COLUMNS
    bread obj 0.5
    bread cal 250
    bread protein 8
    bread calcium 20
    bread beans_share -100
    milk obj 0.8
    milk cal 120
    milk protein 8
    milk calcium 300
    milk beans_share -48
    eggs obj 1.2
    eggs cal 80
    eggs protein 6
    eggs calcium 50
    eggs beans_share -32
    beans obj 0.6
    beans cal 300
    beans protein 15
    beans calcium 100
    beans beans_share 180
RHS
    rhs cal 2000
    rhs protein 60
    rhs calcium 800
BOUNDS
  UP BND bread 8
  UP BND milk 8
  UP BND eggs 8
  UP BND beans 8
ENDATA
