NAME BAKERY
OBJSENSE
    MAX
ROWS
 N obj
 L oven
 L flour
 L rolls_cap
COLUMNS
    MARKER  'MARKER'  'INTORG'
    loaves obj 30
    loaves oven 3
    loaves flour 5
    loaves rolls_cap -1.5
    rolls obj 18
    rolls oven 1
    rolls flour 2
    rolls rolls_cap 1
    cakes obj 45
    cakes oven 4
    cakes flour 6
    MARKER  'MARKER'  'INTEND'
RHS
    rhs oven 40
    rhs flour 70
BOUNDS
  UP BND loaves 1000
  UP BND rolls 1000
  UP BND cakes 4
ENDATA
