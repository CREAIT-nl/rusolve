\ File written by HiGHS .lp file handler
max
 obj: +3 x +2 y -1 z +4 b +0.5 w +1 v  +7
st
 c1: +1 x +1 y +1 z <= +10
 c2lo: +1 x -1 y >= +2
 c2up: +1 x -1 y <= +8
 c3: +1 x +1 b = +3
 c4: +1 y +1 v >= -5
 c5: -1 y +1 z >= -4
bounds
 x <= 6
 -2 <= y <= 5
 z free
 b <= 1
 w = 2
 -inf <= v <= 3
bin
 b
gen
 x
semi
end
