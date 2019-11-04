#!/usr/bin/env python

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

for i in range(start,end):
    if i % 7 != 0:
        print(i)
