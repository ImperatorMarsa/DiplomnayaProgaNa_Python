from multiprocessing.dummy import Pool
from math import *
import time

p=Pool()
def f(a):
	return tan(sin(a)**4-cos(a)**8)**.5, 666

mas=[9,8,7]

k=time.time()
for x in mas:
	s=f(x)
print(time.time()-k)

k=time.time()
rezult=p.map(f,mas)
print(time.time()-k)
#p.close()
#p.join()

print(rezult)