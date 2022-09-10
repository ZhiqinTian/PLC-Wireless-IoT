import numpy as np
import matlab.engine

eng = matlab.engine.start_matlab()
dic1={'name':matlab.double([[1,2,3],[1,2,4]])}
dic2={'name':matlab.double([2,3,4])}
dic=[dic1,dic2]
d=eng.test_m1(dic)

print(d)
# d = eng.test_m1(1.0,5.0,nargout=1)
# c = eng.test_m2()
# d1=np.array(d)
# print(d1)
# print(type(d),d,c)
# eng.exit()