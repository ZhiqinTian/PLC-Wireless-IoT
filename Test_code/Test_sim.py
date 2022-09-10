from Sim.Simpy_Core import *

env=env()
def a(env):
    try:
        yield env.timeout(10)
        print('ok')
    except simpy.Interrupt:
        pass
def b(env,p):
    yield env.timeout(5)
    print('interrupt')
    p.interrupt()
    
p1=env.process(a(env))
p2=env.process(b(env,p1))
env.run(until=15)