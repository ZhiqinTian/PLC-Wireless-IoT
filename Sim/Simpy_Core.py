from collections import deque as Queue
from simpy.util import start_delayed as start_at
from simpy import AnyOf
import random
import simpy
import inspect

def set_random_seed(seed=0):
    random.seed(seed)

def env():
    return simpy.Environment()

def ensure_generator(env,func,*args,**kwargs):
    #Make sure that func is a generator function. If it is not, return a generator wrapper
    if inspect.isgeneratorfunction(func):
        return func(*args,**kwargs)
    else:
        def _wrapper():
            func(*args,**kwargs)
            yield env.timeout(0)
        return _wrapper()
