from .Simpy_Core import *
from .Time import *


# Define the event-driven simulator class
class Simulator:
    env = env()
    time_delay = env.timeout
    time_over = None
    apps = []

    # get the current simulation time
    @classmethod
    def now(cls):
        return MilliSecond(cls.env.now)

    # Creat a process starting after a time delay
    @classmethod
    def schedule(cls, delay, func, *args, **kwargs):
        if delay == 0:
            return cls.schedule_now(func, *args, **kwargs)
        func = ensure_generator(cls.env, func, *args, **kwargs)
        return start_at(cls.env, func, delay=delay)

    # Creat a process starting now
    @classmethod
    def schedule_now(cls, func, *args, **kwargs):
        func = ensure_generator(cls.env, func, *args, **kwargs)
        return cls.env.process(func)

    # Creat an event which can be triggered when a condition is satisfied
    @classmethod
    def create_event(cls):
        return cls.env.event()

    # Add node in the simulator
    @classmethod
    def add_app(cls, app):
        cls.apps.append(app)

    # Start the simulator by execute the network applications
    @classmethod
    def run(cls):
        for app in cls.apps:
            app.run()
        if cls.time_over:
            cls.env.run(until=cls.time_over.sim_time)
        else:
            cls.env.run()

    class timer:
        pass
