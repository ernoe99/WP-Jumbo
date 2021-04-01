class System:
    def __init__(self, components, actors, sensors):
        self.components = components
        self.actors = actors
        self.sensors = sensors

    def calculate(self):
        for component in self.components
            component.setactors()
            component.calculate()
            component.setsensors()

class SimpleSystem (System):
    def __init__(self, components, actors, sensors):
        super().__init__(self, components, actors, sensors)

class Component:
    def __init__(self):
        self.Interface_in()









