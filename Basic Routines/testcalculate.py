# Test von Heatpump calculate

class HPsuper:
    def __init__(self):
       self.pc = 1000.0

    def calc(self, x, e, comp):
        return [x, e, comp, self.pc]


class HP(HPsuper):
    def __init__(self):
        super().__init__()

    def calc(self):
        return super().calc(700.0, 800.0, 900.0)



hp1 = HP()
hp2 = HPsuper()

[a, b, c, d] = hp2.calc(100.0, 200.0, 300.0)

print(a)
print(b)
print(c)
print(d)


[a, b, c, d] = hp1.calc()

print(a)
print(b)
print(c)
print(d)


