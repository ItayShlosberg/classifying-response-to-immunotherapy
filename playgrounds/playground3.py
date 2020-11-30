


class a:
    def __init__(self):
        self.v1 = 4
    def func(self):
        print("func")

class b(a):
    def __init__(self):
        a.__init__(self)
        self.v2 = 5


x = b()
x.func()
print(x.v1)

