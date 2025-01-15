class position:
    x = 0.0
    v = 0.0
    a = 0.0
    mass = 0.0

    def __init__(self, x, v, a, mass):
        self.x = x
        self.v = v
        self.a = a
        self.mass = mass

    def update(self, dt, f):
        self.v += self.a * dt
        self.x += self.v * dt
        self.a = f / self.mass


class spring:
    length = 0.0
    k = 0.0
    pos = position(0.0, 0.0, 0.0, 0.0)

    def __init__(self, length, k, x, mass):
        self.length = length
        self.k = k
        self.pos = position(x, 0.0, 0.0, mass)
        self.x1 = x - self.length / 2
        self.x2 = x + self.length / 2

    def update(self, dt, f):
        self.pos.update(dt, f)

    def force(self, x1, x2) -> float:
        return self.k * (x2 - x1 - self.length)


class chain:
    spring_list = []
    length = 0.0
    k = 0.0
    mass = 0.0
    num = 0
    particle = position(0.0, 0.0, 0.0, 0.0)

    def __init__(self, length, k, mass, num, M, l0):
        self.length = length
        self.k = k
        self.mass = mass
        self.num = num

        k1 = k * num
        l1 = length / num
        m1 = mass / num

        for i in range(num):
            self.spring_list.append(spring(l1, k1, l0 / num * (i + 0.5), m1))

        self.particle = position(l0, 0.0, (l0 - length) * k / M, M)

    def update(self, dt) -> float:
        x = [0.0]
        for i in range(self.num):
            x.append(2 * self.spring_list[i].pos.x - x[-1])

        force = []
        for i in range(self.num):
            force.append(self.spring_list[i].force(x[i], x[i + 1]))

        self.spring_list[0].update(dt, -force[0] + force[1])
        for i in range(1, self.num - 1):
            self.spring_list[i].update(dt, -force[i - 1] + force[i + 1])
        self.spring_list[self.num - 1].update(
            dt, -force[self.num - 2] - self.particle.a * self.particle.mass
        )

        self.particle.update(dt, -force[self.num - 1])
        return self.particle.x


def main():
    length = 1.0
    k = 1.0
    mass = 1.0
    num = 100
    M = 2.0
    l0 = 1.2
    dt = 0.1
    chain_ = chain(length, k, mass, num, M, l0)
    for i in range(100):
        print(chain_.update(dt))


if __name__ == "__main__":
    main()
