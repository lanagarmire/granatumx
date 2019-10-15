class Vertex:

    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy
        self.next_ = None
        self.prev = None


def print_dllist(start: Vertex):
    v = start
    lst = []
    while v is not None:
        lst.append((v.dx, v.dy))
        v = v.next_
    print(lst)


def scan(start: Vertex):
    v = start
    while True:
        vv = v.next_
        if vv is None:
            break

        vvv = vv.next_
        if vvv is None:
            break

        if slope(v) < slope(vv):
            print(f"{(v.dx, v.dy)}'s slope is less than {(vv.dx, vv.dy)}, trying the latter")
            v = vv
            continue

        print(f"{(v.dx, v.dy)}'s slope is *no* less than {(vv.dx, vv.dy)}")

        # v.slope >= vv.slope, need to v with vvv
        v.next_ = vvv
        v.dx += vv.dx
        v.dy += vv.dy
        vvv.prev = v

        if v.prev is None:
            continue
        else:
            v = v.prev
            print(f"trying the previous node {(v.dx, v.dy)}")
            continue


def slope(v):
    if v.dx == 0:
        return math.inf
    else:
        return v.dy / v.dx


def insert(v, dx1, dy1, dx2, dy2):
    ddx = (dx1 + dx2) - v.dx
    ddy = (dy1 + dy2) - v.dy

    vvv = v.next_

    vv = Vertex(dx2, dy2)

    v.next_ = vv
    v.dx = dx1
    v.dy = dy1
    vv.next_ = vvv
    vv.dx = dx2
    vv.dy = dy2

    vvv.prev = vv
    vv.prev = v

    return ddx, ddy


def test2():
    # data = [(0, 0), (1, 3), (3, 4), (8, 6), (9, 7), (10, 12)]
    data = [(1, 3), (2, 1), (5, 2), (1, 1), (1, 5), (0, 1)]

    tx, ty = dx, dy = data[0]
    v = Vertex(dx, dy)
    start = v

    for dx, dy in data[1:]:
        tx += dx
        ty += dy
        vv = Vertex(dx, dy)
        v.next_ = vv
        vv.prev = v
        v = vv

    print_dllist(start)
    scan(start)
    print_dllist(start)

    print(find_lowest(start, tx, ty))


def find_lowest(start: Vertex, tx: float, ty: float):
    x, y = 0, 0
    v = start
    slope = ty / tx
    min_d = 0
    while v is not None:
        d = y - x * slope
        if d < min_d:
            min_d = d

        x += v.dx
        y += v.dy
        v = v.next_

    return min_d
