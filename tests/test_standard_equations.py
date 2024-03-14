from prodes.calculations import standard_equations as se
from prodes.core.point import Point


def test_distance():
    a, b, c, d = Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)

    assert se.distance(a, b) == 1
    assert se.distance(a, c) == 1
    assert se.distance(a, d) == 1


def test_direction():
    """tests if the correct direction is found"""

    a, b, c, d = Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)

    assert se.direction(a, b) == (1, 0, 0)
    assert se.direction(a, c) == (0, 1, 0)
    assert se.direction(a, d) == (0, 0, 1)


def test_find_middel():
    """tests if the correct direction is found"""

    a, b, c, d = Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)

    assert se.find_middel(a, b) == (0.5, 0, 0)
    assert se.find_middel(a, c) == (0, 0.5, 0)
    assert se.find_middel(a, d) == (0, 0, 0.5)


def test_trimean():
    """tests if the correct trimean is returned"""

    test_list = [1, 2, 4, 4, 5, 6, 7, 8, 9, 10, 11]

    assert se.trimean(test_list) == 6.25