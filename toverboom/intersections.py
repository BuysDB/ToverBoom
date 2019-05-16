#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

# Obtained from: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
def on_segment(p, q, r):
    '''Given three colinear points p, q, r, the function checks if
    point q lies on line segment "pr"
    '''
    if (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
        q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1])):
        return True
    return False

def orientation(p, q, r):
    '''Find orientation of ordered triplet (p, q, r).
    The function returns following values
    0 --> p, q and r are colinear
    1 --> Clockwise
    2 --> Counterclockwise
    '''

    val = ((q[1] - p[1]) * (r[0] - q[0]) -
            (q[0] - p[0]) * (r[1] - q[1]))
    if val == 0:
        return 0  # colinear
    elif val > 0:
        return 1   # clockwise
    else:
        return 2  # counter-clockwise


def do_intersect(p1, q1, p2, q2):
    '''Main function to check whether the closed line segments p1 - q1 and p2
       - q2 intersect'''
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if (o1 != o2 and o3 != o4):
        return True

    # Special Cases
    # p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 and on_segment(p1, p2, q1)):
        return True

    # p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 and on_segment(p1, q2, q1)):
        return True

    # p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 and on_segment(p2, p1, q2)):
        return True

    # p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 and on_segment(p2, q1, q2)):
        return True

    return False # Doesn't fall in any of the above cases

def expandStart( x0, y0, x1, y1):
    return ( (x0+(x1-x0)*0.01 , y0+(y1-y0)*0.01), (x1, y1) )

def test_intersect_func():
    p1 = (1, 1)
    q1 = (1, 3)
    p2 = (1, 1)
    q2 = (10, 2)

    p1,q1 = expandStart( *(p1+q1) )
    p2,q2 = expandStart( *(q2+q2) )

    fig, ax = plt.subplots()
    ax.plot([p1[0], q1[0]], [p1[1], q1[1]], 'x-')
    ax.plot([p2[0], q2[0]], [p2[1], q2[1]], 'x-')
    print(do_intersect(p1, q1, p2, q2))

    p1 = (-1, 1)
    q1 = (1, 3)
    p2 = (-1, 1)
    q2 = (10, 2)



    fig, ax = plt.subplots()
    ax.plot([p1[0], q1[0]], [p1[1], q1[1]], 'x-')
    ax.plot([p2[0], q2[0]], [p2[1], q2[1]], 'x-')
    print(do_intersect(p1, q1, p2, q2))
    plt.show()
