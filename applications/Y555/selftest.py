import os
import sys

if not r"." in sys.path:
    sys.path.append(r".")
import inc_p3
import p3_to_mm
import bimm

def test_all():
    inc_p3.test_all()
    p3_to_mm.test_all()
    bimm.test_all()

if __name__ == "__main__":
    test_all()


  