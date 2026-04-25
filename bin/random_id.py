#!/usr/bin/env python
import string
import random
import sys
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

sys.stdout.write(id_generator(9))
