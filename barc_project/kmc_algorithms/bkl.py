#!/usr/bin/env python3

import numpy as np

def bkl(procs, rnd, pred, sys):
    epsilon = 0.0000001

    sys_time = 0

    while pred(sys, sys_time):
        rates = procs.cumulativeRates(sys)
        if not rates or rates[-1] < epsilon:
            return (sys_time, False)

        uQ = rnd() * rates[-1]
        for i, r in enumerate(rates):
            if r > uQ:
                procindex = i
                break

        procs.executeSelected(procindex, sys, uQ)

        dtime = np.log(1 / rnd()) / rates[-1]
        sys.update(sys_time, dtime)

    return (sys_time, True)
