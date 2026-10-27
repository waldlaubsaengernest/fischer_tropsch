def compute_barrier_function(eads):
    barriers = {
               "CO->C"       :  eads["O-H"] + eads["C"]    - eads["CO"]     -    eads["H"],
               "C->CH4"      :  eads["CH3-H_s2"]           - eads["C"]      -  4*eads["H"],
               "C+C->CCH3"   :  eads["CCH-H"]    - 2*eads["C"]    -  2*eads["H"],
               }
    return barriers
