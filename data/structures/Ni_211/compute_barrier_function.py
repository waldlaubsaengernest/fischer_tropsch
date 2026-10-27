def compute_barrier_function(eads):
    barriers = {
               "CO->C"       : eads["C-OH"]     - eads["CO"]     -    eads["H"],
               "C->CH4"      : eads["CH2-H"]    - eads["C"]      -  3*eads["H"],
               "C+C->CCH3"   : eads["C-CH"]     - 2*eads["C"]    -    eads["H"],
               }
    return barriers
