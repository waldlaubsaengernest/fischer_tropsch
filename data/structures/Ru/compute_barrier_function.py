def compute_barrier_function(eads):
    barriers = {
               "CO->C"       :  eads["C-OH"]            - eads["CO"]     -     eads["H"],
               "C->CH4"      :  eads["CH3-H"]           - eads["CH"]     -   3*eads["H"],
               "C+C->CCH3"   :  eads["C-CH"]            + eads["H"]      -   2*eads["CH"]
               }
    return barriers
