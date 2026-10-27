def compute_barrier_function(eads):
    barriers = {
               "CO->C"       :   eads["O-H+H+C+H2O(g)"] - eads["CO+2H+H2O(g)"],
               "C->CH4"      :   eads["CH3-H"]          -   eads["CH3"]    -     eads["H"],
               "C+C->CCH3"   :   eads["C-CH"]           - 2*eads["CH3"]    +   5*eads["H"]
               }
    return barriers
