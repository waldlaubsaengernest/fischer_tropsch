def compute_barrier_function(eads):
    barriers = {
               "CO->C"       :  eads["OH-H+C+H2O(g)"]   - eads["CO+2H+H2O(g)"],
               "C->CH4"      :  eads["CH2-H"]           - eads["CH"]     -   2*eads["H"],
               "C+C->CCH3"   :  eads["CH-CH"]           - 2*eads["CH"]
               }
    return barriers
