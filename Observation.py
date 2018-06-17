from typing import List

import aminoAcidEncoding as Encode

class Observation(object):
    def __init__(self,
                 peptide: str,
                 ic50: float,
                 target: int = None,
                 bit_code: List[bool] = None,
                 blo_map: List[float] = None,
                 blosum : List[int] = None):
        """
        Creates a new instance of the Peptide scoring object
        :param observations:
        :param observations_ic50:
        :param observations_bitcode:
        :param observations_blomap:
        :param target:
        """
        if bit_code is None:
            bit_code = Encode.twentybit_encode(peptide)
        if blo_map is None:
            blo_map = Encode.blomap_encode(peptide)
        if blosum is None:
            blosum = Encode.blosum_encode(peptide)
        self.peptide  : str         = peptide
        self.ic50     : float       = ic50
        self.bit_code : List[bool]  = bit_code
        self.blo_map  : List[float] = blo_map
        self.target   : int         = target
        self.blosum   : List[int]   = blosum
        pass

    def __str__(self):
        return "{} IC50: {},\nbit code: {},\nblomap: {}\nblosum: {}\ntarget: {}".format(
            self.peptide, self.ic50, self.bit_code, self.blo_map, self.blosum, self.target)