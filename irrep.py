from dataclasses import dataclass

@dataclass 
class IrrepNames:
    wp: str          ## Irrep with parity
    np: str          ## Irrep no parity
    ferm: bool          ## Is this a double cover?
    lg: bool            ## LG?
    dim: int            ## dimension
    G: int              ## G-parity

irrep_names_no_par={
  "A1": IrrepNames(wp= "A1", np= "A1", ferm= False, lg= False, dim= 1, G= 0),
  "A2": IrrepNames(wp= "A2", np= "A2", ferm= False, lg= False, dim= 1, G= 0),
  "T1": IrrepNames(wp= "T1", np= "T1", ferm= False, lg= False, dim= 3, G= 0),
  "T2": IrrepNames(wp= "T2", np= "T2", ferm= False, lg= False, dim= 3, G= 0),
  "E":  IrrepNames(wp= "E", np= "E", ferm= False, lg= False, dim= 2, G= 0)
}

