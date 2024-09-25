from dataclasses import dataclass
from typing import List 

class Exotraction:
    channel = str 
    irrep = str 
    nvecs = int
    nt = int 
    t_sources = List[int]
    