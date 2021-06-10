from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import numpy as np
from enum import Enum
from textwrap import dedent, wrap
from progress.bar import Bar
import pickle

class chemical_category(Enum):
    PRODUCT = 0
    REACTANT = 1

# k is default unit unless otherwise specified page (14 pdf) (8 index)

@dataclass(init=True)
class Interval:
    temperature_range: Tuple[float] = field(init=False)      # temperature range K
    t_Exponents      : List[float]  = field(init=False)             # temperature coe
    cp_coefficients  : List[float]  = field(init=False)             # Cp R coe Rankine? assuming kelvin
    b_constants      : List[float]  = field(init=False)             # integration constants
    enthalpy         : float        = field(init=False)             #H(298.15)-H(0)

@dataclass(init=True)
class Species:
    # calcualtes all in quantity / RT
    name                : str = field(init=False) 
    info                : str = field(init=False) 
    source              : str = field(init=False) 
    intervals_cnt       : int = field(init=False) 
    elements            : Dict[str,int] = field(init=False) 
    molecular_weight    : float = field(init=False)     # g/mol
    heat_of_formation   : float = field(init=False)   # Jmol
    phase               : float  = field(init=False)   # zero for gas, non-zero for condensed phase
    species_type        : chemical_category = field(init=False) 
    intervals           : List[Interval] = field(init=False) 

    def _get_interval(self, temperature):
        _in = lambda temp, rng: True if (temp > min(rng) and temp <= max(rng)) else False
        for interval in self.intervals:
            if _in(temperature, interval.temperature_range):
                return interval
        else:
            raise ValueError("Missing Thermochemical Data")

    # https://shepherd.caltech.edu/EDL/PublicResources/sdt/refs/NASA-RP-1311-2.pdf
    # equations from page (80 pdf) (74 index)

    def HeatCapacity(self, temperature, interval=None):
        # Cp (pressure)
        interval = self._get_interval(temperature) if interval == None else interval
        return np.sum(
            [
                a*temperature**exp
                for a, exp in zip(
                    interval.cp_coefficients,
                    interval.t_Exponents
                )
            ]
        )

    def Enthalpy(self, temperature, interval=None):
        interval = self._get_interval(temperature) if interval == None else interval
        a = interval.cp_coefficients
        T = temperature
        b = interval.b_constants
        return -a[0]*T**-2 + a[1]*np.log(T)*T**-1 \
            +a[2] + a[3]*(T/2) + a[4]*(T**2)/3 + a[5]*(T**3)/4 \
                +a[6]*(T**4)/5 + b[0]/T

    def Entropy(self, temperature, interval=None):
        # fixed at P=1atm
        interval = self._get_interval(temperature) if interval == None else interval
        a = interval.cp_coefficients
        T = temperature
        b = interval.b_constants
        return -a[0]*(T**-2)/2 - a[1]*(T**-1) + a[2]*np.log(T) \
            +a[3]*T + a[4]*(T**2)/2 + a[5]*(T**3)/3 + b[1]

    def gibbs_star(self, temperature):
        # eq 16
        interval = self._get_interval(temperature)
        return self.heat_of_formation + (self.Enthalpy(temperature, interval=interval)- interval.enthalpy) + temperature*self.Entropy(temperature, interval=interval)

    def gibbs(self, temperature, N_ratio, pressure_ratio, R):
        # pressure ratio = (species pressure) / (mixture pressure)
        # mol ratio      = (species mols)     / (mixture mols)
        return self.gibbs_star(temperature) + np.log(N_ratio) + np.log(pressure_ratio) 

    def element_cnt(self, elm:str):
        try:
            return self.elements[elm]
        except KeyError as e:
            return 0

from functools import lru_cache

@dataclass
class SpeciesMixture:
    # for increased ratio of fuel to oxy, add species multiple times
    species: List[Species]
    #caching for faster programming
    Nj: Dict[str, float]
    #end composition mixture for products is default all 1 mol at beginning
    Sdict = Dict[str,Species]

    def _build_species_dict(self):
        for s in self.species:
            self.SDict[s.name] = s

    def initial_product_moles(self, guess = None) -> None:
        if guess == None:
            #crude guess
            self.Nj = {s.name:0.1/len(self.species) for s in self.species}
        else:
            self.Nj = guess
        self.N = 0.1

    @lru_cache
    def _get_species(self, species:str):
        for s in self.species:
            if species == s.name:
                return s
            return None
    
    @lru_cache()
    def molecular_weight(self) -> float:
        return np.sum(
            [
                species.molecular_weight * self.Nj*[species.name]/self.N
                for species in self.species
            ]
        )

    def Nj(self) -> dict:
        return self.Nj

    @lru_cache()
    def nfraction(self, species:str) -> float:
        return self.Nj[species]/self.N

    @lru_cache()
    def get_elements_species(self, species:str) -> dict:
        return self._get_elements(species).elements

    @lru_cache()
    def get_elements_global(self):
        elements_global = {}
        for s in self.species:
            for elm, cnt in s.elements.items():
                x = self.mols[s.name] * cnt
                if elm in elements_global.keys():
                    elements_global[elm] += x
                else:
                    elements_global[elm] = x
        else:
            return elements_global

@dataclass(init=True)
class TransportData:
    temperature_range   : Tuple[float] = field(init=False)
    coefficients        : List[float] = field(init=False)

@dataclass(init=True)
class TransportSpecies:
    name        : str = field(init=False)
    secondary   : str = field(init=False)
    info        : str = field(init=False)
    V           : List[TransportData] = field(init=False)
    C           : List[TransportData] = field(init=False)

    #https://shepherd.caltech.edu/EDL/PublicResources/sdt/refs/NASA-RP-1311-2.pdf
    #page (102 pdf) (96 index)

    def _eval(self, T, consts):
        return consts[0]*np.log(T) + consts[1]/T + consts[2]/(T**2) + consts[3]

    def _get_interval(self, ref, temperature):
        _in = lambda temp, rng: True if (temp > min(rng) and temp <= max(rng)) else False
        for interval in ref:
            if _in(temperature, interval.temperature_range):
                return interval
            raise ValueError("Missing Transport Data")

    def Viscosity(self, temperature):
        # ln(Nu)
        interval = self._get_interval(self.V, temperature)
        return self._eval(temperature, interval.coefficients)

    def ThermalConductivity(self, temperature):
        # ln(Lambda)
        interval = self._get_interval(self.C, temperature)
        return self._eval(temperature, interval.coefficients)


@dataclass(init=True)
class ThermoDatabase:
    ThermoINPinfo   : str = field(init=False) 
    TransINPinfo    : str = field(init=False) 
    Products        : Dict[str,Species] = field(init=False) 
    Reactants       : Dict[str,Species] = field(init=False)
    Transport       : Dict[str, List[TransportSpecies]] = field(init=False) 