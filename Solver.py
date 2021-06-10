import numpy as np
from Database import *
from dataclasses import dataclass
# http://courses.washington.edu/mengr524/sp15/handouts/gibbs_min.pdf
# following this generalized approach
from functools import lru_cache


@dataclass
class Problem:
    # organizes the problem
    fixed_pressure: float
    fixed_temperature: float
    Product_Mixture: SpeciesMixture
    Reactant_Mixtture: SpeciesMixture

    def gibbs_product(self, species):
        species = self.Product_Mixture._get_species(species)
        n_ratio = self.Product_Mixture.nfraction(species.name)
        return species.gibbs(self.fixed_temperature, n_ratio, self.fixed_pressure*n_ratio)

    def gibbs_reactant(self, species):
        pass


def optimize_composition(problem: Problem, convergence:float):

    problem.Product_Mixture._build_species_dict()
    problem.Reactant_Mixtture._build_species_dict()
    product_Nj = problem.Product_Mixture.Nj()
    # static data used for matrix construction
    product_elements = [problem.Product_Mixture.get_elements_species(s.name) for s in problem.Product_Mixture.species]
    reactant_elements = [problem.Reactant_Mixture.get_elements_species(s.name) for s in problem.Reactant_Mixture.species]
    n_species = len(problem.Product_Mixture.species)
    global_product_elements = problem.Reactant_Mixtture.get_elements_global()

    # write wrapper to eval this from construct_matrices
    def e(delN: np.array, delNj: np.array, N:np.array, Nj:np.array):
        #underrelaxation parameter
        #Gordon and McBride 1994
    
        e1 = 2 / np.max(5, np.abs(delN), np.abs(delNj))
        e2 = np.min(
            np.abs(
                -np.log(Nj/N) - 9.2103404\
                  /  (delNj - delN)
            )
        )
        return np.min(1, e1, e2)

    def mk_mat(x):

        @lru_cache()
        def A1Mat():
            # identity matrix
            return np.identity(n_species)

        @lru_cache()
        def A2Mat():
            # -a1 -a21 -a32 ... -1 for element count for each species
            return np.block([np.array([
                np.array([0 if k not in list(product_elements[s]) \
                    else -product_elements[s][k] for k, v in global_product_elements.keys()])
                for s in product_elements.keys()
            ]), np.full((n_species, 1), -1)])

        @lru_cache()
        def A3Mat(A2):
            # N1a11 N2a12 N3a13 ...
            Nj = problem.Product_Mixture.Nj()
            mol_row_vec = np.array([
                v for _, v in Nj.items()
            ])  
            return mol_row_vec * A2[:, 0:-1]
            
        @lru_cache()
        def A4Mat():
            # 000 -- 000 -N
            pass


        A2 = A2Mat()
        A = np.block([
            [A1Mat(), A2],
            [A3Mat(A2), A4Mat()]
        ])


        def b():
            gibbs = np.array([
                problem.gibbs_product(s)
                for s in problem.Product_Mixture.SDict.keys()
            ])
            #b0 = 
            
        b = b()

        return A, b

    def update():
        pass

    CONVERGENCE_FLAG = False
    while not CONVERGENCE_FLAG:
        pass
