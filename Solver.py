from functools import lru_cache
from Database import *
from Parser import *
import numpy as np
from utils import *

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

class Solver:
    """
        Free Minimization of Gibbs Energy solver for equilibrium
    using multivariate newton raphson method using an underrelaxation parameter
    for updates (Gordon and McBride 1994)
    """

    def __init__(self, problem: Problem):
        self.problem = problem
        self.Products = self.problem.Product_Mixture
        self.Reactants = self.problem.Reactant_Mixture
        self.product_species_cnt = len(self.Products.species)
        self.product_elements = {s.name:self.Products.get_elements_species(s.name) for s in self.Products.species}
        self.global_product_elements = self.Reactants.get_elements_global()
        self.i_len = 0
        self.j_len = 0

    def _A1Block(self):
        mat =  np.identity(self.product_species_cnt)
        #matprint(mat)
        return mat

    def _A2Block(self):
        mat =  np.block([np.array([
                np.array([0 if k not in list(self.product_elements[s].keys()) \
                    else -self.product_elements[s][k] for k, v in self.global_product_elements.items()])
                for s in self.product_elements.keys()
            ]), np.full((self.product_species_cnt, 1), -1)])
        self.i_len = mat.shape[0]
        self.j_len = mat.shape[1]
        #matprint(mat)
        return mat

    def _A3Block(self, A2):
        # N1a11 N2a12 N3a13 ...
            #print(A2.shape)
            Nj = self.Products._Nj()
            print(Nj)
            mol_row_vec = np.array([[
                v for _, v in Nj.items()
            ] for _ in range(A2.shape[1]-1)
            ])  
            #print(mol_row_vec)
            
            mat =  np.multiply(mol_row_vec, A2[:, 0:-1].T)
            mol_row_vec = np.array([
                v for _, v in Nj.items()
            
            ])
            mat = np.block([
                [mat],
                [mol_row_vec]
            ])
            #matprint(mat)
            return mat

    def _A4Block(self):
        mat = np.zeros(shape=(self.i_len, self.j_len))
        mat[self.i_len-1, self.j_len-1] = -self.Products.N
        #matprint(mat)
        return mat

    def _AMatrix(self):
        A2 = self._A2Block()
        return np.block([
            [self._A1Block(), A2,],
            [self._A3Block(A2), self._A4Block()]
        ])

    def _BCol(self):
        pass

    def _update(self, del_x):
        pass

    def run(self):
        pass