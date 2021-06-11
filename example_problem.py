from Database import *
from Parser import *

global DB
DB = load_database_from_inp()

Reactants = ["O2(L)", "CH4(L)"]
Products = ["CO2", "H2O", "CO"]

eg = Problem(
    fixed_pressure = 6, # bar
    fixed_temperature = 4000, # kelvin
    Product_Mixture = SpeciesMixture(
        species = [DB.Products[i] for i in Products]
    ), 
    Reactant_Mixture = SpeciesMixture(
        species = [DB.Reactants[i] for i in Reactants]
    )
)

eg.Product_Mixture._build_species_dict()
eg.Reactant_Mixture._build_species_dict()
