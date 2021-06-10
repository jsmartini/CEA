from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import numpy as np
from enum import Enum
from textwrap import dedent, wrap
from progress.bar import Bar
import pickle
from Database import *
import re

def ParseThermoINP(fpath: str, database: ThermoDatabase):
    NEW_SPECIES_FLAG = True
    INFO_LINE_FLAG = False
    NEW_INTERVAL_FLAG = False
    TITLE_FLAGS = ["!", "thermo"]
    INTERVAL_CNT = 0
    TARGET_INTERVAL_LINE_CNT = 3 # 3 lines for each interval entry
    PRODUCT_FLAG = ["END PRODUCTS", True]
    REACTANT_FLAG = ["END REACTANTS", False]  
    File_info = [] 
    database.Products = {}
    database.Reactants = {}
    # place holders
    target_species = Species()
    target_species.intervals = []
    target_interval = Interval()
    target_interval.cp_coefficients = []
    target_interval.t_Exponents = []
    target_interval.b_constants = []
    length = 0
    interval1regex = re.compile(r".+\d\d")
    with open(fpath, "r") as f:
        while (l:=f.readline()): length+=1
    with open(fpath, "r") as inp:

        status = Bar("Loading Thermochemical Data", max = length)

        for cnt, line in enumerate(inp):
            status.next()
            #print(f"{cnt}: {line}")
            if cnt < 42:        # parses lines 0-42 for file info
                File_info.append(line)
                continue
            
            elif PRODUCT_FLAG[0] in line:
                PRODUCT_FLAG[1] = False
                continue
            #start a new entry
            elif NEW_SPECIES_FLAG:    # if new species parse
                #print(f"species meta data {cnt}")
                target_species = Species()
                target_species.intervals = []
                NEW_SPECIES_FLAG = False
                INFO_LINE_FLAG = True   # move on to the info line
                tokens_ = line.split()
                target_species.name = tokens_[0]
                target_species.info = "".join(tokens_[0:-1])
                continue
                #go to next line to grab species entry meta-data and phyiscal data
            elif INFO_LINE_FLAG:
                #print(f"info {cnt}")
                INFO_LINE_FLAG = False
                INTERVAL_CNT = int(line[0:3])   # total number of intervals
                target_species.interval_cnt = INTERVAL_CNT
                target_species.source = line[3:10]
                offset = 2
                start = 10
                elements = {}
                while start < 50:
                    elements[line[start:start+2]] = int(float(line[start+4:start+8]))
                    start = start+8
                    #  
                    # 2 tpis96 AL  1.00H   1.00F   2.00    0.00    0.00 0   65.9862844    -765299.182
                    # how the loop progresses through the info, ends on col 51 ( index 50 zero start)
                    #target_species.elements[line[10:12]] = int(line[14:18])
                    #target_species.elements[line[18:20]] = int(line[22:24])
                target_species.elements = elements
                target_species.phase = float(line[51:53])
                target_species.molecular_weight = float(line[55:66])
                target_species.heat_of_formation = float(line[66:])
                continue

            elif INTERVAL_CNT > 0:
                #print(f"interval {INTERVAL_CNT}")
                """
                Logs interval
                Ref.
                index.
                0   200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        12299.975
                0 1.050181431D+05-1.398654202D+03 9.166334840D+00 3.867286870D-03-4.830841090D-06
                0 3.099287181D-09-8.440702990D-13                -8.659048800D+04-2.557299683D+01
                0   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        12299.975
                """

                if TARGET_INTERVAL_LINE_CNT > 0:
                    #print(f"target_interval {TARGET_INTERVAL_LINE_CNT}")
                    if TARGET_INTERVAL_LINE_CNT == 3:
                        # temeprature range and t exp and enthalpy
                        tokens = line.split()
                        try:
                            target_interval.temperature_range = (
                                float(tokens[0]),
                                float(tokens[1])
                            )
                        except:
                            target_interval.temperature_range = (
                                float(tokens[0].replace("D", "E")),
                                float(tokens[1].replace("D", "E"))
                            )

                        target_interval.t_Exponents = tuple(map(float, tokens[2:-2]))
                        target_interval.enthalpy = tokens[-1]
                        TARGET_INTERVAL_LINE_CNT -= 1

                    elif TARGET_INTERVAL_LINE_CNT == 2:
                        #first 5 cp coe
                        offset = 16
                        index = 0
                        while index < len(line)-1:
                            #print(line[index:index+offset])
                            target_interval.cp_coefficients.append(float(line[index:index+offset].replace("D", "E")))
                            index += offset
                        TARGET_INTERVAL_LINE_CNT -= 1

                    elif TARGET_INTERVAL_LINE_CNT == 1:
                        line = "".join(line.split(" ")).strip("\n")
                        consts_list = []
                        offset = 15
                        index = 0
                        #while index < len(line):
                        #    if line[index:index+offset] == '':
                        #        break
                        #    print(line[index:index+offset])
                        #    input()
                        #    consts_list.append(float(line[index:index+offset].replace("D", "E")))
                        #    
                        #    index += offset + 1 
                        consts_list = re.split('(?<=D(\+|\-)\d\d)', line)
                        #print(consts_list)
                        # crap code to fix edge case bug

                        #I hate this so very much
                        #some version error caused this.
                        bad_entries = ['', " ", "+", "-"]
                        for idx, elm in enumerate(consts_list):
                            for bad in bad_entries:
                                if elm == bad:
                                     consts_list.remove(elm)
                                     #print(elm)
                        if consts_list[-1] == '':
                            consts_list = consts_list[:-1]
                        #print(consts_list)

                        cast = lambda x: float(x.replace("D", "E"))
                        consts_list = list(map(cast, consts_list))
                        target_interval.cp_coefficients.extend(consts_list[0:2])
                        target_interval.b_constants = consts_list[2:]
                        print(target_interval.b_constants)
                        target_species.intervals.append(target_interval)
                        target_interval = Interval()
                        target_interval.cp_coefficients = []
                        target_interval.t_Exponents = []
                        target_interval.b_constants = []
                        TARGET_INTERVAL_LINE_CNT = 3
                        INTERVAL_CNT -= 1   

            if INTERVAL_CNT == 0 and not NEW_SPECIES_FLAG and not INFO_LINE_FLAG:
                # reset and start new species
                NEW_SPECIES_FLAG = True
                INFO_LINE_FLAG = False
                NEW_INTERVAL_FLAG = False
                
                if PRODUCT_FLAG[1]: 
                    target_species.species_type = chemical_category.PRODUCT
                    database.Products[target_species.name] = target_species
                if not PRODUCT_FLAG[1]: 
                    target_species.species_type = chemical_category.REACTANT
                    database.Reactants[target_species.name] = target_species
                target_species = Species()
                target_species.intervals = []

        else:
            # turn array of text into paragraph maintaining newline formatting
            File_info = "\n".join(File_info)
            status.finish()
            database.ThermoINPinfo = File_info
            return database

def ParseTransINP(fpath:str, database: ThermoDatabase):
    length = 0
    database.Transport = {}

    def new_entry() -> TransportSpecies:
        t = TransportSpecies()
        t.secondary = ""
        t.V = []
        t.C = []
        return t

    def new_data() -> TransportData:
        t = TransportData()
        t.temperature_range = ()
        t.coefficients = []
        return t

    with open(fpath, "r") as f:
        while (l:=f.readline()): length+=1

    with open(fpath, 'r') as trans:

        status = Bar("Loading Transport Data", max = length)
        target_species = new_entry()
        for cnt, line in enumerate(trans):
            status.next()
            
            if "." not in line:
                try:
                    if target_species.name not in database.Transport:
                        database.Transport[target_species.name] = []
                    database.Transport[target_species.name].append(target_species)
                except:
                    pass
                target_species = new_entry()
                tokens_ = line.split()
                target_species.name = tokens_[0]
                if len(tokens_) > 4:
                    target_species.secondary = tokens_[1]
                    target_species.info = " ".join(tokens_[2:])
                else:
                    target_species.info = " ".join(tokens_[1:])
            else:
                identifier = line[0:3]
                temps = tuple(map(float, line[3:20].split()))
                consts = line[20:]
                index = 0
                offset = 14
                coefficients = []
                while index < len(consts)-1:
                    if consts[index:index+offset] == '':
                        break
                    try:
                        const = float(consts[index:index+offset].replace(" ", ""))
                    except:
                        const = float(consts[index:index+offset])
                        consts.append(const)
                    index += offset + 1
                target_transport_data = new_data()
                target_transport_data.temperature_range = temps
                target_transport_data.coefficients = consts
                if "V" == identifier:
                    target_species.V.append(target_transport_data)
                elif "C" == identifier:
                    target_species.append(target_transport_data)    
        else:
            status.finish()
            return database


def load_database_from_inp(thermo_fpath:str, trans_fpath:str):
    db = ThermoDatabase()
    ParseThermoINP(thermo_fpath, db)
    ParseTransINP(trans_fpath, db)
    return db

def load_database_from_pickle(db_fpath:str):
    with open(db_fpath, "r") as f:
        return pickle.load(f)

if __name__ == "__main__":
    database = ThermoDatabase()
    ParseThermoINP("thermo.inp", database)
    ParseTransINP("trans.inp", database)
    pickle.dump(database, open("database.bin", "wb"))
    from IPython import embed
    embed()