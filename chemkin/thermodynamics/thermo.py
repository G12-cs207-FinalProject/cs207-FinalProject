import os.path
import sqlite3


class ThermoDAO():
    """ Database Access Object for thermodynanmics
    """

    def __init__ (self, db_name):
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        self.db_path = os.path.join(BASE_DIR, db_name)

    def get_coeffs (self, species_name, temp_range):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7
                    FROM {} 
                    WHERE SPECIES_NAME = "{}"'''.format(temp_range.upper(), species_name)
        coeffs = list(cursor.execute(query).fetchall()[0])
        db.commit()
        db.close()
        return coeffs

    def get_species(self, temp, temp_range):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        if temp_range == 'low': # temp_range == 'low'
            query = '''SELECT SPECIES_NAME FROM {} WHERE TLOW < {}'''.format(temp_range.upper(), temp)
        else: # temp_range == 'high'
            query = '''SELECT SPECIES_NAME FROM {} WHERE THIGH > {}'''.format(temp_range.upper(), temp)
        species = []
        for s in cursor.execute(query).fetchall():
            species.append(s[0])
        db.commit()
        db.close()
        return species

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)