###############################################################################
# Tests for chemkin.thermodynamics.thermo module
###############################################################################

from chemkin.thermodynamics.thermo import ThermoDAO

def test_get_coeffs():
	dao = ThermoDAO('NASA_coef.sqlite')
	high_coefs = dao.get_coeffs('H', 'high')
	low_coefs = dao.get_coeffs('H', 'low')
	assert len(high_coefs) == 7
	assert len(low_coefs) == 7


def test_get_species():
	dao = ThermoDAO('NASA_coef.sqlite')
	species_high = dao.get_species(1001, 'high')
	species_low = dao.get_species(999, 'low')
	assert len(species_high) > 0
	assert len(species_low) > 0
