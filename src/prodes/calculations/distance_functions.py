import math

def atom_charge_coulomb(charge):
	return charge * 1.6e-19

def charge_simple(charge, distance, dielectric_constant):
	"""takes the charge in coulomb, distance in meters and dielectric constant to calculate the potential"""

	absolute_permittivity = 8.854e-12
	permittivity = dielectric_constant*absolute_permittivity

	return charge/(permittivity*distance*4*math.pi) 

def potential_multiple_media(charge : float, distances : list, dielectic_constants: list):
	"""takes the charge in coulomb, distances in meters and dielectric constants to calculate the potential at a point in space"""
	
	absolute_permittivity = 8.854e-12
	denominator = 0
	for distance, dielectric_constant in zip(distances, dielectic_constants):
		permittivity = dielectric_constant*absolute_permittivity
		denominator += permittivity*distance

	return (charge/(denominator*4*math.pi))