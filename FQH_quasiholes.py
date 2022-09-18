import FQH_states as FQH
import numpy as np
from itertools import combinations, compress, product
import misc

def raised(partition, index_set):
	return np.array([partition[i]+1 if i in index_set else partition[i] for i in range(len(partition))])

def is_valid(partition, No):
	if partition[-1] >= No:
		return False # Assume partition is sorted
	elif any((partition[1:]-partition[:-1])==0):
		return False
	else:
		return True

def all_valid_raises(partition, k, No):
	Ne = len(partition)
	all_raises = [raised(partition, index_set) for index_set in combinations(range(Ne), k)]
	valid_check = [is_valid(thing, No) for thing in all_raises]
	return [misc.index_to_binary(thing, No, get_string=True) for thing in compress(all_raises, valid_check)]

def FQH_quasihole_poly(Ne, pos, ground_state):
	# note that ground state must also be a pure polynomial
	q = len(pos)
	assert q>=0
	if len(pos)==0:
		return ground_state
	else:
		state_rec = FQH_quasihole_poly(Ne, pos[:-1], ground_state)
		w = pos[-1]
		if w==0:
			new_basis = ["0"+thing for thing in state_rec.basis]
			return FQH.fqh_state((new_basis, state_rec.coef))
		else:
			state_bas = state_rec.get_basis_index()
			dim = len(state_bas)

			state = FQH.fqh_state()
			for (i,j) in product(range(dim), range(Ne+1)):
				basis = all_valid_raises(state_bas[i], Ne-j, 3*Ne-2+q)
				coef  = np.ones(len(basis)) * ((-w/2)**j) * state_rec.coef[i]

				state += FQH.fqh_state((basis, coef))

			return FQH.prune(state)

if __name__=="__main__":
	import time
	from two_qh_state import two_qh_IQH_state, two_qh_IQH_state_false

	Ne=6
	gs = FQH.fqh_state("jack_6e")

	st = time.time()
	s1 = FQH_quasihole_poly(Ne, [2],gs)
	s1.disk_normalize()
	print(f"One quasihole state: {time.time()-st} seconds")
	s1.plot_disk_density("FQH_state_density/Laughlin_1h.pdf")
	s1.printwf("FQH_state_density/Laughlin_1h")

	print("---")

	st = time.time()
	s2 = FQH_quasihole_poly(Ne, [-2,2],gs)
	s2.disk_normalize()
	print(f"Two quasihole state: {time.time()-st} seconds")
	s2.plot_disk_density("IQH_state_density/Laughlin_2h.pdf")
	s2.printwf("FQH_state_density/Laughlin_2h")

	print("---")
	"""
	st = time.time()
	pos = [-1+2j,2,-1-2j]
	s3 = IQH_quasihole_poly(Ne,pos)
	s3.disk_normalize()
	print(f"Three quasihole state: {time.time()-st} seconds")
	s3.plot_disk_density("IQH_state_density/IQH_3h.pdf")
	s3.printwf("IQH_state_density/IQH_3h")


	print("---")
	
	st = time.time()
	pos = [0]
	s3 = FQH_quasihole_poly(6,pos, gs)
	s3.disk_normalize()
	print(f"One quasihole state: {time.time()-st} seconds")
	s3.printwf("FQH_state_density/test")
	s3.plot_disk_density("FQH_state_density/test.pdf")
	"""
