from monodromy import charge_monodromy, flavor_charge_monodromy
import itertools

### How many times we expect an MS wall to jump from a zone to another one
CHARGE_ORBIT_DEPTH = 3

class Zone:
	"""
	The u-plane is divided into zones, these are vertical strips
	between naighboring pairs of branch points, and delimited
	by branch cuts.
	The lower boundary of the zone is a horizontal segment 
	running between the cuts and placed at the height corresponding
	to the imaginary part of the higher branch point.
	The zone between branch points i and i+1 is zone i.
	If we have n branch points, zone n-1 will thus be the last,
	starting from the right of the last branch-cut and including 
	the rest of the u-plane, up to the left of the 0-th branch-cut.

	An example:

			|		|		|		|
			|		|		|		|
		3	|	0	|	1	|	2	|	3
			|		|		|		|
			|		|		|		|
			|-------x-------|		|
			x				|		|
							|-------x
					3		x			

	"""
	
	def __init__(self, number, bp_left, bp_right, is_last_zone=False):
		self.number = number
		self.bp_left = bp_left
		self.bp_right = bp_right

		self.is_last_zone = is_last_zone

		if self.is_last_zone == False:
			self.left_boundary = bp_left.locus.real
			self.right_boundary = bp_right.locus.real
			self.bottom_boundary = min(
									[bp_left.locus.imag, bp_right.locus.imag]
									)

		self.left_neighbor = None
		self.right_neigbor = None
		
		### Will be employed for zones 0,1, ..., n-2
		self.bottom_neighbor = None

		### Will be employed for zone n-1
		self.upper_neighbors = None


	def set_neighbors(self, neighbors):
		"""
		Neigbors of a zone will be used to determine the 
		potential orbit of the charges at an intersection point.
		Notice that from zone n-1, one may hop into any of the 
		other zones without undergoing monodromy.
		However, one may also hop into zones 0, n-2 by undergoing
		a monodromy. That is why the last zone (n-1) has both 
		a left neighbor (n-2), a right neighbor (0), and a set of
		upper_neighbors, which includes all other zones, including 0, n-2.
		When going from n-1 -> 0 or n-2 as right/left neighbors we 
		will take monodromy into account.
		When going from n-1 -> upper_neighbors there will be no monodromy,
		because it means that we move vertically to the neigboring zone.
		"""
		if self.is_last_zone==False:
			left_neighbor, bottom_neighbor, right_neighbor = neighbors
			self.left_neighbor = left_neighbor
			self.right_neighbor = right_neighbor
			self.bottom_neighbor = bottom_neighbor
		
		elif self.is_last_zone==True:
			### By convention, the left neighbor of the last zone
			### is zone n-2, while the right neigbor is zone 0
			left_neighbor, upper_neighbors, right_neighbor = neighbors
			self.left_neighbor = left_neighbor
			self.right_neighbor = right_neighbor
			self.upper_neighbors = upper_neighbors

	def contains_point(self, u):
		if self.is_last_zone == False:
			if self.left_boundary < u.real < self.right_boundary and \
					u.imag > self.bottom_boundary:
				return True

			else:
				return False

		elif self.is_last_zone == True:
			raise ValueError('To find out if a point belongs to the last zone,'
							+' you have to check that it doesnt belong to any '
							+ 'other zone.')



def build_charge_orbit(intersection_point, zones):
	initial_zone = determine_zone(intersection_point.locus, zones)

	current_zone = initial_zone
	current_gauge_charges = intersection_point.gauge_charges
	current_flavor_charges = intersection_point.flavor_charges

	return charges_in_neighboring_zones(
											CHARGE_ORBIT_DEPTH,
											current_zone, 
											current_gauge_charges, 
											current_flavor_charges,
											zones
										)


def charges_in_neighboring_zones(depth, current_zone, gauge_charges, 
														flavor_charges, zones):
	depth_1 = depth - 1
	
	### Given the fgauge charges [g_1, g_2] of the intersection
	### i.e. those of the two parent k-walls, we compute their
	### parallel transport across branch cuts to the left and the right
	### of the current zone.
	### We do the same for the flavor charges.

	gauge_charges_to_the_left = \
		[
			charge_monodromy(gauge_charges[i], current_zone.bp_left, 'ccw') \
			for i in [0, 1]
		]
		
	flavor_charges_to_the_left = \
		[
			flavor_charge_monodromy(
									gauge_charges[i], flavor_charges[i], 
									current_zone.bp_left, 'ccw'
									) \
			for i in [0, 1]
		]

	gauge_charges_to_the_right = \
		[
			charge_monodromy(gauge_charges[i], current_zone.bp_right, 'cw') \
			for i in [0, 1]
		]

	flavor_charges_to_the_right = \
		[
			flavor_charge_monodromy(
									gauge_charges[i], flavor_charges[i],
									current_zone.bp_right, 'cw'
									) \
			for i in [0, 1]
		]

	new_charges = [
						[
							'gauge charges' , gauge_charges_to_the_left, 
							'flavor charges' , flavor_charges_to_the_left
						],
						[
							'gauge charges' , gauge_charges_to_the_right, 
							'flavor charges' , flavor_charges_to_the_right
						],
					]

	if depth == 0:
		return new_charges
	
	else:
		### return these 'neighboring charges', 
		### as well as the 'neighboring charges' 
		### of the neighboring zones (i.e. the)
		### neigbors of 2nd degree.

		if current_zone.is_last_zone == False:
			all_new_charges = new_charges \
							+ charges_in_neighboring_zones(
												depth_1, 
												current_zone.left_neighbor,
												gauge_charges_to_the_left,
												flavor_charges_to_the_left,
												zones
												) \
							+ charges_in_neighboring_zones(
												depth_1, 
												current_zone.right_neighbor,
												gauge_charges_to_the_right,
												flavor_charges_to_the_right,
												zones
												) \
							+ charges_in_neighboring_zones(
												depth_1, 
												current_zone.bottom_neighbor,
												gauge_charges,
												flavor_charges,
												zones
												)
		
		elif current_zone.is_last_zone == True:
			all_new_charges =  new_charges \
							+ charges_in_neighboring_zones(
												depth_1, 
												current_zone.left_neighbor,
												gauge_charges_to_the_left,
												flavor_charges_to_the_left,
												zones
												) \
							+ charges_in_neighboring_zones(
												depth_1, 
												current_zone.right_neighbor,
												gauge_charges_to_the_right,
												flavor_charges_to_the_right,
												zones
												)
			
			for upper_zone in current_zone.upper_neighbors:
				all_new_charges += charges_in_neighboring_zones(
												depth_1, 
												upper_zone,
												gauge_charges,
												flavor_charges,
												zones
												)

		### Now we remove the duplicates
		all_new_charges.sort()
		unique_new_charges = \
							list(all_new_charges for all_new_charges,_ \
							in itertools.groupby(all_new_charges) \
							)
		
		return unique_new_charges



def determine_zone(u, zones):
	found_zone = False
	for z in zones[:-1]:
		if z.contains_point(u) == True:
			found_zone = True
			return z

	if found_zone == False:
		### If it doesn't belong to zones 0,..,n-2 then it belongs 
		### to the last one.
		return zones[-1]
	

def orbits_coincide(charge_orbit_1, charge_orbit_2):
	answer = False
	for x in charge_orbit_1:
		for y in charge_orbit_2:
			if x == y:
				answer = True

	return answer

def orbit_is_contained(orbit_list, orbit):
	answer = False
	for i, orbit_i in enumerate(orbit_list):
		if orbits_coincide(orbit, orbit_i):
			answer = i
			break

	return answer




