import sage
from sage.rings.ring import Algebra
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import AlgebraElement
from sage.combinat.free_module import CombinatorialFreeModuleElement

class KLRAlgebraElement(CombinatorialFreeModuleElement):
	
	def __init__(self, *args, **kwds):
		
		CombinatorialFreeModuleElement.__init__(self, *args, **kwds)
		
		P = self.parent()
		self._n = P.positive_root().height()
		self.rename(self._get_print_name())
	
	def _get_print_name(self):
		mcs = self.monomial_coefficients()
		if len(mcs) != 0:
			new_name = ""
			for item in mcs.keys():
				X = ""
				for i in range(1, self._n+1):
					if item[0][i-1] == 1:
						X += "x_{%s}"%i
					elif item[0][i-1] != 0:
						X += "x_{%s}^{%s}"%(i, item[0][i-1])
				if item[1].is_one():
					T = ""
				else:
					T = "t_{%s}"%item[1]
				E = "e_{%s}"%item[2]
				coeff = mcs[item]
				if coeff == 1:
					new_name += "%s%s%s + "%(X,T,E)
				else:
					new_name += "%s*%s%s%s + "%(coeff,X,T,E)
			new_name = new_name[:-3]
			return new_name
		return "0"
	
	def _tikz_monomial(self, x_sep=0.5, y_sep=0.3, coeff=1, first_term=True, coeff_height=None):
		rw = self._get_reduced_word()
		iv = self.get_integer_vector()
		perm = self.get_permutation()
		n = self._n
		l = len(rw)
		max_dots = max(iv)
		dot_space = ceil((max_dots+2)/3.0)*y_sep
		if coeff_height == None:
			coeff_height = (y_sep*(l+2)+dot_space)/2
		front_matter = "\\begin{tikzpicture}\n\t\\foreach \\x in {1, ..., %s}\n\t\t\\node (\\x_-1) at (%s*\\x, %s){};\n\t\\foreach \\x in {1, ..., %s}\n\t\t\\foreach \\y in {0, ..., %s}\n\t\t\t\\node (\\x_\\y) at (%s*\\x, %s*\\y){};\n"%(n, x_sep, -dot_space, n, l+2, x_sep, y_sep)
		if coeff == 1:
			coeff = ""
		else:
			coeff = "\\,\\,\\left(%s\\right)\\,\\,"%latex(coeff)
			
			
		if first_term:
			end_matter = "\t\\node[on grid, above=%s of 1_-1] (temp) {};\n\t\\node[left=0 of temp]{$%s$};\n\\end{tikzpicture}"%(coeff_height, coeff)
		else:
			end_matter = "\t\\node[on grid, above=%s of 1_-1] (temp) {};\n\t\\node[left=0 of temp]{$+%s$};\n\\end{tikzpicture}"%(coeff_height, coeff)
		
		
		
		simple_roots = self.parent().root_system().root_lattice().simple_roots()
		
		
		
#			color_list = rainbow(len(simple_roots))
#			colors = {i:color_list[i-1] for i in range(len(simple_roots))}
		
		
		color_list = rainbow(len(set(perm)))
		colors = {list(set(perm))[i]:color_list[i] for i in range(len(list(set(perm))))}
		
		
		color_string = ""
		for i in colors:
			color_string += "\t\\definecolor{color%s}{HTML}{%s}\n"%(i, colors[i][1:].upper())
		
		strands = {i:"\t\\draw [color=color%s] (%s_-1.center)--(%s_0.center)--(%s_1.center)"%(perm[i-1], i, i, i) for i in range(1,n+1)}
		for i in range(l):
			strands[rw[i]], strands[rw[i]+1] = (strands[rw[i]+1]+"--(%s_%s.center)"%(rw[i], i+2), strands[rw[i]]+"--(%s_%s.center)"%(rw[i]+1, i+2))
			for j in range(1, n+1):
				if j != rw[i] and j != rw[i]+1:
					strands[j] += "--(%s_%s.center)"%(j, i+2)
		for j in range(1, n+1):
			strands[j] += "--(%s_%s.center);\n"%(j, l+2)
		
		dots = {i:"\t\\foreach \\y in {1, ..., %s}\n\t\t\\node[on grid, circle, fill, inner sep=1pt, above=%s/%s*\\y of %s_-1]{};\n"%(iv[i-1], dot_space, iv[i-1]+1, i) for i in range(1, n+1)}
		
		
		the_string = front_matter+color_string
		for i in strands:
			the_string += strands[i]
		for i in dots:
			if iv[i-1] != 0:
				the_string += dots[i]
		the_string += end_matter
		return the_string
		
		
		
		
	def tikz(self, x_sep=0.5, y_sep=0.3):
		the_string = ""
		mcs = self.monomial_coefficients()
		first_term=True
		min_height = Infinity
		for monomial in mcs:
			my_mon = self.parent().monomial(monomial)
			rw = my_mon._get_reduced_word()
			iv = my_mon.get_integer_vector()
			l = len(rw)
			max_dots = max(iv)
			dot_space = ceil((max_dots+2)/3.0)*y_sep
			this_height = (y_sep*(l+2)+dot_space)/2
			if this_height < min_height:
				min_height = this_height
		the_terms = []
		for monomial in mcs:
			the_terms += ["\n"+self.parent().monomial(monomial)._tikz_monomial(x_sep=x_sep, y_sep=y_sep, coeff=mcs[monomial], first_term=first_term, coeff_height=min_height)]
			first_term=False
		for term in the_terms:
			the_string += term
		return the_string[1:]

	def _get_reduced_word(self):
		my_parent = self.parent()
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return my_parent.reduced_word(self.get_Sn_element())
		
#	def __repr__(self):
#		return "HAHAHA GOTCHA"
	
	def is_monomial(self):
		if len(self.monomial_coefficients()) == 1:
			return True
		return False
	
	# determines the idempotent e
	def get_permutation(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return self.leading_support()[2]
		
	# determines the t
	def get_Sn_element(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return self.leading_support()[1]
	
	# determines the X
	def get_integer_vector(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return self.leading_support()[0]
		
	def get_first_nonzero_index(self):
		v = self.get_integer_vector()
		for i in range(len(v)):
			if v[i] != 0:
				return i
		return -1

	def is_x(self):
		for i in range(1, self._n+1):
			if self == self.parent().x(i):
				return i
		return False		
	
	def is_t(self):
		for i in range(1, self._n):
			if self == self.parent().t(i):
				return i
		return False
	
	def is_e(self):
		for p in self.parent()._P:
			if self == self.parent().e(p):
				return p
		return False
	
	""" DO I NEED THIS? NOT FINISHED...
	def x_part(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		mcs = self.monomial_coefficients()
	"""

class KLRAlgebra(CombinatorialFreeModule):
	r"""
	Create the given KLR algebra.
	
	INPUT:
	
	- "base_ring", a ring
	- "root_system", a root system (see RootSystem)
	- "positive_root", an element of the positive root lattice of root_system
	- "reduced_words", a function whose inputs are elements of the symmetric group S_n (where n is the height of the positive root) and whose outputs are a distinguished choice of reduced expressions for such an S_n element. The default is "lambda g: g.reduced_word()"
	- "signs", a function which takes integers i and j as inputs and returns +- 1 with the stipulation that the image of (i, j) is opposite the image of (j, i). This function need not be defined for all pairs of integers, only for pairs integers whose corresponding entry in the Cartan matrix of root_system is negative
	
	OUTPUT:
	
	- a "KLRAlgebra" instance
	
	"""
	
	Element = KLRAlgebraElement
	
	def __init__(self, base_ring, root_system, positive_root, reduced_words=None, signs=None, **kwds):
		r"""
		See "KLRAlgebra" for full documentation.
		"""
		self._root_system = root_system
		
		
		# Get the cartan matrix. If we are not of finite type, this won't work so I have
		# to make a lazy family. The only supported infinite type is ["A", ZZ] right now.
		# To add support for other infinite types, one needs only to define the proper
		# function giving the appropriate entry in the infinite-dimensional cartan matrix
		# for that infinite-type.
		try:
			self._cartan_matrix = root_system.cartan_matrix()
		except AttributeError as e:
			if root_system.cartan_type() == CartanType(["A", ZZ]):
				self._cartan_matrix = Family(ZZ, lambda i: Family(ZZ, lambda j: 2*kronecker_delta(i, j) - kronecker_delta(i-1, j) - kronecker_delta(i, j-1)))
			else:
				raise e
		
		self._positive_root = positive_root
		self._n = positive_root.height()
		self._Sn = SymmetricGroup(self._n)
		self._name = "KLR Algebra over %s for %s with positive root %s"%(repr(base_ring), repr(self._root_system), repr(self._positive_root))
		if reduced_words != None:
			self._name += " with custom reduced words"
			if signs != None:
				self._name += " and sign conventions"
		elif signs != None:
			self._name += " with custom sign conventions"
		
		#self._reduced_words = reduced_words or {g:g.reduced_word() for g in self._Sn}
		
		if reduced_words == None:
			self._custom_words = False
			temp_reduced_words = lambda g: g.reduced_word()
			self._reduced_words = lambda g: temp_reduced_words(self._Sn(g))
		else:
			self._custom_words = True
			self._reduced_words = reduced_words
		
		if signs == None:
			self._custom_signs = False
			self._signs = lambda i, j: sign(j-i)
		else:
			self._custom_sings = True
			self._signs = signs
		
		from sage.combinat.integer_vector import IntegerVectors_k
		self._IVn = IntegerVectors_k(self._n)
		
		IndexList = []
		mcs = positive_root.monomial_coefficients()
		for key in mcs:
			IndexList += mcs[key]*[key]
					
		self._P = Permutations(IndexList)
		
		self._CP = cartesian_product([self._IVn, self._Sn, self._P])
		
		CombinatorialFreeModule.__init__(self, base_ring, self._CP, category=GradedAlgebrasWithBasis(base_ring), **kwds)
		
	def _repr_(self):
		r"""
		Return the string representation of self.
		"""
		return self._name
	
	#def base_ring(self):
	#	return self.base()
	
	def change_ring(self, ring=None):
		r"""
		Return the KLR algebra identical to this one, ove the specified ring
		"""
		if ring == None:
			ring = self.base_ring()
		if self._custom_words:
			reduced_words = self._reduced_words
		else:
			reduced_words = None
		if self._custom_signs:
			signs = self._signs
		else:
			signs = None
		return KLRAlgebra(ring, self._root_system, self._positive_root, reduced_words=reduced_words, signs=signs)
	
	def characteristic(self):
		r"""
		Return the characteristic of this KLR Algebra.
		"""
		return self.base_ring().characteristic()
		
	def is_finite(self):
		r"""
		Return whether or not this KLR algebra is finite. It isn't unless its base ring is finite and its associated positive root is zero.
		"""
		return self.base_ring().is_finite() and self._positive_root == 0
				
	def root_system(self):
		r"""
		Return the root system associated to this KLR Algebra.
		"""
		return self._root_system
	
	def cartan_matrix(self):
		r"""
		Return the cartan matrix of the root system associated to this KLR Algebra.
		"""
		return self._cartan_matrix
	
	def positive_root(self):
		r"""
		Return the positive root associated to this KLR Algebra (as an element of the root lattice of the root system associated to this KLR Algebra).
		"""
		return self._positive_root
	
	def symmetric_group(self):
		r"""
		Return the symmetric group associated to this KLR Algebra. This is the symmetric group S_n where n is the height of the positive root associated to this KLR Algebra.
		"""
		return self._Sn
		
	def reduced_word(self, g=None):
		r"""
		Return the distinguished reduced word for g.
		
		INPUT:
		
		- "g", an element of the symmetric group associated to this KLR Algebra.
		
		OUTPUT:
		
		- If g != None, return the distinguished reduced word of g.
		- If g == None, return the function which yields the distinguished reduced words.
		"""
		if g == None:
			return self._reduced_words
		if g not in self._Sn:
			raise ValueError("%s not an element of %s"%(g, self._Sn))
		return self._reduced_words(self._Sn(g))
	
	def signs(self, pair=None):
		r"""
		Return the sign associated to pair. 
		
		INPUT:
		
		- "pair", a pair of integers (i, j) whose corresponding entry in the cartan matrix associated to this KLR Algebra is negative.
		
		OUTPUT:
		
		- If pair != None, return the sign (1 or -1) associated to the pair (i, j). This should always return sign((i, j)) == -sign((j, i)).
		- If pair == None, return the function which yields the signs associated to pairs of integers.
		"""
		if pair == None:
			return self._signs
		return self._signs(pair[0], pair[1])
		
	def x(self, i=None):
		r"""
		Return the polynomial associated to i.
		
		INPUT:
		
		- "i", either an integer between 1 and n, or a length n integer vector.
		
		OUTPUT:
		
		- If i is an integer between 1 and n, returns the KLR algebra generator x_i.
		- If i is a length n integer vector, returns associated product of algebra generators.
		- If i == None, return a dictionary whose keys are integers 1 through n and values the associated polynomial generators of this KLR algebra.
		"""
		if i == None:
			return {j:self.x(j) for j in range(1, self._n+1)}
		if i not in range(1, self._n+1) and i not in self._IVn:
			raise ValueError("%s not an integer between 1 and %s, nor an element of %s"%(i, self._n, self._IVn))
		#B = self.basis()
		if i in self._IVn:
			v = self._IVn(i)
		else:
			v = self._IVn((0,)*(i-1)+(1,)+(0,)*(self._n-i))
		g = self._Sn(())
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
	
	def polynomial(self, v=None):
		r"""
		Deprecated. Use x instead
		"""
		if v == None:
			raise NotImplementedError("Must specify an integer vector for now")
		if v not in self._IVn:
			raise ValueError("%s not an appropriate integer vector"%v)
		v = self._IVn(v)
		g = self._Sn(())
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
		
	def t(self, g=None):
		r"""
		Return the crossing associated to g.
		
		INPUT:
		
		- "g", either an integer between 1 and n-1, or an element of the symmetric group S_n.
		
		OUTPUT:
		
		- If g is an integer between 1 and n-1, returns a crossing of the gth and (g+1)th strands.
		- If g is an element of S_n, returns a braid corresponding to the distinguished reduced word of g.
		- If g == None, returns a dictionary whose keys are the integers between 1 and n-1, and values are the associated crossings.
		"""
		if g == None:
			return {j:self.t(j) for j in range(1, self._n)}
		if g in range(1, self._n):
			g = self._Sn((g, g+1))
		elif g in self._Sn:
			g = self._Sn(g)
		else:
			raise ValueError("%s not an integer between 1 and %s, nor an element of %s"%(g, self._n-1, self._Sn))
		#B = self.basis()
		v = self._IVn((0,)*self._n)
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
	
	def e(self, p=None):
		r"""
		Return the idempotent corresponding to p
		
		INPUT:
		
		- "p", a permutation of length n whose entries correspond to the coefficients of the positive root associated to this KLR algebra.
		
		OUTPUT:
		
		- If p != None, returns the idempotent corresponding to p.
		- If p == None, returns a dictionary whose keys are permutations as above, and values are the corresponding idempotents.
		"""
		if p == None:
			return {tuple(j):self.e(j) for j in self._P}
		p = self._P(p)
		if p not in self._P:
			raise ValueError("%s not in %s"%(p, self._P))
		#B = self.basis()
		v = self._IVn((0,)*self._n)
		g = self._Sn(())
		p = self._P(p)
		return self.monomial(self._CP((v, g, p)))
		
	def one(self):
		r"""
		Returns the identity element of this KLR algebra. This is the sum of all idempotents corresponding to length n permutations whose entries correspond to the coefficients of the positive root associated with this KLR algebra.
		"""
		g = self._Sn.one()
		v = self._IVn((0,)*self._n)
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])

#	I'm not sure what this is supposed to do, or if it is necessary.		
#	def gen(self, i=0):
		
		
	def algebra_generators(self):
		r"""
		Returns the algebra generators. See x(), t(), e().
		"""
		return [self.x(i) for i in range(1, self._n+1)]+[self.t(i) for i in range(1, self._n)]+[self.e(p) for p in self._P]

	
	
	
	
	def _shuffle(self, word, p, Graph=None):
		r"""
		Returns an element of this KLR algebra which is equal to the product of the generators corresponding to the reduced word "word" and the idempotent determined by p
		
		INPUT:
		
		- "word", a list representing a reduced word for some S_n element
		- "p", a permutation corresponding to some idempotent e_p
		- "Graph", the reduced word graph containing the reduced word "word"
		
		OUTPUT:
		
		- An element of this KLR algebra (writting in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k]
		"""
		
		g = self._g_ify(word)
		good_word = self.reduced_word(g)
		
		if word == good_word:
			v = self._IVn((0,)*self._n)
			return self.monomial(self._CP((v, g, p)))
		
		if Graph == None:
			G = g.reduced_word_graph()
		else:
			G = Graph
		
		SP = G.shortest_path(tuple(word), tuple(good_word))
		i = 0
		try:
			while G.edge_label(SP[i], SP[i+1]) == 2:
				i += 1
			current_word = list(SP[i])
		except IndexError:
			v = self._IVn((0,)*self._n)
			return self.monomial(self._CP((v, g, p)))
			
			
			
		# Now there IS a next word, and the current word and next word differ by a braid relation.
		# Find the index at which this braid relation occurs. i is the index along SP at which 
		# current_word lies.
		next_word = list(SP[i+1])
		j = 0
		while current_word[j] == next_word[j]:
			j += 1
		
		# Now j is the first index at which current_word and next_word differ
		
		
		# end --> j+3, word --> current_word
		
		new_p_m = self._P(self._g_ify(current_word[j+3:])(tuple(p)))
		if current_word[j+1] == current_word[j+2] + 1:
			k = current_word[j+2]
			new = self._shuffle(current_word[:j]+[current_word[j+1], current_word[j+2], current_word[j+1]]+current_word[j+3:], p, Graph=G)
			c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
			if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
				new_g_m = self._Sn(())
				new_p_r = p
				new_p_l = new_p_m
				for r in range(-c_ij):
					new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
					new -= self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(current_word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(current_word[j+3:], new_p_r))
			return new
			
		if current_word[j+1] == current_word[j+2] - 1:
			k = current_word[j+1]
			new = self._shuffle(current_word[:j]+[current_word[j+1], current_word[j+2], current_word[j+1]]+current_word[j+3:], p, Graph=G)
			c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
			if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
				new_g_m = self._Sn(())
				new_p_r = p
				new_p_l = new_p_m
				for r in range(-c_ij):
					new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
					new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(current_word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(current_word[j+3:], new_p_r))
			return new
		
		
		
	
	
	def _reduce(self, word, p):
		r"""
		Returns an element of this KLR algebra which is equal to the product of the generators corresponding to "word" and the idempotent determined by "p"
		
		INPUT:
		
		- "word", any list containing positive integer entries between 1 and n
		- "p", a permutation corresponding to some idempotent e_p
		
		OUTPUT:
		
		- An element of this KLR algebra (writting in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k]
		"""
		# Base case: if my word is already reduced.
		# NOTE: THIS IS GOING TO GIVE A WRONG ANSWER! BECAUSE I HAVE NOT CHECKED TO SEE
		# IF I HAVE THE DISTINGUISED REDUCED WORD! I'M GOING TO PUT SOME EXTRA STUFF IN!
		g = self._g_ify(word)
		if g.length() == len(word):
			#v = self._IVn((0,)*self._n)
			#return self.monomial(self._CP((v, g, p)))
			
			#Uncomment this, and comment above for magic to happen.
			return self._shuffle(word, p)
	
	
		# I am looking for a minimal unreduced subword. Start at the beginning.
		i = 0
		# Increase the index until the subword starting at the beginning is unreduced.
		while self._g_ify(word[:i+1]).length() == i+1:
			i += 1
	
	
		# Now go backwards...
		j = i
		# Decrease this second index until the subword starting at j ending at i is unreduced.
		while self._g_ify(word[j:i+1]).length() == i+1-j:
			j -= 1
		
	
	
		# Another semi-base case: If my minimal unreduced subword has length 2 (and so it is a
		# double-crossing), apply the appropriate relation to reduce it.
		if i - j == 1:
			k = word[i]
			middle_p = self._P(self._g_ify(word[i+1:])(tuple(p)))
			if middle_p[k-1] == middle_p[k]:
				return self.zero()
			if self.cartan_matrix()[middle_p[k-1]-1][middle_p[k]-1] < 0:
				
				first_v = self._IVn((0,)*(k-1)+(-self.cartan_matrix()[middle_p[k-1]-1][middle_p[k]-1],)+(0,)*(self._n - k))
				first_g = self._Sn(())
				first_p = middle_p
				second_v = self._IVn((0,)*k+(-self.cartan_matrix()[middle_p[k]-1][middle_p[k-1]-1],)+(0,)*(self._n - k - 1))
				second_g = self._Sn(())
				second_p = middle_p
				
				first = self.monomial(self._CP((first_v, first_g, first_p)))
				second = self.monomial(self._CP((second_v, second_g, second_p)))
				
				middle = self.signs(pair=(middle_p[k-1]-1, middle_p[k]-1))*(first - second)
				
				new_p_l = self._P(self._g_ify(word[j:])(tuple(p)))
				new_p_r = p
				return self._reduce(word[:j], new_p_l) * (middle * self._reduce(word[i+1:], new_p_r))

			return self._reduce(word[:j]+word[i+1:], p)
	
	
	
		# Another semi-base case: If my minimal unreduced word can be reduced by swapping
		# distant transpositions, do such a swap.
		if abs(word[i] - word[i-1]) > 1:
			return self._reduce(word[:i-1]+[word[i], word[i-1]]+word[i+1:], p)
	
	
	
		# Now I have guaranteed that my minimal unreduced word ends in two adjacent
		# transpositions, and a consequence of minimality is that the two strands which appear
		# in the last two transpositions but don't cross there, DO cross somewhere else in the
		# minimal unreduced word, before position 0. Straighten this piece.	
		
		return self._straighten(word, j+1, i+1, p)


	def _straighten(self, word, start, end, p):
		r"""
		Applies appropriate relations to rewrite the product of some generator
		
		INPUT:
		
		- "word", any list of positive integer entries between 1 and n
		- "start", together with "end" specifies a slice of the list "word." This slice should represent a braid diagram where the last pair of crossings are adjacent, the length is at least 3, and the two strands which appear in the last two crossings but DON'T cross there, DO cross somewhere else in the slice.
		- "end", together with "start" specifies a slice of the list "word" subject to the criteria above.
		- "p", a permutation corresponding to some idempotent e_p
		
		OUTPUT:
		
		- An element of this KLR algebra (writting in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k], but each term is closer to being reduced
		"""
		# Base case: If I have something like (i)(i+1)(i) at the end, then apply the braid relation.
		if word[end-1] == word[end-3]:
			new_p_m = self._P(self._g_ify(word[end:])(tuple(p)))
			if word[end-2] == word[end-1] + 1:
				k = word[end-1]
				new = self._reduce(word[:end-3]+[word[end-2], word[end-1], word[end-2]]+word[end:], p)
				c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
					new_g_m = self._Sn(())
					new_p_r = p
					new_p_l = new_p_m
					for r in range(-c_ij):
						new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
						new -= self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(word[:end-3], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[end:], new_p_r))
				return new
				
			if word[end-2] == word[end-1] - 1:
				k = word[end-2]
				new = self._reduce(word[:end-3]+[word[end-2], word[end-1], word[end-2]]+word[end:], p)
				c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
					new_g_m = self._Sn(())
					new_p_r = p
					new_p_l = new_p_m
					for r in range(-c_ij):
						new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
						new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * self._reduce(word[:end-3], new_p_l) * (self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[end:], new_p_r))
				return new
				
	
	
		# I am guaranteed to have those two distinguished strands cross SOMEWHERE, find that
		# location
		j = self._find_crossing(word[:end])
	
	
		# Another semi-base case: If that crossing can be moved upward in the diagram, then
		# do it.
		if abs(word[j] - word[j+1]) > 1:
			return self._straighten(word[:j]+[word[j+1], word[j]]+word[j+2:], start, end, p)
	
	
		# Now we are guaranteed to be in an upside-down situation from what we started with
		# BUT WITH SHORTER LENGTH. So straighten this shorter subword
		
		return self._straighten_reverse(word, j, end, p)
	
	
	# NEED TO FLIP EVERYTHING AROUND
	def _straighten_reverse(self, word, start, end, p):
		r"""
		Same as _straighten, but the algorithm is reversed. See documentation for _straighten.
		"""
		# Base case: If I have something like (i)(i+1)(i), then apply the braid relation.
		if word[start] == word[start+2]:
			new_p_m = self._P(self._g_ify(word[start+3:])(tuple(p)))
			if word[start+1] == word[start] + 1:
				k = word[start]
				new = self._reduce(word[:start]+[word[start+1], word[start], word[start+1]]+word[start+3:], p)
				c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
					new_g_m = self._Sn(())
					new_p_l = new_p_m
					new_p_r = p
					for r in range(-c_ij):
						new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
						new -= self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * self._reduce(word[:start], new_p_l) * (self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[end:], new_p_r))
				return new
				
			if word[start+1] == word[start] - 1:
				k = word[start+1]
				new = self._reduce(word[:start]+[word[start+1], word[start], word[start+1]]+word[start+3:], p)
				c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
					new_g = self._Sn(())
					new_p_r = p
					new_p_l = new_p_m
					for r in range(-c_ij):
						new_v = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
						new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * self._reduce(word[:start], new_p_l) * (self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[end:], new_p_r))
				return new
				
	
		# I am guaranteed to have those two distinguished strands cross SOMEWHERE, find that
		# location
		j = start+self._find_crossing_reverse(word[start:])
	
	
		# Another semi-base case: If that crossing can be moved downward in the diagram, then
		# do it.
		if abs(word[j] - word[j-1]) > 1:
			return self._straighten_reverse(word[:j-1]+[word[j], word[j-1]]+word[j+1:], start, end, p)
	
	
		# Now we are guaranteed to be in an upside-down situation from what we started with
		# BUT WITH SHORTER LENGTH. So straighten this shorter subword
		
		return self._straighten(word, start, j+1, p)

	'''
	Input should be the same as in the method "_straighten" above. Returns the index at which
	those two distinguished strands cross.
	'''
	def _find_crossing(self, word):
		r"""
		Returns the index at which a particular crossing occurs in the braid diagram associated with "word"
		
		INPUT:
		
		- "word", should represent a braid diagram where the last pair of crossings are adjacent, the length is at least 3, and the two strands which appear in the last two crossings but DON'T cross there, DO cross somewhere else.
		
		OUTPUT:
		
		- The index at which the two distinguished strands above cross
		"""
		# Make a word with an extra crossing at the top, in particular a crossing of the two
		# strands of interest
		temp = word+[word[-2]]
	
	
		# We'll be looking for the index at which those two strands cross again, which is the
		# same as the latest index at which we have an unreduced subword on the tail.
		j = len(word)
		while self._g_ify(temp[j:]).length() == len(temp[j:]):
			j -= 1
		return j
	
	
	
	'''
	Input should be the same as in the method "_reverse_straighten" above. Returns the index at which
	those two distinguished strands cross.
	'''
	def _find_crossing_reverse(self, word):
		r"""
		Same as _find_crossing, but the algorithm is reversed. See documentation for _find_crossing.
		"""
		# Make a word with an extra crossing at the bottom, in particular a crossing of the two
		# strands of interest
		temp = [word[1]]+word
	
	
		# We'll be looking for the index at which those two strands cross again, which is the
		# same as the earliest index at which we have an unreduced subword on the head.
		j = 1
		while self._g_ify(temp[:j+1]).length() == len(temp[:j+1]):
			j += 1
		return j
		

	def _TT(self, left_g, right_g, p):
		r"""
		Returns the product of the KLR algebra elements corresponding to "left_g" and "right_g" with 'matching' idempotents corresponding to "p"
		
		INPUT:
		
		- "left_g", an element of S_n
		- "right_g", an element of S_n
		- "p", a permutation corresponding to some idempotent e_p
		
		OUTPUT:
		
		- The product of t_{left_g}e_{right_g(p)} and t_{right_g}e_p written in the distinguished basis
		"""
		
		word = self.reduced_word(left_g)+self.reduced_word(right_g)
		return self._reduce(word, p)
	
	
	
	'''
	called when we need to swap t with x
	'''
	def _switch(self, word, v, p):
		r"""
		Returns the KLR algebra element corresponding to the product of appropriate generators determined by "word", "v" and "p"
		
		INPUT:
		
		- "word", a list of positive integers between 1 and n
		- "v", an integer vector in self._IVn
		- "p", a permutation corresponding to some idempotent e_p
		
		OUTPUT:
		
		- The product t_{s_1}t_{s_2}...t_{s_k}X^v e_p where word = [s_1, s_2, ..., s_k] written in the distinguished basis
		"""
		
		if len(word) == 0 or v == self._IVn((0,)*self._n):
			new_g = self._g_ify(word)
			new_v = v
			new_p = p
			return self.monomial(self._CP((new_v, new_g, new_p)))
		k = word[-1]
		i = self._get_first_nonzero_index(v)
		l = i + 1
		
		
		
		left_word = word[:-1]
		new_index = self._g_ify([k])(l)
		left_v = self._IVn((0,)*(new_index-1) + (1,) + (0,)*(self._n-new_index))
		left_p = self._P(self._g_ify([k])(tuple(p)))
		
		right_v = []
		for j in range(len(v)):
			if j == i:
				right_v += [v[j]-1]
			else:
				right_v += [v[j]]
		right_v = self._IVn(right_v)
		
		right_word = [word[-1]]
		right_p = p
		
		
		first = self._switch(left_word, left_v, left_p) * self._switch(right_word, right_v, right_p)
		
		
		if p[k-1] == p[k] and (l==k+1 or l==k):
			second_v = right_v
			second_word = word[:-1]
			if l == k+1:
				second = self._switch(second_word, second_v, p)
			elif l == k:
				second = -1*self._switch(second_word, second_v, p)
		else:
			second = self.zero()
		
		return first + second
		
		
	'''
	lazy function
	'''
	def _TX(self, left_g, right_v, p):
		r"""
		Returns the product of the KLR algebra elements corresponding to "left_g" and "right_v" with matching idempotents e_p
		
		INPUT:
		
		- "left_g", an element of S_n
		- "right_v", an integer vector of length n in self._IVn
		- "p", a permutation corresponding to some idempotent e_p
		
		OUTPUT:
		
		- The product t_{left_g} * X^v * e_p
		"""
		word = self.reduced_word(left_g)
		return self._switch(word, right_v, p)
	
		
	def product_on_basis(self, left, right):
		r"""
		Returns the product of two basis elements.
		
		INPUT:
		
		- "left", the basis element on the left
		- "right", the basis element on the right
		
		OUTPUT:
		
		- The product left * right written in the distinguished basis
		
		"""
		left_p = left.get_permutation()
		right_p = right.get_permutation()
		left_g = left.get_Sn_element()
		right_g = right.get_Sn_element()
		left_v = left.get_integer_vector()
		right_v = right.get_integer_vector()
		
		# If idempotents don't match
		if right_g.inverse()(tuple(left_p)) != tuple(right_p):
			return self.zero()
		
		# If left doesn't have any t's
		if left_g == self._Sn(()):
			new_v = self._add_integer_vectors(left_v, right_v)
			new_g = right_g
			new_p = right_p
			return self.monomial(self._CP((new_v, new_g, new_p)))
		
		# If right doesn't have any x's
		if right_v == self._IVn((0,)*self._n):
			new_v = left_v
			new_g = self._Sn(())
			new_p = self._P((left_g*right_g)(tuple(right_p)))
			return self.monomial(self._CP((new_v, new_g, new_p)))*self._TT(left_g, right_g, right_p)
		
		# Finally, assume left has some t's, right has some x's.
		new_v_l = left_v
		new_g_l = self._Sn(())
		new_p_l = self._P(left_g(tuple(left_p)))
		new_v_r = self._IVn((0,)*self._n)
		new_g_r = right_g
		new_p_r = right_p
		
		return self.monomial(self._CP((new_v_l, new_g_l, new_p_l)))*(self._TX(left_g, right_v, left_p)*self.monomial(self._CP((new_v_r, new_g_r, new_p_r))))
		
	def product(self, left, right):
		r"""
		Return the product of two elements of this KLR algebra
		
		INPUT:
		
		- "left", the element on the left
		- "right", the element on the right
		
		OUTPUT:
		
		- The product left * right
		
		"""
		sum = self.zero()
		for left_term in left.monomial_coefficients().keys():
			for right_term in right.monomial_coefficients().keys():
				sum += left[left_term]*(right[right_term]*self.product_on_basis(self.monomial(left_term), self.monomial(right_term)))
		return sum
	
	#WORKING ON THIS TOO
	def degree_on_basis(self, elt):
		r"""
		Returns the degree of the KLR algebra element associated to this basis element
		"""
		
		rw = self.reduced_word(elt[1])
		iv = elt[0]
		perm = elt[2]
		
		alpha = self.root_system().root_lattice().basis()
		
		sum = 0
		for k in rw[::-1]:
			sum += -alpha[perm[k-1]].symmetric_form(alpha[perm[k]])
			perm = self._P(self._Sn((k, k+1))(tuple(perm)))
		for k in range(1, len(iv)+1):
			sum += iv[k-1]*(alpha[perm[k-1]].symmetric_form(alpha[perm[k-1]]))
		
		return sum
	
	
	def _g_ify(self, word):
		r"""
		Returns the element of S_n corresponding to word
		"""
		return self._Sn.prod([self._Sn((i, i+1)) for i in word])
	
	
	def _add_integer_vectors(self, first, second):
		r"""
		Returns the point-wise sum of the two integer vectors
		"""
		return self._IVn((first[i]+second[i] for i in range(len(first))))
	
	def _get_first_nonzero_index(self, v):
		r"""
		Returns the first index at which the integer vector "v" is non-zero
		"""
		for i in range(len(v)):
			if v[i] != 0:
				return i
		return -1
	
# '''
# Takes a list of positive integers and turns it into an S_n element by multiplying the
# appropriate transpositions.
# '''
# def g_ify(word, Sn=None):
# 	
# 	if Sn == None:
# 		# Make a symmetric group of the appropriate size.
# 		n = 0
# 		if word != []:
# 			n = max(word)
# 		Sn = SymmetricGroup(n+1)
# 	
# 	# Multiply the appropriate transpositions and return the result.
# 	return Sn.prod([Sn((i, i+1)) for i in word])
# 
# def add_integer_vectors(first, second):
# 	if first.parent() != second.parent():
# 		raise TypeError("%s has parent %s but %s has parent %s"%(first, first.parent(), second, second.parent()))
# 	if len(first) != len(second):
# 		raise ValueError("%s has length %s but %s has length %s"%(first, len(first), second, len(second)))
# 	P = first.parent()
# 	return P((first[i]+second[i] for i in range(len(first))))
# 
# def get_first_nonzero_index(v):
# 	for i in range(len(v)):
# 		if v[i] != 0:
# 			return i
# 	return -1


'''
# FOR TESTING:
RS = RootSystem(["A", 6])
RL = RS.root_lattice()
alpha = RL.basis()
a = alpha[1]+alpha[2]
K = KLRAlgebra(ZZ, RS, a)
x,t,e = K.x,K.t,K.e

'''


# GOOD FOR NIL-HECKE STUFF:

RS = RootSystem(["A",Infinity])
RL = RS.root_lattice()
alpha = RL.basis()
m = 5
a = m*alpha[1]
omega = [(i, m-i+1) for i in range(1, m/2+1)]
K = KLRAlgebra(ZZ, RS, a)
x, t, e = (K.x, K.t, K.e)
T = t(omega)
X = K.prod([x(i)**(i-1) for i in range(2, m+1)])

'''
# OTHER TESTING:

#a = alpha[1]+2*alpha[3]
#a = alpha[1]+2*alpha[3]+alpha[2]+alpha[5]
a = alpha[1]+2*alpha[2]+alpha[3]+alpha[4]
K = KLRAlgebra(ZZ, RS, a)
x, t, e = (K.x, K.t, K.e)
z = K.an_element()

elt1=x(2)*(x(3)*(x(3)*(x(4)*(x(4)*(x(4)*(x(5)*(x(5)*(x(5)*(x(5)*e((2,1,2,3,4)))))))))))
elt2=x(2)*(x(3)*(x(3)*(x(4)*(x(4)*(x(4)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*e((2,1,2,3,4))))))))))))
elt3=x(2)*(x(3)*(x(3)*(x(4)*(x(4)*(x(4)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*e((2,1,2,3,4)))))))))))))
elt4=x(2)*(x(3)*(x(3)*(x(4)*(x(4)*(x(4)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*(x(5)*e((2,1,2,3,4))))))))))))))
elt=elt1+elt2-elt3+elt4
#a = K.monomial(K._CP((K._IVn((0, 0, 0)), K._Sn((1, 2)), K._P((1, 3, 3)))))
#b = K.monomial(K._CP((K._IVn((1, 0, 0)), K._Sn(()), K._P((1, 3, 3)))))
#c = K.monomial(K._CP((K._IVn((0, 1, 0)), K._Sn((1, 2)), K._P((1, 3, 3)))))
'''