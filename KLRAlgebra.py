import sage
from sage.rings.ring import Algebra
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
#from sage.structure.element import AlgebraElement
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.bindable_class import BindableClass
from sage.combinat.integer_vector import IntegerVectors_k
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_algebras import GradedAlgebras
from sage.misc.cachefunc import cached_method

#I am putting a superfluous comment here

class BasisElement(CombinatorialFreeModule.Element):
	r"""
	Generic basis element. Only things we do here are define some printing methods.
	"""
	def __init__(self, *args, **kwds):
		r"""
		Initialize self.
		"""
		CombinatorialFreeModule.Element.__init__(self, *args, **kwds)

		P = self.parent().realization_of()
		self._n = P.positive_root().height()
		self.rename(self._get_print_name())

	def _get_print_name(self):
		r"""
		Return a string which represents self.

		This cannot be implemented at this level of generality.
		"""
		raise NotImplementedError('This should be implemented by subclasses')

	def tikz(self, x_sep=0.5, y_sep=0.3):
		r"""
		Return tikz code to which draws a picture of this element.

		INPUT:

		- `x_sep` -- float; the distance between strands

		- `y_sep` -- float; the vertical distance between crossings

		OUTPUT: The tikz code which generates a picture of this element

		NOTE:

			Requires the following in your preamble:
			`\usepackage{tikz}`
			`\usetikzlibrary{positioning}`

		Asks subclass for tikzmonomials (basis elements), receives them and lines them up with '+' in between them
		"""
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

class XTElement(BasisElement):
	'''
	An element of the basis of a KLR algebra which has dots at the top of diagrams
	and crossings at the bottom
	'''

	#Correct order
	def _get_print_name(self):
		'''

		'''
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

	#TODO: Fix coefficient height
	def _tikz_monomial(self, x_sep=0.5, y_sep=0.3, coeff=1, first_term=True, coeff_height=None):
		rw = self._get_reduced_word()[::-1]
		iv = self.get_integer_vector()
		perm = self.get_permutation()
		n = self._n
		l = len(rw)
		max_dots = max(iv)
		if max_dots == 0:
			dot_space = 0
		else:
			dot_space = ceil((max_dots+2)/3.0)*y_sep
		if coeff_height == None:
			coeff_height = (y_sep*(l+2)+dot_space)/2
		front_matter = "\\begin{tikzpicture}\n\t\\foreach \\x in {1, ..., %s}\n\t\t\\node (\\x_-1) at (%s*\\x, %s){};\n\t\\foreach \\x in {1, ..., %s}\n\t\t\\foreach \\y in {0, ..., %s}\n\t\t\t\\node (\\x_\\y) at (%s*\\x, %s*\\y){};\n"%(n, x_sep, y_sep*(l+2)+dot_space, n, l+2, x_sep, y_sep)
		if coeff == 1:
			coeff = ""
		else:
			coeff = "\\,\\,\\left(%s\\right)\\,\\,"%latex(coeff)


		if first_term:
			end_matter = "\t\\node[on grid, below=%s of 1_-1] (temp) {};\n\t\\node[left=0 of temp]{$%s$};\n\\end{tikzpicture}"%(coeff_height, coeff)
		else:
			end_matter = "\t\\node[on grid, below=%s of 1_-1] (temp) {};\n\t\\node[left=0 of temp]{$+%s$};\n\\end{tikzpicture}"%(coeff_height, coeff)



		simple_roots = self.parent().realization_of().root_system().root_lattice().simple_roots()


		color_list = rainbow(len(set(perm)))
		colors = {list(set(perm))[i]:color_list[i] for i in range(len(list(set(perm))))}


		color_string = ""
		for i in colors:
			color_string += "\t\\definecolor{color%s}{HTML}{%s}\n"%(i, colors[i][1:].upper())

		strands = {i:"\t\\draw [color=color%s] (%s_0.center)--(%s_1.center)"%(perm[i-1], i, i) for i in range(1,n+1)}
		for i in range(l):
			strands[rw[i]], strands[rw[i]+1] = (strands[rw[i]+1]+"--(%s_%s.center)"%(rw[i], i+2), strands[rw[i]]+"--(%s_%s.center)"%(rw[i]+1, i+2))
			for j in range(1, n+1):
				if j != rw[i] and j != rw[i]+1:
					strands[j] += "--(%s_%s.center)"%(j, i+2)
		for j in range(1, n+1):
			strands[j] += "--(%s_%s.center)--(%s_-1.center);\n"%(j, l+2, j)

		dots = {i:"\t\\foreach \\y in {1, ..., %s}\n\t\t\\node[on grid, circle, fill, inner sep=1pt, below=%s/%s*\\y of %s_-1]{};\n"%(iv[i-1], dot_space, iv[i-1]+1, i) for i in range(1, n+1)}


		the_string = front_matter+color_string
		for i in strands:
			the_string += strands[i]
		for i in dots:
			if iv[i-1] != 0:
				the_string += dots[i]
		the_string += end_matter
		return the_string

class TXElement(BasisElement):
	'''
	An element of the basis of a KLR algebra which has crossings at the top of
	diagrams and dots at the bottom
	'''

	#Correct Order
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
					new_name += "%s%s%s + "%(T,X,E)
				else:
					new_name += "%s*%s%s%s + "%(coeff,T,X,E)
			new_name = new_name[:-3]
			return new_name
		else:
			return "0"

	#Correct Order
	def _tikz_monomial(self, x_sep=0.5, y_sep=0.3, coeff=1, first_term=True, coeff_height=None):
		rw = self._get_reduced_word()[::-1]
		iv = self.get_integer_vector()
		perm = self.get_permutation()
		n = self._n
		l = len(rw)
		max_dots = max(iv)
		if max_dots == 0:
			dot_space = 0
		else:
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



		simple_roots = self.parent().realization_of().root_system().root_lattice().simple_roots()


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

class KLRAlgebra(UniqueRepresentation, Parent):
	def __init__(self, base_ring, positive_root, signs=None, *args, **kwargs):
		'''
		base_ring is the ring over which to define this KLR algebra

		positive_root is an element of the positive root lattice of some root
		system. It knows which root system it came from.

		signs is a function which takes as input two integers i and j, and
		outputs either plus 1 or minus 1. Real people almost never need to
		specify a non-default function.

		there's a mysterious _custom_signs instance variable
		'''
		category = GradedAlgebras(base_ring)
		self._base = base_ring
		Parent.__init__(self, category=category.WithRealizations())

		self._R = base_ring
		self._root_system = positive_root.parent().root_system

		# Get the cartan matrix. If we are not of finite type, this won't work so I have
		# to make a lazy family. The only supported infinite type is ["A", ZZ] right now.
		# To add support for other infinite types, one needs only to define the proper
		# function giving the appropriate entry in the infinite-dimensional cartan matrix
		# for that infinite-type. Should also implement "symmetric_form" below......
		try:
			self._cartan_matrix = self._root_system.cartan_matrix()
		except AttributeError as e:
			if self._root_system.cartan_type() == CartanType(["A", ZZ]):
				self._cartan_matrix = Family(ZZ, lambda i: Family(ZZ, lambda j: 2*kronecker_delta(i, j) - kronecker_delta(i-1, j) - kronecker_delta(i, j-1)))
			else:
				raise e

		self._positive_root = positive_root

		# I forget where I need to use the _custom_signs instance variable
		if signs == None:
			self._custom_signs = False
			self._signs = lambda i, j: sign(j-i) #sign is a built-in sage function
		else:
			self._custom_signs = True
			self._signs = signs

		XT = KLRAlgebra.XT(self)
		TX = KLRAlgebra.TX(self)

		my_category = self.GradedBases()

		TX_to_XT = TX.module_morphism(XT._TX_to_XT_on_basis, category=my_category, codomain=XT)
		XT_to_TX = XT.module_morphism(TX._XT_to_TX_on_basis, category=my_category, codomain=TX)

		TX_to_XT.register_as_coercion()
		XT_to_TX.register_as_coercion()


		# PLAYING AROUND WITH PRINTING STUFF.........

		# defaults = ('dot_symbol':'x', 'crossing_symbol':'t', 'strand_symbol':'e')
		# for arg in kwargs:
		# 	if arg not in defaults:
		# 		raise ValueError('%s is not a valid print option'%arg)
		# self._print_options = {key:kwargs.get(key, defaults[key]) for key in defaults}
		#
		# for realization in self.realizations():

	def _repr_(self):
		r"""
		Return the string representation of self.
		"""
		name = "KLR Algebra over %s for %s with positive root %s"%(repr(self._R), repr(self._root_system), repr(self._positive_root))
		if self._custom_signs:
			name += " with custom sign conventions"
		return name

	def change_ring(self, ring=None):
		r"""
		Return the KLR algebra identical to this one, over the specified ring
		"""
		if ring == None:
			ring = self.base_ring()
		if self._custom_signs:
			signs = self._signs
		else:
			signs = None
		return KLRAlgebra(ring, self._positive_root, signs=signs)

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

	def a_realization(self):
		return self.XT()

	def algebra_generators(self):
		return self.a_realization().algebra_generators()

	monoid_generators = algebra_generators

	class GradedBases(Category_realization_of_parent):
		def super_categories(self):
			A = self.base()
			category = GradedAlgebrasWithBasis(A.base_ring())
			return [A.Realizations(), category.Realizations()]

		class ElementMethods:
			def _get_reduced_word(self):
				my_parent = self.parent()
				if not self.is_monomial():
					raise ValueError("%s is not a monomial"%self)
				return my_parent.reduced_word(self.get_Sn_element())

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

			def tikz(self, x_sep=0.5, y_sep=0.3):
				the_string = ""
				mcs = self.monomial_coefficients()
				min_height = Infinity
				for monomial in mcs:
					my_mon = self.parent().monomial(monomial)
					rw = self.parent().reduced_word(monomial[1])
					iv = monomial[0]
					l = len(rw)
					max_dots = max(iv)
					dot_space = ceil((max_dots+2)/3.0)*y_sep
					this_height = (y_sep*(l+2)+dot_space)/2
					if this_height < min_height:
						min_height = this_height
				the_terms = []
				first_term=True
				for monomial in mcs:
					the_terms += ["\n"+self.parent().monomial(monomial)._tikz_monomial(x_sep=x_sep, y_sep=y_sep, coeff=mcs[monomial], first_term=first_term, coeff_height=min_height)]
					first_term=False
				for term in the_terms:
					the_string += term
				return the_string[1:]

			def plot(self, coefficient_padding=1, **kwargs):
				cp = coefficient_padding
				par = self.parent()
				mcs = self.monomial_coefficients()
				width = par._n
				heights = {mc:par.plot_on_basis(mc, height_only=True) for mc in mcs}
				pieces = []
				h_shift = 0
				first_term = True
				for mc in mcs:
					if first_term:
						if mcs[mc] == par.base_ring().one():
							h_shift += 0
						else:
							coeff = text("$\\left(%s\\right)$"%latex(mcs[mc]), (h_shift, 0), color='black', **kwargs)
							pieces += [coeff]
							h_shift += cp
					else:
						if mcs[mc] == par.base_ring().one():
							coeff = text("$+$", (h_shift, 0), color='black', **kwargs)
						else:
							coeff = text("$+\\left(%s\\right)$"%latex(mcs[mc]), (h_shift, 0), color='black', **kwargs)
						pieces += [coeff]
						h_shift += cp
					first_term = False
					pieces += [par.plot_on_basis(mc, h_shift=h_shift, v_shift=-heights[mc]/2.0, **kwargs)]
					h_shift += width + cp - 1
				G = sum(pieces)
				G.axes(False)
				G.set_aspect_ratio(1)
				return G

		class ParentMethods:
			def _name(self):
				if self._custom_words:
					if self._A._custom_signs:
						return str(self._A) + " and custom reduced words"
					else:
						return str(self._A) + " with custom reduced words"
				else:
					return str(self._A)

			def _repr_(self):
				return "%s in %s basis"%(self._name(), self._realization_name())

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

			def x(self, i=None):
				r"""
				Return the monomial associated to i.

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

			def longest_crossing(self):
				return self.t(self.symmetric_group().long_element())

			# Do I have the best way of inputting parameters? It seems clunky
			# and unintuitive... Should I handle errors?
			def divided_power_idempotent(self, *args):
				the_is, the_ms = [arg[0] for arg in args], [arg[1] for arg in args]
				my_groups = []
				v = []
				start = 1
				for m in the_ms:
					end = start + m
					my_groups += [SymmetricGroup(range(start, end))]
					v += range(end-start)
					start = end

				g = self._Sn.prod([self._Sn(G.long_element()) for G in my_groups])
				p = self._P(sum([(the_is[i],)*the_ms[i] for i in range(len(the_is))], ()))
				v = self._IVn(tuple(v))

				return self.monomial(self._CP((v, g, p)))

			def one(self):
				r"""
				Returns the identity element of this KLR algebra. This is the sum of all idempotents corresponding to length n permutations whose entries correspond to the coefficients of the positive root associated with this KLR algebra.
				"""
				g = self._Sn.one()
				v = self._IVn((0,)*self._n)
				return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])

			def algebra_generators(self):
				r"""
				Returns the algebra generators. See x(), t(), e().
				"""
				return [self.x(i) for i in range(1, self._n+1)]+[self.t(i) for i in range(1, self._n)]+[self.e(p) for p in self._P]

			def cartan_matrix(self):
				return self.realization_of().cartan_matrix()

			def signs(self, pair=None):
				return self.realization_of().signs(pair=pair)

			def symmetric_group(self):
				r"""
				Return the symmetric group associated to this KLR Algebra. This is the symmetric group S_n where n is the height of the positive root associated to this KLR Algebra.
				"""
				return self._Sn

			def symmetric_form(self, left, right):
				try:
					return left.symmetric_form(right)
				except AttributeError as e:
					if self._root_system.cartan_type() == CartanType(["A", ZZ]):
						mcl = left.monomial_coefficients()
						mcr = right.monomial_coefficients()
						return sum([mcl[i]*mcr[j]*self.cartan_matrix()[i-1][j-1] for i in mcl for j in mcr])
					else:
						raise e

			def degree_on_basis(self, elt):
				r"""
				Returns the degree of the KLR algebra element associated to this basis element
				"""
				rw = self.reduced_word(elt[1])
				iv = elt[0]
				perm = elt[2]

				alpha = self._A.root_system().root_lattice().basis()

				sum = 0
				for k in rw[::-1]:

					sum -= self.symmetric_form(alpha[perm[k-1]],alpha[perm[k]])
					perm = self._P(self._Sn((k, k+1))(tuple(perm)))
				for k in range(1, len(iv)+1):
					sum += iv[k-1]*(self.symmetric_form(alpha[perm[k-1]],alpha[perm[k-1]]))

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
				Returns the first index at which the integer vector "v" is non-zero, -1 if never.
				"""
				for i in range(len(v)):
					if v[i] != 0:
						return i
				return -1

			def _get_last_nonzero_index(self, v):
				r"""
				Returns the last index at which the integer vector "v" is non-zero, -1 if never.
				"""
				for i in range(len(v))[::-1]:
					if v[i] != 0:
						return i
				return -1

			@cached_method
			def _get_move(self, word, target):
				if word == target:
					return ("n", 0)


				assert len(word) == len(target)

				# Check to see how far we match the target
				i = 0
				try:
					while word[i] == target[i]:
						i += 1
				# Error is thrown if nothing matches. We don't care.
				except IndexError:
					pass

				# If anything matches, just check the leftovers.
				if i > 0:
					my_move = list(self._get_move(word[i:], target[i:]))
					my_move[1] += i
					return tuple(my_move)

				# Now the beginning doesn't match. Find where the target first crossing
				# crosses in word.
				temp = (target[0],) + word
				j = 0
				while self._g_ify(temp[:j]).length() == j:
					j += 1
				j -= 2
				# j is now the index at which the distinguished strands cross in word.
				return self._float_move(word[:j+1])

			'''
			give a move which floats the bottom crossing to the top
			'''
			@cached_method
			def _float_move(self, word):
				# if the two bottom crossings are far apart, just move them past
				# each other for free. Do this until you can't anymore or until the crossing
				# is all the way at the top.
				j = len(word)-1
				if j > 0 and abs(word[j] - word[j-1]) > 1:
					return ('s', j-1)

				# Now we can assume that the crossing of interest is as far up as possible.
				# If it's at the top, then we are done with this part, just go back to
				# reduce_full
				# if j == 0:
				# 	return word

				# Now we can assume it got stuck on some adjacent transposition.
				#print word, j
				i = self._find_crossing(word)
				if i == j-2:
					return ('b', j-2)
				my_move = list(self._sink_move(word[i:]))
				my_move[1] += i
				return tuple(my_move)

			'''
			give a move which floats the top crossing to the bottom
			'''
			@cached_method
			def _sink_move(self, word):
				my_move = list(self._float_move(word[::-1]))
				if my_move[0] == "s":
					my_move[1] = len(word) - 2 - my_move[1]
				elif my_move[0] == "b":
					my_move[1] = len(word) - 3 - my_move[1]
				return tuple(my_move)

			# ALMOST READY TO REWRITE THIS ACCORDING TO Reduce.py SO THE
			# COMPUTATIONS WILL ACTUALLY TERMINATE!
			@cached_method
			def _shuffle(self, word, p):
				r"""
				Returns an element of this KLR algebra which is equal to the product of the generators corresponding to the reduced word "word" and the idempotent determined by p

				INPUT:

				- "word", a list representing a reduced word for some S_n element
				- "p", a permutation corresponding to some idempotent e_p
				- "Graph", the reduced word graph containing the reduced word "word"

				OUTPUT:

				- An element of this KLR algebra (written in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k]
				"""

				g = self._g_ify(word)
				good_word = tuple(self.reduced_word(g))

				if word == good_word:
					v = self._IVn((0,)*self._n)
					return self.monomial(self._CP((v, g, p)))

				my_move = self._get_move(word, good_word)
				word = list(word)
				while my_move[0] == 's':
					i = my_move[1]
					word[i], word[i+1] = word[i+1], word[i]
					my_move = self._get_move(tuple(word), good_word)
				word = tuple(word)

				if my_move[0] == 'n':
					v = self._IVn((0,)*self._n)
					return self.monomial(self._CP((v, g, p)))



				# Now there IS a next word, and the current word and next word differ by a braid relation.
				# Find the index at which this braid relation occurs. i is the index along SP at which
				# current_word lies.
				j = my_move[1]
				next_word = word[:j] + (word[j+1], word[j], word[j+1]) + word[j+3:]

				# Now j is the first index at which current_word and next_word differ


				# end --> j+3, word --> current_word

				#OLLLLLLDDDDD
				# new_p_m = self._P(self._g_ify(current_word[j+3:])(tuple(p)))
				# if current_word[j+1] == current_word[j+2] + 1:
				# 	k = current_word[j+2]
				# 	new = self._shuffle(current_word[:j]+(current_word[j+1], current_word[j+2], current_word[j+1])+current_word[j+3:], p)
				# 	c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				# 	if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
				# 		new_g_m = self._Sn(())
				# 		new_p_r = p
				# 		new_p_l = new_p_m
				# 		for r in range(-c_ij):
				# 			new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
				# 			new -= self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(current_word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(current_word[j+3:], new_p_r))
				# 	return new
				#
				# if current_word[j+1] == current_word[j+2] - 1:
				# 	k = current_word[j+1]
				# 	new = self._shuffle(current_word[:j]+(current_word[j+1], current_word[j+2], current_word[j+1])+current_word[j+3:], p)
				# 	c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
				# 	if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
				# 		new_g_m = self._Sn(())
				# 		new_p_r = p
				# 		new_p_l = new_p_m
				# 		for r in range(-c_ij):
				# 			new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
				# 			new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(current_word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(current_word[j+3:], new_p_r))
				# 	return new

				# NEWWWWWWWWWW
				new_p_m = self._P(self._g_ify(word[j+3:])(tuple(p)))
				if word[j+1] == word[j+2] + 1:
					k = word[j+2]
					new = self._shuffle(word[:j]+(word[j+1], word[j+2], word[j+1])+word[j+3:], p)
					c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
					if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
						new_g_m = self._Sn(())
						new_p_r = p
						new_p_l = new_p_m
						for r in range(-c_ij):
							new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
							new -= self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[j+3:], new_p_r))
					return new

				if word[j+1] == word[j+2] - 1:
					k = word[j+1]
					new = self._shuffle(word[:j]+(word[j+1], word[j+2], word[j+1])+word[j+3:], p)
					c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
					if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
						new_g_m = self._Sn(())
						new_p_r = p
						new_p_l = new_p_m
						for r in range(-c_ij):
							new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
							new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * (self._reduce(word[:j], new_p_l) * self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[j+3:], new_p_r))
					return new

			@cached_method
			def _reduce(self, word, p):
				r"""
				Returns an element of this KLR algebra which is equal to the product of the generators corresponding to "word" and the idempotent determined by "p"

				INPUT:

				- "word", any list containing positive integer entries between 1 and n
				- "p", a permutation corresponding to some idempotent e_p

				OUTPUT:

				- An element of this KLR algebra (written in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k]
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
					return self._reduce(word[:i-1]+(word[i], word[i-1])+word[i+1:], p)



				# Now I have guaranteed that my minimal unreduced word ends in two adjacent
				# transpositions, and a consequence of minimality is that the two strands which appear
				# in the last two transpositions but don't cross there, DO cross somewhere else in the
				# minimal unreduced word, before position 0. Straighten this piece.

				return self._straighten(word, j+1, i+1, p)

			@cached_method
			def _straighten(self, word, start, end, p):
				r"""
				Applies appropriate relations to rewrite the product of some generators

				INPUT:

				- "word", any list of positive integer entries between 1 and n
				- "start", together with "end" specifies a slice of the list "word." This slice should represent a braid diagram where the last pair of crossings are adjacent, the length is at least 3, and the two strands which appear in the last two crossings but DON'T cross there, DO cross somewhere else in the slice.
				- "end", together with "start" specifies a slice of the list "word" subject to the criteria above.
				- "p", a permutation corresponding to some idempotent e_p

				OUTPUT:

				- An element of this KLR algebra (written in the distinguished basis) which is equal to t_{s_1}t_{s_2} ... t_{s_k}e_p where word = [s_1, ..., s_k], but each term is closer to being reduced
				"""
				# Base case: If I have something like (i)(i+1)(i) at the end, then apply the braid relation.
				if word[end-1] == word[end-3]:
					new_p_m = self._P(self._g_ify(word[end:])(tuple(p)))
					if word[end-2] == word[end-1] + 1:
						k = word[end-1]
						new = self._reduce(word[:end-3]+(word[end-2], word[end-1], word[end-2])+word[end:], p)
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
						new = self._reduce(word[:end-3]+(word[end-2], word[end-1], word[end-2])+word[end:], p)
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
					return self._straighten(word[:j]+(word[j+1], word[j])+word[j+2:], start, end, p)


				# Now we are guaranteed to be in an upside-down situation from what we started with
				# BUT WITH SHORTER LENGTH. So straighten this shorter subword

				return self._straighten_reverse(word, j, end, p)

			# NEED TO FLIP EVERYTHING AROUND
			@cached_method
			def _straighten_reverse(self, word, start, end, p):
				r"""
				Same as _straighten, but the algorithm is reversed. See documentation for _straighten.
				"""
				# Base case: If I have something like (i)(i+1)(i), then apply the braid relation.
				if word[start] == word[start+2]:
					new_p_m = self._P(self._g_ify(word[start+3:])(tuple(p)))
					if word[start+1] == word[start] + 1:
						k = word[start]
						new = self._reduce(word[:start]+(word[start+1], word[start], word[start+1])+word[start+3:], p)
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
						new = self._reduce(word[:start]+(word[start+1], word[start], word[start+1])+word[start+3:], p)
						c_ij = self.cartan_matrix()[new_p_m[k-1]-1][new_p_m[k]-1]
						if c_ij < 0 and new_p_m[k-1] == new_p_m[k+1]:
							new_g_m = self._Sn(())
							new_p_r = p
							new_p_l = new_p_m
							for r in range(-c_ij):
								new_v_m = self._IVn((0,)*(k-1) + (r,) + (0,) + (-1-c_ij-r,) + (0,)*(self._n-k-2))
								new += self.signs(pair=(new_p_m[k-1]-1, new_p_m[k]-1)) * self._reduce(word[:start], new_p_l) * (self.monomial(self._CP((new_v_m, new_g_m, new_p_m))) * self._reduce(word[end:], new_p_r))
						return new


				# I am guaranteed to have those two distinguished strands cross SOMEWHERE, find that
				# location
				j = start+self._find_crossing_reverse(word[start:])


				# Another semi-base case: If that crossing can be moved downward in the diagram, then
				# do it.
				if abs(word[j] - word[j-1]) > 1:
					return self._straighten_reverse(word[:j-1]+(word[j], word[j-1])+word[j+1:], start, end, p)


				# Now we are guaranteed to be in an upside-down situation from what we started with
				# BUT WITH SHORTER LENGTH. So straighten this shorter subword

				return self._straighten(word, start, j+1, p)

			'''
			Input should be the same as in the method "_straighten" above. Returns the index at which
			those two distinguished strands cross.
			'''
			@cached_method
			def _find_crossing(self, word):
				r"""
				Returns the index at which a particular crossing occurs in the braid diagram associated with "word"

				INPUT:

				- "word", should represent a braid diagram where the last pair of crossings are adjacent, the length is at least 3, and the two strands which appear in the last two crossings but DON'T cross there, DO cross somewhere else.

				OUTPUT:

				- The index at which the two distinguished strands above cross
				"""
				# Make a word with an extra crossing at the bottom, in particular a crossing of the two
				# strands of interest
				temp = word+(word[-2],)


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
			@cached_method
			def _find_crossing_reverse(self, word):
				r"""
				Same as _find_crossing, but the algorithm is reversed. See documentation for _find_crossing.
				"""
				# Make a word with an extra crossing at the bottom, in particular a crossing of the two
				# strands of interest
				temp = (word[1],)+word


				# We'll be looking for the index at which those two strands cross again, which is the
				# same as the earliest index at which we have an unreduced subword on the head.
				j = 1
				while self._g_ify(temp[:j+1]).length() == len(temp[:j+1]):
					j += 1
				return j

			@cached_method
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

				word = tuple(self.reduced_word(left_g))+tuple(self.reduced_word(right_g))
				return self._reduce(word, p)

	class Basis(CombinatorialFreeModule, BindableClass):
		def __init__(self, A, reduced_words=None):
			r"""
			See "KLRAlgebra" for full documentation.
			"""

			self._A = A

			self._n = self._A.positive_root().height()
			self._Sn = SymmetricGroup(self._n)

			self._IVn = IntegerVectors_k(self._n)
			IndexList = []
			mcs = self._A.positive_root().monomial_coefficients()
			for key in mcs:
				IndexList += mcs[key]*[key]
			self._P = Permutations(IndexList)
			self._CP = cartesian_product([self._IVn, self._Sn, self._P])

			if reduced_words == None:
				self._custom_words = False
				self._reduced_words = lambda g: self._Sn(g).reduced_word()
			else:
				self._custom_words = True
				self._reduced_words = reduced_words

			CombinatorialFreeModule.__init__(self, A._R, self._CP, category=A.GradedBases())

			if self._custom_words:
				my_category = self._A.GradedBases()
				if self._realization_name() == 'XT':
					XT = KLRAlgebra.XT(self._A)
					XTR = self
					XT_to_XTR = XT.module_morphism(XTR._XT_to_XTR_on_basis, category=my_category, codomain=XTR)
					XTR_to_XT = XTR.module_morphism(XTR._XTR_to_XT_on_basis, category=my_category, codomain=XT)
					XT_to_XTR.register_as_coercion()
					XTR_to_XT.register_as_coercion()

				elif self._realization_name() == 'TX':
					TX = KLRAlgebra.TX(self._A)
					TXR = self
					TX_to_TXR = TX.module_morphism(TXR._TX_to_TXR_on_basis, category=my_category, codomain=TXR)
					TXR_to_TX = TXR.module_morphism(TXR._TXR_to_TX_on_basis, category=my_category, codomain=TX)
					TX_to_TXR.register_as_coercion()
					TXR_to_TX.register_as_coercion()

	class XT(Basis):
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

		Element = XTElement

		def _TX_to_XT_on_basis(self, elt):
			word = tuple(self.realization_of().TX().reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			return self._switch(word, v, p)

		def _XT_to_XTR_on_basis(self, elt):
			XT = self.realization_of().XT()
			word = tuple(XT.reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			left_p = self._P(elt[1](tuple(p)))
			return self.monomial(self._CP((v, self._Sn(()), left_p))) * self._shuffle(word, p)

		def _XTR_to_XT_on_basis(self, elt):
			XT = self.realization_of().XT()
			word = tuple(self.reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			left_p = XT._P(elt[1](tuple(p)))
			return XT.monomial(XT._CP((v, XT._Sn(()), left_p))) * XT._shuffle(word, p)


		'''
		called when we need to swap t with x
		'''
		@cached_method
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

			right_word = (word[-1],)
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
		@cached_method
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
			word = tuple(self.reduced_word(left_g))
			return self._switch(word, right_v, p)

		@cached_method
		def product_on_basis(self, left, right):
			r"""
			Returns the product of two basis elements.

			INPUT:

			- "left", the basis element on the left
			- "right", the basis element on the right

			OUTPUT:

			- The product left * right written in the distinguished basis

			"""
			left_p = left[2]
			right_p = right[2]
			left_g = left[1]
			right_g = right[1]
			left_v = left[0]
			right_v = right[0]

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

		def plot_on_basis(self, elt, v_shift=None, h_shift=None, height_only=False, **kwargs):
			if v_shift == None:
				v = 0
			else:
				v = v_shift
			if h_shift == None:
				h = 0
			else:
				h = h_shift

			rw = self.reduced_word(elt[1])[::-1]
			iv = elt[0]
			perm = elt[2]
			n = self._n
			l = len(rw)

			max_dots = max(iv)
			# if max_dots == 0:
			# 	dot_space = 0
			# else:
			# 	dot_space = ceil((max_dots+2)/3.0)
			# #temp
			if rw == []:
				dot_space = max(1, max_dots/5.0)
			else:
				dot_space = max_dots/5.0

			if height_only:
				return dot_space+l

			color_list = rainbow(len(set(perm)))
			colors = {list(set(perm))[i]:color_list[i] for i in range(len(list(set(perm))))}

			strands = {i:[(i-1+h, v)] for i in range(1,n+1)}
			for i in range(l):
				strands[rw[i]], strands[rw[i]+1] = strands[rw[i]+1]+[(rw[i]-1+h, i+1+v)], strands[rw[i]]+[(rw[i]+h, i+1+v)]
				for j in range(1, n+1):
					if j != rw[i] and j != rw[i]+1:
						strands[j] += [(j-1+h, i+1+v)]
			for j in range(1, n+1):
				strands[j] += [(j-1+h, l+v+dot_space)]


			# strands = {i:"\t\\draw [color=color%s] (%s_-1.center)--(%s_0.center)--(%s_1.center)"%(perm[i-1], i, i, i) for i in range(1,n+1)}
			# for i in range(l):
			# 	strands[rw[i]], strands[rw[i]+1] = (strands[rw[i]+1]+"--(%s_%s.center)"%(rw[i], i+2), strands[rw[i]]+"--(%s_%s.center)"%(rw[i]+1, i+2))
			# 	for j in range(1, n+1):
			# 		if j != rw[i] and j != rw[i]+1:
			# 			strands[j] += "--(%s_%s.center)"%(j, i+2)
			# for j in range(1, n+1):
			# 	strands[j] += "--(%s_%s.center);\n"%(j, l+2)

			the_lines = [line(strands[i], color=colors[perm[int(strands[i][0][0]-h)]], **kwargs) for i in range(1, n+1)]


			dots = []
			for i in range(n):
				n_dots = iv[i]
				for j in range(n_dots):
					dots += [circle((i+h, -(j+1)*dot_space/(n_dots+1.0) + v + l + dot_space), 0.05, fill=True, color='black', **kwargs)]



			G = sum(the_lines)+sum(dots)
			G.set_aspect_ratio(1)
			G.axes(False)
			return G

	class TX(Basis):
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

		Element = TXElement

		def _XT_to_TX_on_basis(self, elt):
			word = tuple(self.realization_of().XT().reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			return self._switch(word, v, p)

		def _TX_to_TXR_on_basis(self, elt):
			TX = self.realization_of().TX()
			word = tuple(TX.reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			return self._shuffle(word, p) * self.monomial(self._CP((v, self._Sn(()), p)))

		def _TXR_to_TX_on_basis(self, elt):
			TX = self.realization_of().TX()
			word = tuple(self.reduced_word(elt[1]))
			v = elt[0]
			p = elt[2]
			return TX._shuffle(word, p) * TX.monomial(TX._CP((v, TX._Sn(()), p)))

		'''
		called when we need to swap t with x
		'''
		@cached_method
		def _switch(self, word, v, p):
			r"""
			Returns the KLR algebra element corresponding to the product of appropriate generators determined by "word", "v" and "p"

			INPUT:

			- "word", a list of positive integers between 1 and n
			- "v", an integer vector in self._IVn
			- "p", a permutation corresponding to some idempotent e_p

			OUTPUT:

			- The product X^v t_{s_1}t_{s_2}...t_{s_k} e_p where word = [s_1, s_2, ..., s_k] written in the distinguished basis
			"""

			if len(word) == 0 or v == self._IVn((0,)*self._n):
				new_g = self._g_ify(word)
				new_v = v
				new_p = p
				return self.monomial(self._CP((new_v, new_g, new_p)))

			k = word[0]
			i = self._get_last_nonzero_index(v)
			l = i + 1

			right_word = word[1:]
			new_index = self._g_ify([k])(l)
			right_v = self._IVn((0,)*(new_index-1) + (1,) + (0,)*(self._n-new_index))
			left_p = self._P(self._g_ify(word[1:])(tuple(p))) #correct

			left_v = []
			for j in range(len(v)):
				if j == i:
					left_v += [v[j]-1]
				else:
					left_v += [v[j]]
			left_v = self._IVn(left_v)

			left_word = (word[0],)
			right_p = p

			first = self._switch(left_word, left_v, left_p) * self._switch(right_word, right_v, right_p)

			new_p = left_p
			if new_p[k-1] == new_p[k] and (l==k+1 or l==k):
				second_v = left_v
				second_word = word[1:]
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
		#Should be fixed... not tested yet
		@cached_method
		def _XT(self, right_g, left_v, p):
			r"""
			Returns the product of the KLR algebra elements corresponding to "left_g" and "right_v" with matching idempotents e_p

			INPUT:

			- "right_g", an element of S_n
			- "left_v", an integer vector of length n in self._IVn
			- "p", a permutation corresponding to some idempotent e_p

			OUTPUT:

			- The product X^v * t_{right_g} * e_p
			"""
			word = tuple(self.reduced_word(right_g))
			return self._switch(word, left_v, p)

		#Should be fixed... not tested yet
		@cached_method
		def product_on_basis(self, left, right):
			r"""
			Returns the product of two basis elements.

			INPUT:

			- "left", the basis element on the left
			- "right", the basis element on the right

			OUTPUT:

			- The product left * right written in the distinguished basis

			"""
			left_p = left[2]
			right_p = right[2]
			left_g = left[1]
			right_g = right[1]
			left_v = left[0]
			right_v = right[0]

			# If idempotents don't match
			if right_g.inverse()(tuple(left_p)) != tuple(right_p):
				return self.zero()

			# If right doesn't have any t's
			if right_g == self._Sn(()):
				new_v = self._add_integer_vectors(left_v, right_v)
				new_g = left_g
				new_p = right_p
				return self.monomial(self._CP((new_v, new_g, new_p)))

			# If left doesn't have any x's
			if left_v == self._IVn((0,)*self._n):
				new_v = right_v
				new_g = self._Sn(())
				new_p = right_p
				return self._TT(left_g, right_g, right_p) * self.monomial(self._CP((new_v, new_g, new_p)))

			# Finally, assume right has some t's, left has some x's.
			new_v_l = self._IVn((0,)*self._n)
			new_g_l = left_g
			new_p_l = left_p
			new_v_r = right_v
			new_g_r = self._Sn(())
			new_p_r = right_p


			return self.monomial(self._CP((new_v_l, new_g_l, new_p_l)))*(self._XT(right_g, left_v, right_p)*self.monomial(self._CP((new_v_r, new_g_r, new_p_r))))

		def plot_on_basis(self, elt, v_shift=None, h_shift=None, height_only=False, **kwargs):
			if v_shift == None:
				v = 0
			else:
				v = v_shift
			if h_shift == None:
				h = 0
			else:
				h = h_shift

			rw = self.reduced_word(elt[1])[::-1]
			iv = elt[0]
			perm = elt[2]
			n = self._n
			l = len(rw)

			max_dots = max(iv)
			# if max_dots == 0:
			# 	dot_space = 0
			# else:
			# 	dot_space = ceil((max_dots+2)/3.0)
			# #temp
			if rw == []:
				dot_space = max(1, max_dots/5.0)
			else:
				dot_space = max_dots/5.0

			if height_only:
				return dot_space+l

			color_list = rainbow(len(set(perm)))
			colors = {list(set(perm))[i]:color_list[i] for i in range(len(list(set(perm))))}

			strands = {i:[(i-1+h, v), (i-1+h, v+dot_space)] for i in range(1,n+1)}
			for i in range(l):
				strands[rw[i]], strands[rw[i]+1] = strands[rw[i]+1]+[(rw[i]-1+h, i+1+v+dot_space)], strands[rw[i]]+[(rw[i]+h, i+1+v+dot_space)]
				for j in range(1, n+1):
					if j != rw[i] and j != rw[i]+1:
						strands[j] += [(j-1+h, i+1+v+dot_space)]


			# strands = {i:"\t\\draw [color=color%s] (%s_-1.center)--(%s_0.center)--(%s_1.center)"%(perm[i-1], i, i, i) for i in range(1,n+1)}
			# for i in range(l):
			# 	strands[rw[i]], strands[rw[i]+1] = (strands[rw[i]+1]+"--(%s_%s.center)"%(rw[i], i+2), strands[rw[i]]+"--(%s_%s.center)"%(rw[i]+1, i+2))
			# 	for j in range(1, n+1):
			# 		if j != rw[i] and j != rw[i]+1:
			# 			strands[j] += "--(%s_%s.center)"%(j, i+2)
			# for j in range(1, n+1):
			# 	strands[j] += "--(%s_%s.center);\n"%(j, l+2)

			the_lines = [line(strands[i], color=colors[perm[int(strands[i][0][0]-h)]], **kwargs) for i in range(1, n+1)]


			dots = []
			for i in range(n):
				n_dots = iv[i]
				for j in range(n_dots):
					dots += [circle((i+h, -(j+1)*dot_space/(n_dots+1.0) + v + dot_space), 0.05, fill=True, color='black', **kwargs)]



			G = sum(the_lines)+sum(dots)
			G.set_aspect_ratio(1)
			G.axes(False)
			return G


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



# FOR TESTING:

# RS = RootSystem(["A", 2])
# RL = RS.root_lattice()
# alpha = RL.basis()
# a = 3*(alpha[1]+alpha[2])
# A = GaussianIntegers()
# j = A.basis()[1]
# K = KLRAlgebra(ZZ, a)
# KXT = K.XT()
# KTX = K.TX()
# #x,t,e = KTX.x,KTX.t,KTX.e
# x,t,e = KXT.x,KXT.t,KXT.e
#
# elt = 33817*x((0,1,2,3,8,30))*t(((1,6,4),(2,3)))*e((1,1,1,2,2,2))+x((1,2,0,0,0,1))*e((1,1,2,2,1,2))-4002*t(1)*e((2,1,2,1,2,1))

#elt = x((1,2,3,4))*t((1,2,3,4))*e((1,1,2,3))
#elt = t((1,2,3,4))*e((1,1,2,3))




# GOOD FOR NIL-HECKE STUFF:
'''
RS = RootSystem(["A",1])
RL = RS.root_lattice()
alpha = RL.basis()
m = 5
a = m*alpha[1]
K = KLRAlgebra(ZZ, a)
TX = K.TX()
x, t, e, f = TX.x, TX.t, TX.e, TX.divided_power_idempotent
elt = f((1,m))
'''


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
