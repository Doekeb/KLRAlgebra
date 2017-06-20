import sage
from sage.rings.ring import Algebra
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import AlgebraElement
from sage.combinat.free_module import CombinatorialFreeModuleElement

def add_integer_vectors(first, second):
	if first.parent() != second.parent():
		raise TypeError("%s has parent %s but %s has parent %s"%(first, first.parent(), second, second.parent()))
	if len(first) != len(second):
		raise ValueError("%s has length %s but %s has length %s"%(first, len(first), second, len(second)))
	P = first.parent()
	return P((first[i]+second[i] for i in range(len(first))))

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
	
	def tikz(self, x_sep=0.5, y_sep=0.3, coeff=1, first_term=True):
		if self.is_monomial():
			rw = self._get_reduced_word()
			n = self._n
			l = len(rw)
			front_matter = "\\begin{tikzpicture}\n\t\\foreach \\x in {1, ..., %s}\n\t\t\\foreach \\y in {0, ..., %s}\n\t\t\t\\node (\\x_\\y) at (%s*\\x, %s*\\y){};\n"%(n, l+2, x_sep, y_sep)
			if coeff == 1:
				coeff = ""
			else:
				coeff = "\\,\\,%s"%coeff
			if first_term:
				end_matter = "\t\\node[left=4pt of %s]{$%s$};\n\\end{tikzpicture}"%("1_%s"%(floor((l+2)/2)), coeff)
			else:
				end_matter = "\t\\node[left=4pt of %s]{$+%s$};\n\\end{tikzpicture}"%("1_%s"%(floor((l+2)/2)), coeff)
			
			simple_roots = self.parent().root_system().root_lattice().simple_roots()
			color_list = rainbow(len(simple_roots))
			colors = {i:color_list[i-1] for i in range(len(simple_roots))}
			perm = self.get_permutation()
			color_string = ""
			for i in colors:
				color_string += "\t\\definecolor{color%s}{HTML}{%s}\n"%(i+1, colors[i][1:].upper())
			
			strands = {i:"\t\\draw [color=color%s] (%s_0.center)--(%s_1.center)"%(perm[i-1], i, i) for i in range(1,n+1)}
			for i in range(l):
				strands[rw[i]], strands[rw[i]+1] = (strands[rw[i]+1]+"--(%s_%s.center)"%(rw[i], i+2), strands[rw[i]]+"--(%s_%s.center)"%(rw[i]+1, i+2))
				for j in range(1, n+1):
					if j != rw[i] and j != rw[i]+1:
						strands[j] += "--(%s_%s.center)"%(j, i+2)
			for j in range(1, n+1):
				strands[j] += "--(%s_%s.center);\n"%(j, l+2)
			
			the_string = front_matter+color_string
			for i in strands:
				the_string += strands[i]
			the_string += end_matter
			return the_string
		the_string = ""
		mcs = self.monomial_coefficients()
		first_term=True
		for monomial in mcs:
			the_string += "\n"+self.parent().monomial(monomial).tikz(coeff=mcs[monomial], first_term=first_term)
			first_term=False
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
	
	def get_permutation(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return self.leading_support()[2]
		
	def get_Sn_element(self):
		if not self.is_monomial():
			raise ValueError("%s is not a monomial"%self)
		return self.leading_support()[1]
	
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
	Element = KLRAlgebraElement
	def __init__(self, base_ring, root_system, positive_root, reduced_words=None, signs=None, **kwds):
		self._root_system = root_system
		self._positive_root = positive_root
		self._n = positive_root.height()
		self._Sn = SymmetricGroup(self._n)
		self._name = "KLR Algebra over %s for %s with positive root %s"%(repr(self.base()), repr(self._root_system), repr(self._positive_root))
		if reduced_words != None:
			self._name += " with custom reduced words"
			if signs != None:
				self._name += " and sign conventions"
		elif signs != None:
			self._name += " with custom sign conventions"
		
		#self._reduced_words = reduced_words or {g:g.reduced_word() for g in self._Sn}
		
		temp_reduced_words = reduced_words or (lambda g: g.reduced_word())
		self._reduced_words = lambda g: temp_reduced_words(self._Sn(g))
		
		self._signs = signs or (lambda i, j: sign(j-i))
		
		from sage.combinat.integer_vector import IntegerVectors_k
		self._IVn = IntegerVectors_k(self._n)
		
		L = positive_root.dense_coefficient_list()
		IndexList = []
		for i in range(len(L)):
			IndexList += [i+1]*L[i]
		
		self._P = Permutations(IndexList)
		
		self._CP = cartesian_product([self._IVn, self._Sn, self._P])
		
		CombinatorialFreeModule.__init__(self, base_ring, self._CP, category=GradedAlgebrasWithBasis(base_ring), **kwds)
	
	def _repr_(self):
		return self._name
	
	def base_ring(self):
		return self.base().base_ring()
	
	def characteristic(self):
		return self.base_ring().characteristic()
				
	def root_system(self):
		return self._root_system
	
	def positive_root(self):
		return self._positive_root
	
	def symmetric_group(self):
		return self._Sn
		
	def reduced_word(self, g=None):
		if g == None:
			return self._reduced_words
		if g not in self._Sn:
			raise ValueError("%s not an element of %s"%(g, self._Sn))
		return self._reduced_words(self._Sn(g))
	
	def x(self, i=None):
		if i == None:
			return {j:self.x(j) for j in range(1, self._n+1)}
		if i not in range(1, self._n+1):
			raise ValueError("%s not an integer between 1 and %s"%(i, self._n))
		#B = self.basis()
		v = self._IVn((0,)*(i-1)+(1,)+(0,)*(self._n-i))
		g = self._Sn(())
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
	
	def polynomial(self, v=None):
		if v == None:
			raise NotImplementedError("Must specify an integer vector for now")
		if v not in self._IVn:
			raise ValueError("%s not an appropriate integer vector"%v)
		v = self._IVn(v)
		g = self._Sn(())
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
		
	def t(self, g=None):
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
		#B = self.basis()
		g = self._Sn(())
		v = self._IVn((0,)*self._n)
		return self.sum([self.monomial(self._CP((v, g, p))) for p in self._P])
		
	def gen(self, i=0):
		pass
		
	def algebra_generators(self):
		return [self.x(i) for i in range(1, self._n+1)]+[self.t(i) for i in range(1, self._n)]+[self.e(p) for p in self._P]
	
	def _t_times_x(self, t, x):
		k = t.is_t()
		l = x.is_x()
		if not (k and l):
			raise ValueError("%s should be a t and %s should be a x"%(t, x))
		sum = self.zero()
		#print "%s first"%sum
		for i in self._P:
			#print i
			g = self._Sn((k, k+1))
			identity = self._Sn(())
			zero = self._IVn((0,)*self._n)
			#Read: i[k] == i[k+1]
			if i[k-1] == i[k] and l == k+1:
				v = self._IVn((0,)*(k-1) + (1,) + (0,)*(self._n-k))
				sum += self.monomial(self._CP((v, g, i))) + self.monomial(self._CP((zero, identity, i)))
				#print "%s second"%sum
			elif i[k-1] == i[k] and l == k:
				v = self._IVn((0,)*(k) + (1,) + (0,)*(self._n-k-1))
				sum += self.monomial(self._CP((v, g, i))) - self.monomial(self._CP((zero, identity, i)))
				#print "%s third"%sum
			else:
				v = self._IVn((0,)*(l-1)+(1,)+(0,)*(self._n-l))
				sum += self.monomial(self._CP((v, g, i)))
				#print "%s fourth"%sum
		return sum
	
	#I THINK THIS IS THE LAST PIECE I NEED TO WRITE			
	def _simplify_t(self, word):
		g = prod([self._Sn((i, i+1)) for i in word])
		if word in g.reduced_words():
			#print "here I am"
			# !!!!!!WARNING!!!!!!THIS GIVES THE WRONG ANSWER!!!!!!IT'S JUST A PLACEHOLDER!!!!!!
			# I still need to move toward the DISTINGUISHED reduced word.
			return self.t(g)
		i = 2
		sursw = word[-i:] #shortest un-reduced sub-word (starting from the beginning)
		reduced_words = prod([self._Sn((i,i+1)) for i in sursw]).reduced_words()
		while True:
			if sursw not in reduced_words:
				break
			i += 1
			sursw = word[-i:]
			reduced_words = prod([self._Sn((i,i+1)) for i in sursw]).reduced_words()
		#THIS IS WHERE I LEFT OFF
		
	#SHOULD WORK ONCE _simplify_t IS WRITTEN
	def product_on_basis(self, left, right):
		# If idempotents don't match
		if self._P(right.get_Sn_element().inverse()(tuple(left.get_permutation()))) != right.get_permutation():
			return self.zero()
		# Now we may assume idempotents match.
		# If left doesn't have a tau
		if left.get_Sn_element() == self._Sn(()):
			vl = left.get_integer_vector()
			vr = right.get_integer_vector()
			#new_v = self._IVn([vl[i]+vr[i] for i in range(self._n)])
			new_v = add_integer_vectors(vl, vr)
			return self.monomial(self._CP((new_v, right.get_Sn_element(), right.get_permutation())))
		# Now we may assume left has a tau
		# If right doesn't have any x's
		if right.get_first_nonzero_index() == -1:
			# If also right doesn't have any taus (so right is just an e)
			if right.get_Sn_element() == self._Sn(()):
				return left
			
			# Now we may assume right has some taus (so does left, but right still has no x's)
			left_word = left._get_reduced_word()
			right_word = right._get_reduced_word()
			v = left.get_integer_vector()
			new_left = self.sum([self.monomial(self._CP((v, self._Sn(()), p))) for p in self._P])
			return new_left*self._simplify_t(right_word+left_word)*self.e(right.get_permutation())
					
		# Now we may assume right has some x's
		old_left_g = left.get_Sn_element()
		first_reflection_index = self.reduced_word(old_left_g)[0]
		first_reflection_t = self.t(first_reflection_index)
		new_left_g = prod([self._Sn((i, i+1)) for i in self.reduced_word(old_left_g)[1:]])
		left_v = left.get_integer_vector()
		
		new_left = self.sum([self.monomial(self._CP((left_v, new_left_g, p))) for p in self._P])
		
		old_right_v = right.get_integer_vector()
		index = right.get_first_nonzero_index()
		new_right_v = []
		for i in range(self._n):
			if i == index:
				new_right_v += [old_right_v[i]-1]
			else:
				new_right_v += [old_right_v[i]]
		new_right_v = self._IVn(new_right_v)
		right_g = right.get_Sn_element()
		right_p = right.get_permutation()
		
		new_right = self.monomial(self._CP((new_right_v, right_g, right_p)))
		
		return new_left*self._t_times_x(first_reflection_t, self.x(index+1))*new_right
		
	def product(self, left, right):
		sum = self.zero()
		for left_term in left.monomial_coefficients().keys():
			for right_term in right.monomial_coefficients().keys():
				sum += left[left_term]*right[right_term]*self.product_on_basis(self.monomial(left_term), self.monomial(right_term))
		return sum
	
	#WORKING ON THIS TOO
	def degree_on_basis(self, elt):
		return 0
	
#	def ngens(self):
#		return 3

#class KLRAlgebraElement(AlgebraElement):
	
#	def __init__(self)
		
# DO I NEED THIS? IF SO, DO "import __main__"
# __main__.KLRAlgbra=KLRAlgebra

# FOR TESTING:
RS = RootSystem(["A", 20])
RL = RS.root_lattice()
alpha = RL.basis()
a = 7*alpha[1]
#a = alpha[1]+2*alpha[3]
#a = alpha[1]+2*alpha[3]+alpha[2]+alpha[5]
#a = alpha[1]+alpha[2]+alpha[3]+alpha[4]+alpha[5]+alpha[6]
K = KLRAlgebra(ZZ, RS, a)
x, t, e = (K.x, K.t, K.e)
z = K.an_element()
#a = K.monomial(K._CP((K._IVn((0, 0, 0)), K._Sn((1, 2)), K._P((1, 3, 3)))))
#b = K.monomial(K._CP((K._IVn((1, 0, 0)), K._Sn(()), K._P((1, 3, 3)))))
#c = K.monomial(K._CP((K._IVn((0, 1, 0)), K._Sn((1, 2)), K._P((1, 3, 3)))))
