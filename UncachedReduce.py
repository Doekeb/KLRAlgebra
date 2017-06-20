def reduced_words(g):
	return tuple(g.reduced_word())

def reduce_full(word):
	return reduce(word+reduced_words(g_ify(word))[::-1]+reduced_words(g_ify(word)))

'''
Input can be any list containing positive integer entries. Any such list can be converted
to an element of S_n (for appropriate n) by multiplying the corresponding transpositions.
The result should be another such list which is a minimal length word for the same S_n
element.
'''
def reduce(word):
	
	# Base case: if my word is already reduced.
	if g_ify(word).length() == len(word):
		return word
	
	
	# I am looking for a minimal unreduced subword. Start at the beginning.
	i = 0
	# Increase the index until the subword starting at the beginning is unreduced.
	while g_ify(word[:i+1]).length() == len(word[:i+1]):
		i += 1
	
	
	# Now go backwards...
	j = i
	# Decrease this second index until the subword starting at j ending at i is unreduced.
	while g_ify(word[j:i+1]).length() == len(word[j:i+1]):
		j -= 1
	
	
	
	# Another semi-base case: If my minimal unreduced word has length 2 (and so it is a
	# double-crossing), just get rid of this piece.
	if i - j == 1:
		return reduce(word[:j]+word[i+1:])
	
	
	
	# Another semi-base case: If my minimal unreduced word can be reduced by swapping
	# distant transpositions, do such a swap.
	if abs(word[i] - word[i-1]) > 1:
		return reduce(word[:i-1]+(word[i], word[i-1])+word[i+1:])
	
	
	
	# Now I have guaranteed that my minimal unreduced word ends in two adjacent
	# transpositions, and a consequence of minimality is that the two strands which appear
	# in the last two transpositions but don't cross there, DO cross somewhere else in the
	# minimal unreduced word, before position 0. Straighten this piece.
	return reduce(word[:j+1] + straighten(word[j+1:i+1]) + word[i+1:])




'''
Input should represent a braid diagram where the last pair of transpositions are adjacent,
the length is at least 3, and the two strands which appear in the last two transpositions
but don't cross there, DO cross somewhere else.
'''
def straighten(word):
	
	# Base case: If I have something like (i)(i+1)(i), then apply the braid relation.
	if word[-1] == word[-3]:
		return word[:-3]+(word[-2], word[-1], word[-2])
	
	
	# I am guaranteed to have those two distinguished strands cross SOMEWHERE, find that
	# location
	j = find_crossing(word)
	
	
	# Another semi-base case: If that crossing can be moved upward in the diagram, then
	# do it.
	if abs(word[j] - word[j+1]) > 1:
		return straighten(word[:j]+(word[j+1], word[j])+word[j+2:])
	
	
	# Now we are guaranteed to be in an upside-down situation from what we started with
	# BUT WITH SHORTER LENGTH. So flip this shorter subword, straighten it, then flip it
	# back.
	return word[:j]+straighten(word[j:][::-1])[::-1]



'''
Input should be the same as in the method "straighten" above. Returns the index at which
those two distinguished strands cross.
'''
def find_crossing(word):

	# Make a word with an extra crossing at the top, in particular a crossing of the two
	# strands of interest
	temp = word+(word[-2],)
	
	
	# We'll be looking for the index at which those two strands cross again, which is the
	# same as the latest index at which we have an unreduced subword on the tail.
	j = len(word)
	while g_ify(temp[j:]).length() == len(temp[j:]):
		j -= 1
	return j


'''
Takes a list of positive integers and turns it into an S_n element by multiplying the
appropriate transpositions.
'''
def g_ify(word):
	
	# Make a symmetric group of the appropriate size.
	n = 0
	if word != ():
		n = max(word)
	Sn = SymmetricGroup(n+1)
	
	# Multiply the appropriate transpositions and return the result.
	return Sn.prod([Sn((i, i+1)) for i in word])