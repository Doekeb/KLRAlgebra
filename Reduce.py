@cached_function
def reduced_words(g):
	return tuple(g.reduced_word())

# @cached_function
# def reduce_full(word):
# 	return reduce(word+reduced_words(g_ify(word))[::-1]+reduced_words(g_ify(word)))

'''
Input can be any list containing positive integer entries. Any such list can be converted
to an element of S_n (for appropriate n) by multiplying the corresponding transpositions.
The result should be another such list which is a minimal length word for the same S_n
element.
'''
@cached_function
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

def reduce_full(word):
	# If we aren't reduced yet, reduce first
	if len(word) != g_ify(word).length():
		return reduce_full(reduce(word))

	target = reduced_words(g_ify(word))

	word = list(word)
	my_move = get_move(tuple(word), target)
	print my_move
	while my_move[0] != "n":
		if my_move[0] == "s":
			i = my_move[1]
			assert abs(word[i]-word[i+1]) > 1
			word[i], word[i+1] = word[i+1], word[i]
			#print tuple(word)
		elif my_move[0] == "b":
			i = my_move[1]
			assert word[i] == word[i+2] and abs(word[i+1] - word[i]) == 1
			word[i], word[i+1], word[i+2] = word[i+1], word[i], word[i+1]
			#print tuple(word)
		my_move = get_move(tuple(word), target)
		print my_move

	return tuple(word)

def get_move(word, target):
	assert g_ify(word) == g_ify(target)
	if word == target:
		return ["n", 0]

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
		my_move = get_move(word[i:], target[i:])
		my_move[1] += i
		return my_move

	# Now the beginning doesn't match. Find where the target first crossing
	# crosses in word.
	temp = (target[0],) + word
	j = 0
	while g_ify(temp[:j]).length() == j:
		j += 1
	j -= 2
	# j is now the index at which the distinguished strands cross in word.
	return float_move(word[:j+1])

'''
give a move which floats the bottom crossing to the top
'''
def float_move(word):
	# if the two bottom crossings are far apart, just move them past
	# each other for free. Do this until you can't anymore or until the crossing
	# is all the way at the top.
	j = len(word)-1
	if j > 0 and abs(word[j] - word[j-1]) > 1:
		return ['s', j-1]

	# Now we can assume that the crossing of interest is as far up as possible.
	# If it's at the top, then we are done with this part, just go back to
	# reduce_full
	# if j == 0:
	# 	return word

	# Now we can assume it got stuck on some adjacent transposition.
	#print word, j
	i = find_crossing(word)
	if i == j-2:
		return ['b', j-2]
	my_move = sink_move(word[i:])
	my_move[1] += i
	return my_move

'''
give a move which floats the top crossing to the bottom
'''
def sink_move(word):
	my_move = float_move(word[::-1])
	if my_move[0] == "s":
		my_move[1] = len(word) - 2 - my_move[1]
	elif my_move[0] == "b":
		my_move[1] = len(word) - 3 - my_move[1]
	return my_move

'''
floats the bottom crossing to the top
'''
@cached_method
def float(word):
	# if the two bottom crossings are far apart, just move them past
	# each other for free. Do this until you can't anymore or until the crossing
	# is all the way at the top.
	j = len(word)-1
	while j > 0 and abs(word[j] - word[j-1]) > 1:
		word = list(word)
		word[j], word[j-1] = word[j-1], word[j]
		word = tuple(word)
		j -= 1

	# Now we can assume that the crossing of interest is as far up as possible.
	# If it's at the top, then we are done with this part, just go back to
	# reduce_full
	if j == 0:
		return word

	# Now we can assume it got stuck on some adjacent transposition.
	#print word, j
	i = find_crossing(word[:j+1])
	assert i < j-1
	word = word[:i] + sink(word[i:j-1]) + word[j-1:]
	assert word[j] == word[j-2]
	assert word[j-1] == word[j]-1 or word[j-1] == word[j]+1
	word = list(word)
	word[j], word[j-1], word[j-2] = word[j-1], word[j], word[j-1]
	word = tuple(word)

	return float(word[:j-1])+word[j-1:]


'''
preamble is just for use in printing the entire algorithm
'''
@cached_method
def OLDreduce_full(word, target=None):

	# If we aren't reduced yet, reduce first
	if len(word) != g_ify(word).length():
		return reduce_full(reduce(word))

	# Now we may assume that word is reduced (but maybe not the distinguished
	# one)

	# If a target isn't specified, use the one given by the hardcoded function
	if target == None:
		target = reduced_words(g_ify(word))

	# Make sure word and target represent the same Sn element
	assert g_ify(word) == g_ify(target)

	# Base case... nothing to do (note: this handles the case where the
	# corresponding Sn element is the identity. Both word and target are empty
	# in that case)
	if word == target:
		return word

	# Check to see how far we match the target
	i = 0
	try:
		while word[i] == target[i]:
			i += 1
	# Error is thrown if nothing matches. We don't care.
	except IndexError:
		pass

	# If we match anything, just reduce the leftovers. The rest of the method
	# will be a long case corresponding to no matches.
	if i > 0:
		return word[:i] + reduce_full(word[i:], target=target[i:])

	# Now we assume that the first (top) index DOESN'T match.

	# Look at target to see which strands should cross first,
	# find where they cross.
	temp = (target[0],) + word
	j = 0
	while g_ify(temp[:j]).length() == j:
		j += 1

	j -= 1
	# Use this other method to move that crossing all the way up
	return reduce_full(float(word[:j]) + word[j:], target=target)


'''
sinks the top crossing to the bottom
'''
@cached_method
def sink(word):
	return float(word[::-1])[::-1]



'''
Input should represent a braid diagram where the bottom pair of transpositions are adjacent,
the length is at least 3, and the two strands which appear in the bottom two transpositions
but don't cross there, DO cross somewhere else.
'''
@cached_function
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
@cached_function
def find_crossing(word):

	# Make a word with an extra crossing at the bottom, in particular a crossing of the two
	# strands of interest
	temp = word+(word[-2],)


	# We'll be looking for the index at which those two strands cross again, which is the
	# same as the latest index at which we have an unreduced subword on the tail.
	j = len(word)
	while g_ify(temp[j:]).length() == len(temp[j:]):
		j -= 1
	return j

'''
Input should represent a braid diagram where the top pair of transpotitions are
adjacent and the two strands in the top pair which DON'T cross there, DO cross
somewhere else.
'''
def find_crossing_reverse(word):
	temp = (word[1],)+word
	j = 0
	while g_ify(temp[:j]).length() == j:
		j += 1
	return j-2

'''
Takes a list of positive integers and turns it into an S_n element by multiplying the
appropriate transpositions.
'''
@cached_function
def g_ify(word):

	# Make a symmetric group of the appropriate size.
	n = 0
	if word != ():
		n = max(word)
	Sn = SymmetricGroup(n+1)

	# Multiply the appropriate transpositions and return the result.
	return Sn.prod([Sn((i, i+1)) for i in word])
