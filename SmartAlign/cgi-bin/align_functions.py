
import string
import operator
import math 


polarDict = { 'G':'NP', 'A':'NP', 'V':'NP', 'L': 'NP','I':'NP', 'M':'NP', 'W':'NP','F':'NP', 'P':'NP',
'S':'P', 'T':'P', 'C':'P', 'Y':'P', 'N':'P', 'Q':'P', 
'K':'PC', 'R':'PC', 'H':'PC',
'D':'NC', 'E':'NC', '-': '-'
}
polar_ops = 4


essentialDict = {'A':'N', 'R':'E', 'N':'N', 'D':'N',
'C':'N', 'E':'N', 'Q':'N', 'G':'N',
'H':'E', 'I':'E', 'L':'E', 'K':'E',
'M':'E', 'F':'E', 'P':'N', 'S':'N', 
'T':'E', 'W':'E', 'Y':'N', 'V':'E', '-':'-'}
ess_ops = 2

ketoGlycloDict = { 'A':'G', 'R':'G', 'N':'G', 'D':'G',
'C':'G', 'E':'G', 'Q':'G', 'G':'G',
'H':'G', 'I':'B', 'L':'K', 'K':'K',
'M':'G', 'F':'B', 'P':'G', 'S':'G',
'T':'B', 'W':'B', 'Y':'B', 'V':'G', '-': '-' }
kg_ops = 3

def align(file):
	char_dict = ketoGlycloDict
	ops = kg_ops

	seq = readFastaFile(file)
	seqList = seq.values()
	#characterize list
	char_lst = characterize(seqList, char_dict)
	#index list
	sorted_list = index_list(char_lst)
	
	sorted_SVG = createSortedSVG(sorted_list, char_dict, ops)
	#create list of index_lists

	#polar_SVG = createPolarSVG(polarList)
	return sorted_SVG

def index_list(input):
	num_seqs = len(input)
	len_seq = len(input[0])
	i_list = []
	for i in xrange(len_seq):
		i_elt = {}
		for y in xrange(num_seqs):
			if input[y][i] in i_elt:
				i_elt[ input[y][i] ] += 1
			else:
				i_elt[ input[y][i]] = 1 
		sorted_elt = sorted(i_elt.iteritems(), key=operator.itemgetter(1))
		i_list.append(sorted_elt)
	return i_list


def readFastaFile(filename):
    '''Read sequences from fasta file and returns a dictionary with the sequence names as keys and the sequences as values
    Input should be the name of the file. Note that this function will open the file.'''
   
    filein = open(filename,"r")
    
    # We'll read the first line and check if it starts with '>' if it
    # does we assign this to the first name and continue. If not we return
    # an error message.
    line = filein.readline()
    if line[0] != '>':
        return "File should start with a sequence in the first line.\n"
    seqs = {}    # create an empty dictionary
    while line and line[0] == '>':
        name = line.lstrip(">").strip()  # remove '>' in beginning of sequence could also be name = line[1:]
        # the second strip removes newlines and whitespaces
        line = filein.readline()
        # we have the name of the first sequence, let's get the sequence
        # we can read the file, line-by-line until we reach a line that
        # starts with '>' that signifies the end of the sequence
        seq = ''
        while line and line[0] != '>':
            seq = seq + line.strip()
            line = filein.readline()
        seqs[name]=seq
    return seqs

def characterize(seqList, char_dict):
	lsts =[]
	for sequence in seqList:
		translated_seq = []
		for character in sequence:
			translated_seq.append(char_dict[character])
		lsts.append(translated_seq)
	return lsts

def createPolarSVG(polar_list):
    '''Function to create an SVG image representing the
    open reading frame in a sequence.
    Input: length of the sequence and the result from the findORF function.
    Output: the code for the SVG image.'''

    #height of query and frames
    box_len = 5

    canvas_height = len(polar_list)*box_len + 4000
    canvas_width = len(polar_list[0])*box_len + 50
 

    #start SVG string and draw the query rectangle with start and end coordinates
    figureSVG = '''
<svg height="{height}" width="{width}" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">'''.format(height = canvas_height, width = canvas_width)
    #figureSVG = figureSVG + 
    x_pos = 0
    y_pos = 0
    seq_SVG = ''

    for sequence in polar_list:
    	for character in sequence:

    		if character == 'NP':
    			seq_SVG = seq_SVG + """<rect fill = "blue" height = "{box_len}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "blue" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (40 + y_pos*box_len), box_len = box_len)
    		elif character == 'P':
    			seq_SVG = seq_SVG + """<rect fill = "red" height = "{box_len}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "red" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (40 + y_pos*box_len), box_len = box_len)
    		elif character == 'PC':
    			seq_SVG = seq_SVG + """<rect fill = "orange" height = "{box_len}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "orange" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (40 + y_pos*box_len), box_len = box_len)
    		elif character == 'NC':
    			seq_SVG = seq_SVG + """<rect fill = "green" height = "{box_len}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "green" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (40 + y_pos*box_len), box_len = box_len)
    		elif character == '-':
    			seq_SVG = seq_SVG + """<rect fill = "black" height = "{box_len}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (40 + y_pos*box_len), box_len = box_len)
    		x_pos += 1

    	x_pos = 0
    	y_pos += 1
	figureSVG = figureSVG + seq_SVG

	text_SVG = ''

	for i in xrange(len(polar_list[0])):
		if i % 10 == 0:
			text_SVG = text_SVG + """<text x ="{text_X}" y ="{y_pos}" font-family = "Serif" font-size = "8" fill = "black"> {num} </text>""".format(text_X = (20 + i*box_len + box_len/(4.0)), y_pos = (30 + (y_pos + 1)*box_len), num = i )

	figureSVG = figureSVG + text_SVG

    figureSVG += '</svg>'
    return figureSVG

def createSortedSVG(input, char_dict, numOps):

	#height of query and frames
	box_len = 15
	box_height = 100

	canvas_height = len(input[0])*box_len + 400
	canvas_width = len(input)*box_len + 100
 

    #start SVG string and draw the query rectangle with start and end coordinates
	figureSVG = '''
<svg height="{height}" width="{width}" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">'''.format(height = canvas_height, width = canvas_width)
	x_pos = 0
	y_pos = canvas_height 
	seq_SVG = ''
	for index in input: 
		#print 'index:', index
		for character in index:
			if character[0] == 'NP':
				seq_SVG = seq_SVG + """<rect fill = "blue" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'P':
				seq_SVG = seq_SVG + """<rect fill = "red" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'PC':
				seq_SVG = seq_SVG + """<rect fill = "orange" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'NC':
				seq_SVG = seq_SVG + """<rect fill = "green" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'E':
				seq_SVG = seq_SVG + """<rect fill = "red" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'N':
				seq_SVG = seq_SVG + """<rect fill = "blue" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'G':
				seq_SVG = seq_SVG + """<rect fill = "blue" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'B':
				seq_SVG = seq_SVG + """<rect fill = "green" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			elif character[0] == 'K':
				seq_SVG = seq_SVG + """<rect fill = "yellow" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character, numOps)), box_len = box_len, height = box_height*find_height(char_dict, index, character, numOps))
				y_pos = y_pos - box_height*find_height(char_dict, index, character, numOps)
			#elif character[0] == '-':
				#seq_SVG = seq_SVG + """<rect fill = "black" height = "{height}" width = "{box_len}" x ="{xcoord}" y ="{ycoord}" stroke = "black" />""".format(xcoord = (20 + x_pos*box_len), ycoord = (y_pos - box_height*find_height(char_dict, index, character)), box_len = box_len, height = box_height*find_height(char_dict, index, character))
				#y_pos = y_pos - box_height*find_height(char_dict, index, character)
		y_pos = canvas_height
		x_pos += 1
	
	figureSVG = figureSVG + seq_SVG
	# text_SVG = ''

	# for i in xrange(len(polar_list[0])):
	# 	if i % 10 == 0:
	# 		text_SVG = text_SVG + """<text x ="{text_X}" y ="{y_pos}" font-family = "Serif" font-size = "8" fill = "black"> {num} </text>""".format(text_X = (20 + i*box_len + box_len/(4.0)), y_pos = (30 + (y_pos + 1)*box_len), num = i )

	# figureSVG = figureSVG + text_SVG

	figureSVG += '</svg>'
	return figureSVG

def find_height(char_dict, index, character, numOps):
	max_height =math.log( numOps, 2)
	num_seq = 0
	#print index
	for symbol in index:
		num_seq += symbol[1]

	rel_freq = float(character[1])/float(num_seq)
	entro = rel_freq*math.log(rel_freq, 2.0)
	R = rel_freq*(max_height + entro)
	return R

	#return max_height - float(character[1])/float(num_seq) * (-1) *float(math.log(float(character[1])/float(num_seq), 2))
