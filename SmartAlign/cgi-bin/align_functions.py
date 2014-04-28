
import string

polarDict = { 'G':'NP', 'A':'NP', 'V':'NP', 'L': 'NP','I':'NP', 'M':'NP', 'W':'NP','F':'NP', 'P':'NP',
'S':'P', 'T':'P', 'C':'P', 'Y':'P', 'N':'P', 'Q':'P', 
'K':'PC', 'R':'PC', 'H':'PC',
'D':'PC', 'E':'PC', '-': '-'
}






def align(file):
	seq = readFastaFile(file)
	seqList = seq.values()
	polarList = polar_characterize(seqList)
	polar_SVG = createPolarSVG(polarList)
	return polar_SVG



def readFasta(file):

	lines = string.split(file, '\n')

	sequences = []
	for line in lines:

		if line[0] != '>' and line[0] != '#':
			sequences.append(line)

	return sequences

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

def polar_characterize(seqList):
	polar_lists =[]
	for sequence in seqList:
		translated_seq = []
		for character in sequence:
			translated_seq.append(polarDict[character])
		polar_lists.append(translated_seq)
	return polar_lists

def createPolarSVG(polar_list):
    '''Function to create an SVG image representing the
    open reading frame in a sequence.
    Input: length of the sequence and the result from the findORF function.
    Output: the code for the SVG image.'''

    #height of query and frames
    box_len = 20

    canvas_height = len(polar_list)*box_len + 100
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
