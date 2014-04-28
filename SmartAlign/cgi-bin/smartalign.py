import cgi
import cgitb
import align_functions


cgitb.enable()

print '''Content-Type: text/html

'''
form = cgi.FieldStorage()


fasta_inputs = form.getvalue("input_alignments")


fasta_file = form.getvalue("datafile")



#query = gbfunctions.queryDB(genename,chrom,coord)
    
figure = align_functions.align(fasta_file)

print '''
<!DOCTYPE html>

<html>
   <head>
       <title>Smart Align</title>
       
       <link type="text/css" rel="stylesheet" href="/default.css" media="all" />  
       <link href="/fonts.css" rel="stylesheet" type="text/css" media="all" />

   </head>
   <body bgcolor = "purple">
           
           <center>
           	<h2> Side Chain Charge Characterization </h2>
           	Red = Polar | Blue = Non-Polar  | Orange = Positively Charged | Green = Negatively Charged 
			<div style="width:80%;height:150px;line-height:3em;overflow: auto; padding:5px; background-color:#FFFFFF;border:20px double #EE82EE">
					{figure} 
			</div>

          </center>
   </body>
</html>
'''.format(figure=figure)

