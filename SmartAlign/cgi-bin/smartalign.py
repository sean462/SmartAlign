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
    
figure_lst = align_functions.align(fasta_file)


output =  '''
<!DOCTYPE html>

<html>
   <head>
       <title>Smart Align</title>
       
       <link type="text/css" rel="stylesheet" href="/default.css" media="all" />  
       <link href="/fonts.css" rel="stylesheet" type="text/css" media="all" />

   </head>
   <body>
           
           <center>
           	<h2> Side Chain Charge Characterization </h2>
           	Red = Polar | Blue = Non-Polar  | Orange = Positively Charged | Green = Negatively Charged 
			<div style="width:80%;height:350px;line-height:3em;overflow: auto; overflow-y: hidden;
 padding:5px; background-color:#FFFFFF;border:20px double #000000;">
					{figure1} 
          </div>

          <p>
          <h2> Essential vs. Nonessential Characterization </h2>
            Red = Essential | Blue = Nonessential 
          <div style="width:80%;height:350px;line-height:3em;overflow: auto; overflow-y: hidden;
 padding:5px; background-color:#FFFFFF;border:20px double #000000">
          {figure2}
          </div>

          <p>
          <h2> Ketogenic vs. Glycogenic Characterization </h2>
            Blue = Glycogenic | Yellow = Ketogenic | Green = Both   
          <div style="width:80%;height:350px;line-height:3em;overflow: auto; overflow-y: hidden;
 padding:5px; background-color:#FFFFFF;border:20px double #000000">
          {figure3}
			</div>

          </center>
   </body>
</html>
'''.format(figure1=figure_lst[0], figure2 = figure_lst[1], figure3 = figure_lst[2])

print output

