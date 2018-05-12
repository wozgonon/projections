#!/bin/python
import fileinput

print("""
var c=document.getElementById("canvas");
var ctx=c.getContext("2d");
ctx.beginPath();
ctx.save(); 
ctx.strokeStyle = 'blue';
ctx.lineWidth = 1.0;
ctx.beginPath();
""")

for line in fileinput.input():
    columns=line.rstrip().split ('\t')
    if len(columns) != 3:
        print("Expect three columns, not: " + line)
    else:
        command=columns[0]
        if command == "drawto":
            command = "ctx.lineTo"
        elif command == "moveto":
            command = "ctx.moveTo"
        
        xx=columns[1]
        yy=columns[2]
        if xx != "inf" and yy !="inf" and xx != "-inf" and yy !="-inf":
            print("%s(%s*3000-8800,%s*3000-5340);" % (command, xx, yy))
    
print("ctx.stroke();")
print("ctx.restor();")

      
