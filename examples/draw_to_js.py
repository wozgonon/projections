#!/bin/python
import fileinput

print("""
var c=document.getElementById("canvas");
var ctx=c.getContext("2d");
ctx.beginPath();
ctx.save(); 
ctx.strokeStyle = 'blue';
ctx.lineWidth = 1.0;


function draw_all (min_x, min_y, max_x, max_y) {

ww = max_y-min_y;
hh = max_x-max_y;

length =  canvas.width > canvas.height ? canvas.height : canvas.width;
x_scale  = length / ww;
y_scale  = length / hh;
scale = x_scale > y_scale ? x_scale : y_scale;

x_offset = min_x * scale;
y_offset = min_y * scale;
//x_offset = (canvas.width  - ww * scale)/2
//y_offset = (canvas.height - hh * scale)/2

ctx.beginPath();

ctx.fillText("scale=" + scale, (canvas.width / 2) - 17, (canvas.height / 2) - 18);
ctx.fillText("x_scale=" + y_scale, (canvas.width / 2) - 17, (canvas.height / 2) +20 );
ctx.fillText("y_scale=" + x_scale, (canvas.width / 2) - 17, (canvas.height / 2) +60);
ctx.fillText("x_offset=" + y_offset, (canvas.width / 2) - 17, (canvas.height / 2) +100);
ctx.fillText("y_offset=" + y_offset, (canvas.width / 2) - 17, (canvas.height / 2) +140);

""")

for line in fileinput.input():
    columns=line.rstrip().split ('\t')
    if len(columns) != 3:
        print("Expect three columns, not: " + line)
    else:
        command=columns[0]
        xx=columns[1]
        yy=columns[2]
        if command == "drawto":
            command = "ctx.lineTo"
        elif command == "moveto":
            command = "ctx.moveTo"
        elif command == "topleft":
            command = "ctx.moveTo"
            min_x=xx
            min_y=yy
        elif command == "bottomright":
            command = "ctx.moveTo"
            max_x=xx
            max_y=yy
        
        if xx != "inf" and yy !="inf" and xx != "-inf" and yy !="-inf":
#            print("%s(%s*3000-8800,%s*3000-5340);" % (command, xx, yy))
            print("%s(%s*scale-x_offset,%s*scale-y_offset);" % (command, xx, yy))

print("""
}
""")

print("draw_all(%s,%s,%s,%s)" % (min_x,min_y,max_x,max_y))

print("""
         
ctx.stroke();
ctx.restore();
""")

      
