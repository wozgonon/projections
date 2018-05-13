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

    ww = max_x-min_x;
    hh = max_y-min_y;

    length =  canvas.width > canvas.height ? canvas.height : canvas.width;
    x_scale  = length / ww;
    y_scale  = length / hh;
    scale = x_scale < y_scale ? x_scale : y_scale;

    xo=(canvas.width - ww * scale)/2;
    yo=(canvas.height - hh * scale)/2;

    x_offset = min_x*scale - xo;
    y_offset = min_y*scale - yo;

    ctx.beginPath();

    ctx.fillText("ww*scale=" + ww*scale,       20, canvas.height - 165);
    ctx.fillText("hh*scale=" + hh*scale,       20, canvas.height - 150);
    ctx.fillText("xo="       + xo,             20, canvas.height - 135);
    ctx.fillText("yo="       + yo,             20, canvas.height - 120);
    ctx.fillText("width="    + canvas.width,   20, canvas.height - 105);
    ctx.fillText("height="   + canvas.height,  20, canvas.height - 90);
    ctx.fillText("scale="    + scale,          20, canvas.height - 75);
    ctx.fillText("x_scale="  + y_scale,        20, canvas.height - 60);
    ctx.fillText("y_scale="  + x_scale,        20, canvas.height - 45);
    ctx.fillText("x_offset=" + x_offset,       20, canvas.height - 30);
    ctx.fillText("y_offset=" + y_offset,       20, canvas.height - 15);

    now = new Date().toUTCString ();
    title="Projections in rust (" + now + ")";
    titleWidth=ctx.measureText(title + " ").width
    titleHeight=ctx.measureText("M").width
    ctx.fillText(title, canvas.width - titleWidth, canvas.height - titleHeight);

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

      
