Code name: InterNewton.c
Author: Marco Antonio Esquivel Basaldua
Date: 15/10/19

Code description:
    This code captures n x position values in "Enter points.txt" corresponding to the interpolation points to work with.
    According to user selection, points in "Enter points.txt" are interpolated in a given funtion by Newton method.
    The original function, the interpolated function and the interpolation points are ploted on the same graph,
    data to plot are stored in exit .txt files.
    
Inputs:
- "Enter points.txt" containing the points to use as interpolation points, first line refers to the number of points to use in the file,
    USER CAN EDIT THIS FILE, make sure these values are in the interval for the chosen function.
- "graph.txt" instructions to plot graphs and points. DO NOT EDIT THIS FILE.
    
Outputs:
- "Exit points.txt" containing the ponits to graph the chosen function and the interpolating polinomyal the given display window.
    First column refers to x coordinate, second column refers to y coordinate of original chosen function, third column refers to y coordinate of interpolating polinomyal.
- "Interpolation points.txt" x and y coordinates of interpolating points read from "Enter points.txt"
- "InterNewton Interpolation.png" this image file contents the grpah of the chosen function, the graph of interpolating polinomyal and the interpolation points

Compiler instructions:
- By positioning in the folder where "InterNewton.c" and "Enter points.txt" are stored, use the Makefile placed in this folder, just write "make" and press enter in the terminal

- Results:
    Interpolation polinomyal using InterNewton method given points in "Enter points.txt"
        
