# Numerical Methods
This repository contains source codes for the functionalities which are used in the course CH-2-16(MO) i.e., Numerical Methods of our curriculum.

<h2>Solutions for Non-linear equations</h2>
<h3>Contains:</h3>
    <b>1.</b> Bisection Method
<br><b>2.</b> Secant Method
<br><b>3.</b> Newton's Method
<br><b>4.</b> Fixed Point Iteration Method
<br>

<h2>Solutions for set of linear Equations</h2>
<h3>Contains:</h3>
<b>1.</b> Direct Methods<br>
<b>a.</b> Gauss Elimination -> considers the diagonal elements as not equal to zero. Can be solved if we use pivot finding.<br>
<b>b.</b> Gauss-Jordan Elimination -> yet to be done<br>
<b>c.</b> Solving Using Crout's Decomposition<br>
<br>
<b>2.</b> Iterative Methods<br>
<b>a.</b> Jacobi Iteration -> not completed<br>
<b>b.</b> Gauss-Siedel Iteration -> takes input of number of rows and columns multiple times due to my coding style
<h2>Solutions for set of non-linear equations by Newton-Raphson procedure</h2>Function f is the system of equations -> all combined into a single function. Hence, the use of int (i).
Each index (i) corresponds to a different non-linear equation.
<br></br>
<h2>Interpolation</h2>
Lagrange Interpolation is done!
Spline interpolation's implementation not done!
<br></br>
More might be added soon!
<br>
<br>
<br>
<h1>IMPORTANT NOTE</h1>
<br>
<a href = "https://www.google.com/search?q=ternary+operators+in+cpp&sca_esv=579920261&sxsrf=AM9HkKnf-FRDuZ0gGItDxA63hqDSNTaETw%3A1699305222620&ei=BldJZZTEJZKa4-EPlfK92As&ved=0ahUKEwiU36KhpbCCAxUSzTgGHRV5D7sQ4dUDCBA&uact=5&oq=ternary+operators+in+cpp&gs_lp=Egxnd3Mtd2l6LXNlcnAiGHRlcm5hcnkgb3BlcmF0b3JzIGluIGNwcDIFEAAYgAQyBhAAGBYYHjIGEAAYFhgeMgYQABgWGB4yCBAAGIoFGIYDSMAlUNsGWPYgcAJ4AZABAJgB-AGgAbMNqgEFMC4xLje4AQPIAQD4AQHCAgoQABhHGNYEGLADwgIKEAAYigUYsAMYQ8ICBxAAGIoFGEPCAgoQABiABBgUGIcCwgIHEAAYDRiABMICCRAAGA0YgAQYCsICChAAGBYYHhgPGAriAwQYACBBiAYBkAYK&sclient=gws-wiz-serp">
DO NOT USE TERNARY OPERATORS IN THIS REPOSITORY. YOU WANNA KNOW WHY? ASK THAT GODDAMNED OPERATOR WHY IT DOESN'T WORK.</a>
<a> Design decision: i did't know to link header and its implementation(cpp) files, so i just worte the implementation in the header files iteslf. Pls do not kill me! 
<!-- # Workflow
if ( function.works() && function.belongs_to(matrix) )<br>
{
<br>(matrix.h) <-- function<br>
}<br>
else<br>
{
<br>(others.h) <-- function<br>
}<br> -->