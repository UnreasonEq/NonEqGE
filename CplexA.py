# Title: CplexA (Python implementation)
# Description: CplexA is a software package (available for Mathematica, Python, and Matlab) to compute probabilities and average properties of macromolecular assembly and its effects in gene regulation
# Authors:
#    Jose M. G. Vilar, IKERBASQUE & Biophysics Unit (CSIC-UPV/EHU), University of the Basque Country
#    Leonor Saiz, Department of Biomedical Engineering, University of California-Davis
# Copyright 2012, Jose M. G. Vilar and Leonor Saiz


from sympy import *

def F(fu,DG,S,te):
    retexp=fu*exp(-DG/te)
    for i in S:
    	retexp=retexp.subs(i,0)+retexp.subs(i,1)
    return expand(retexp)

def AveConf(fu,DG,S,te="RT"):
    if isinstance(te,str):
        te=sympify(te)
        #te=Symbol(tes,positive=True)
    if isinstance(fu,str):
        fu=sympify(fu)    
    if isinstance(DG,str):
        DG=sympify(DG)    
        #DG=parse_expr(DG,{tes:te})    
    if isinstance(S,str):
        S=symbols(S+" ")
        if isinstance(S,tuple):
            S=list(S)
        else:
            S=[S]    
    #retexp=F(fu,DG,S,te).evalf()/F(1,DG,S,te).evalf()
    retexp=F(fu,DG,S,te)/F(1,DG,S,te)
    return limit(retexp, Symbol("inf"), te*oo)
    
def eval_nd(fu):
    num,dem=fu.as_numer_denom()
    #return num.evalf()/dem.evalf()
    return num/dem
    
def subs(expr, replacee, replacer):
    if isinstance(replacee,str):
        replacee=sympify(replacee)    
    if isinstance(replacer,str):
        replacer=sympify(replacer)
    return expr.subs(replacee,replacer)
	
def pythonFunction(exp, vars):
    return eval("lambda "+str(vars)+': '+str(exp))