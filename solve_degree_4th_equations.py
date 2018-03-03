#!/usr/bin/python
# -*- coding: utf-8 -*-



import numpy as np
import cmath as mth
import math
import scipy
from cmath import *
from collections import Counter


def find_roots(polynome):
        """ #this function solves real coefficient polynomial equations of degree smaller than 4;
        #for degrees strictly higher than 5, there are no closed form, so we do not venture into these territories
        
        #------------------
        
        #example of call
        find_roots(polynome)
        
        #------------------
        
        #syntax of input
        polynome=[a,b,c,d,e] 
        # this is for polynomial
        #ax^4+bx^3+cx^2+dx+e
        
        #------------------
        
        #outputs are of type  
        roots,'txt giving the type of roots / how to interpret the roots'
        
        #------------------
        
        #examples of inputs
        #polynome=[1,-4,6,-4,1] # (x-1)^4
        #polynome=[1,0,-2,0,1] # (x^2-1)^2
        #polynome=[1,0,-3,0,2] # (x^2-1)(x^2-2) (+/-1 and +/- sqrt(2) i.e. almost +/- 1.41
        #polynome=[1,0,0,0,1] # roots are +/-exp(i*pi/4) and +/-exp(-i*pi/4)

        #polynome=[0,1,8,46,68] # roots are -2, -(3+5i) and the conjuagte
        #polynome=[1,-3,-2,12,-8] # roots are 1,2,+2,-2 ---> NOT DONE AS YET / test passed 
        #polynome=[1,-12,86,-268,533] # roots are 2+3i and 4+5i conjugates
        #polynome=[1,-3,-9,63,-52] # roots: -2i-3 and conjugate and also +1 and -4 
        #polynome=[1,2,-13,-14,24] # roots:1,-2,3,-4 

        #polynome=[1,-8,21.5,-22,0.75*8.75] # 4 real roots:1+/-0.5 and 3+/-0.5
        #polynome=[1,- 25.15 ,195.665 , - 474.791 , + 355.224] # expand (x-1.6)(x-1.9)(x-10.25)(x-11.4) --> mathworld! wolfram rocks ?
        
        
        #polynome=[0,0,0,0,0] # nul polynomial
        #polynome=[0.0,0.0,0.0,0.0,0.0] # nul polynomial
        #polynome=[0.0,0.0,1,0,1] # x^2+1
        #polynome=[0.0,0.0,2,2*(np.sqrt(2)-1),-2*np.sqrt(2)] # 2(x+sqrt(2)) (x-1)
        #polynome=[0.0,0.0,1,2,2] # (x+(1+i)) (x+(1-i))
        #polynome=[0.0,1,-1,1,-1] # (x-1)(x^2+1)
        #polynome=[0.0,1,np.sqrt(3)+1,np.sqrt(3)-2,-2*np.sqrt(3)] # -1,2,sqrt(3)
        #polynome=[0.0,1,-5,3,9] # 1 and -3 double
        #polynome=[0.0,2,-10,6,18] # 1 and -3 double
        #polynome=[1,4,6,4,1] # -1 multiplicity 4
        #polynome=[0.0,1,0,-3,-2] # (x-2)(x+1)^2
        #polynome=[0.0,1,-3,0,4] # (x-2)^2.(x+1)
        #polynome=[0.0,1,5,9,5] # (x+1)(x+2-i)(x+2+i)

        """
        [a,b,c,d,e]=polynome # ax^4+bx^3+cx^2+dx+e
        default_answer=[[],'default_answer']
        if a==0:
                if b==0:
                        if c==0:
                                if d==0:
                                        if e==0:
                                                return [],'zero polynomials -- any complex number is a root'
                                        else:
                                                return [],'constant non zero polynomials -- no complex number is a root'
                                else:
                                        return [-1.0*e/d],' one root for degree1 polynomials'

                        else:
                                delta=d**2-4*c*e
                                if delta>0:
                                        return [(-1.0*d+np.sqrt(delta))/(2.*c),(-1.0*d-np.sqrt(delta))/(2.*c)],' two real roots for degree2 polynomial'
                                elif delta==0:
                                        return [(-1.0*d)/(2.*c)],' one double real root for degree2 polynomial'
                                elif delta<0:
                                        return [(-1.0*d+complex(0,1)*np.sqrt(-1.0*delta))/(2*c)],' two complex conjugate roots for degree2 polynomial - we give one'
                                else:
                                        return [],' too strange degree two polynomial'
                else:
                        #degree 3 polynomials
                        # we compute the depressed cubic: x^3+px+q=0
                        #but first we divide by b, it is easy and painless
                        e=(1.0*e)/(1.*b)
                        d=(1.0*d)/(1.*b)
                        c=(1.0*c)/(1.*b)
                        b=1.0
                        
                        #we have X^3+cX^2+dX+e=0 to solve
                        #we suppose X=x+P
                        P=-c/3.0
                        #we solve for x and then add x to P to find X -- thus solving for X

                        p=3.*P**2+2.*P*c+d
                        q=P**3+P**2*c+d*P+e
                        
                        delta=-(4.*p**3+27.*q**2)
                        if delta>0:
                                
                                z=complex(-q/2.0,np.sqrt(delta/27.0)/2.0)
                                
                                theta=phase(z)
                                rho=abs(z)
                                
                                crude_roots=np.array([(2*exp((k*2*pi/3+theta/3.0)*complex(0,1))*np.power(rho,1/3)).real for k in [0,1,2]])
                                
                                # other possbility below
                                #crude_roots=np.array([(2*rect(1,(k*2*pi/3+theta/3.0))*np.power(rho,1/3)).real for k in [0,1,2]])

                                crude_roots=crude_roots+P
                                return crude_roots, 'three distinct real roots'
                        elif delta==0:
                                #this should be replaced by a isclose instance --
                                #though there are questions on whether the results would then be consistent for very small input coefficients... i leave that for future generations!
                                
                                z=-q/2.0
                                
                                theta=phase(z)
                                rho=abs(z)

                                crude_roots=np.array([(2*exp((k*2*pi/3+theta/3.0)*complex(0,1))*np.power(rho,1/3)).real for k in [0,1,2]])
                                crude_roots=crude_roots+P

                                #test below
                                #crude_roots=[1,2,5]
                                #crude_roots=[-1,-10,5]
                                #crude_roots=[50,-10,5]
                                
                                argmin_roots=np.argmin([abs(crude_roots[0]-crude_roots[1]),\
                                                        abs(crude_roots[0]-crude_roots[2]),\
                                                        abs(crude_roots[1]-crude_roots[2])])
                                multiplicity=[[[0,2],[2,1]],\
                                              [[0,1],[2,1]],\
                                              [[1,0],[2,1]],\
                                ]
                                #always 2,1 for the multiplicities! --- this is assumed in the 4th degree polynomial solver so do not mess with it.

                                this_case=multiplicity[argmin_roots]
                                crude_roots=[crude_roots[ind] for ind in this_case[0] ]

                                return crude_roots, 'two real roots - multiplicity-'+str(this_case[1])
                        else: #delta<0 --- should do some other condition but laziness is overwhelming # by the way life is too short for safely removing usb key 
                                z_plus=-q/2.0+np.sqrt(-delta/27.0)/2.0
                                z_minus=-q/2.0-np.sqrt(-delta/27.0)/2.0

                                theta_plus=phase(z_plus)
                                rho_plus=abs(z_plus)

                                theta_minus=phase(z_minus)
                                rho_minus=abs(z_minus)
                                
                                
                                cubic_root_plus=np.sign(z_plus)*np.power(rho_plus,1/3)
                                cubic_root_minus=np.sign(z_minus)*np.power(rho_minus,1/3)
                                

                                j_cub_root_unit=exp(2*pi/3*complex(0,1))
                                
                                # the complete thing would be the following -- but let us not go into too many details
                                # crude_roots=[P+cubic_root_plus+cubic_root_minus,\
                                #              P+j_cub_root_unit*cubic_root_plus+j_cub_root_unit**2*cubic_root_minus,\
                                #              P+j_cub_root_unit**2*cubic_root_plus+j_cub_root_unit*cubic_root_minus]
                                crude_roots=[P+cubic_root_plus+cubic_root_minus,\
                                             P+j_cub_root_unit*cubic_root_plus+j_cub_root_unit**2*cubic_root_minus]

                                
                                return crude_roots, 'one real root - and two complex conjugate roots -- we give only one complex root'
        
        else: # a!=0
                # we should solve a degree 4th equation ax^4+bx^3+cx^2+dx+e=0 --> let's get to work
                #so we start by x=X+A with A=-b/4a
                A=-b/(4.0*a)

                #we try to find the depressed form
                #alpha.X^4+beta.X^2+gamma.X+delta=0
                # and solve for X... and then add A to X to retrieve x
                
                alpha=a
                beta=6*A**2*a+3*A*b+c
                gamma=4*A**3*a+3*A**2*b+2*A*c+d
                delta=a*A**4+b*A**3+c*A**2+d*A+e

                # and now we strive to solve the resolvent equation: (bombelli's resolvent below)
                #resolvent_polynomial=[0,1,beta/(alpha*1.0),(beta**2/(4.0*alpha**2)-delta/(alpha*1.0)),((gamma/(alpha*1.0))**2)/8.0]
                # but hey this is only a joke --> we will use Euler's resolvent actually
                
                # euler resolvent below (more straight forward than bombelli's in my biased perspective) -> but possibly numeric methods
                #are not best suited to find whether roots are real or complex or have a precise value
                #all these are approximations
                
                
                resolvent_polynomial=[0,1,beta/(alpha*2.0),beta**2/(16.0*alpha**2)-delta/(alpha*4.0),-(gamma/(alpha*8.0))**2]
                
                
                resolvent_output=find_roots(resolvent_polynomial)
                if 0:
                        Zx=resolvent_output[0][0]
                        res_should_be_zero=Zx**3*resolvent_polynomial[1]+Zx**2*resolvent_polynomial[2]+Zx*resolvent_polynomial[3]+resolvent_polynomial[4]
                        print(res_should_be_zero)
                        print("res_should_be_zero--0")

                        Zx=resolvent_output[0][1]
                        res_should_be_zero=Zx**3*resolvent_polynomial[1]+Zx**2*resolvent_polynomial[2]+Zx*resolvent_polynomial[3]+resolvent_polynomial[4]
                        print(res_should_be_zero)
                        print("res_should_be_zero--1 \n ligne 265")
                THE_product=-1.0*gamma/(alpha*8.0)
                
                if gamma==0:
                        #print('gamma==0 ---> this is a bicubic equation , we should take great care because using euler"s root method would produce difficult to disentangle results')
                        #alpha.X^4+beta.X^2+delta=0
                        deg_2_polynomial=[0,0,alpha,beta,delta]
                        deg_2_polynomial_solutions=find_roots(deg_2_polynomial)
                        
                        # two cases:
                        #real roots ---> no pb
                        #complex roots --> beware
                        
                        if 'real' in deg_2_polynomial_solutions[1]:
                                if 'one' in deg_2_polynomial_solutions[1]:
                                        # there is only one double root --> therefore solutions are doubled
                                        roots=find_roots([0,0,1,0,-deg_2_polynomial_solutions[0][0]])

                                        crude_roots=roots[0]
                                        crude_roots=[A+el for el in crude_roots]
                                        
                                        txt=roots[1]+' -- doubled --> 4 roots'
                                        return crude_roots+crude_roots,txt
                                else:
                                        
                                        #there are two roots-- the equation is of the form (X^2-r1)(X^2-r2) with r1!=r2
                                        roots0=find_roots([0,0,1,0,-deg_2_polynomial_solutions[0][0]])
                                        roots1=find_roots([0,0,1,0,-deg_2_polynomial_solutions[0][1]])
                                        crude_roots0=roots0[0]
                                        crude_roots1=roots1[0]
                                        
                                        crude_roots=crude_roots0+crude_roots1
                                        crude_roots=[A+el for el in crude_roots]

                                        txt=roots0[1]+' AND '+roots1[1]+ ' : makes 4 roots'
                                        return crude_roots,txt
                        elif 'complex' in deg_2_polynomial_solutions[1]:
                                
                                # in this case only one of the two complex conjugate root is given

                                root_res=deg_2_polynomial_solutions[0][0]

                                theta_root_res=phase(root_res)
                                abs_root_res=abs(root_res)

                                complex_sqrt=rect(np.sqrt(abs_root_res),theta_root_res/2.0)
                                
                                return [complex_sqrt],'we give only one complex a+ib and the supplementary solutions are a-ib,-a+ib,-a-ib, which makes 4 roots'
                        
                elif 'two r' in resolvent_output[1]:
                        # we are not in the case of a bicubic -- the resolvent has two real roots one of these is a double root
                        epsilons=[[a,b,c] for a in [-1,1] for b in [-1,1] for c in [-1,1]]

                        three_sqrt_of_res_root=[]
                        for root_res in resolvent_output[0]:
                                sqrt_root_res=find_roots([0,0,1,0,-1.0*root_res])
                                three_sqrt_of_res_root.append(sqrt_root_res[0][0])
                        #print (three_sqrt_of_res_root) # multiplicity is 2,1, there are still 2 roots here actually :( -- contrary to what the name let people imagine
                        #print ('three_sqrt_of_res_root \n ')
                        
                        three_sqrt_of_res_root=[three_sqrt_of_res_root[0],three_sqrt_of_res_root[0],three_sqrt_of_res_root[1]]
                        this_product=-1.0*gamma/(alpha*8.0)

                        crude_roots=[]

                        #differences_between_what_it_should_be_and_what_it_is
                        DBWISBAWII=[]
                        
                        
                        
                        for ind_eps in range(len(epsilons)):
                                
                                DBWISBAWII.append(this_product-np.product(np.multiply(epsilons[ind_eps],three_sqrt_of_res_root)))
                        
                                
                                crude_roots.append(np.sum(np.multiply(epsilons[ind_eps],three_sqrt_of_res_root)))
                        
                        #let us take the 4 smallest differences

                        good_root_indices=np.abs(np.array(DBWISBAWII)).argsort()[:4]
                        crude_roots=[crude_roots[el]+A for el in good_root_indices]
                        DBWISBAWII=[DBWISBAWII[el] for el in good_root_indices]
                        
                        return crude_roots,'four roots, and the following nbs should be ##very small## '+str(DBWISBAWII)
                
                elif 'one r' in resolvent_output[1]:
                        
                        # there is only one real root to the resolvent
                        
                        # this is always real and then complex [real, complex] not the other way around!
                        
                        # below are computed the signs necessary for computing the roots of the original equation from the euler's resolvent's roots
                        epsilons=[[a,b,c] for a in [-1,1] for b in [-1,1] for c in [-1,1]]

                        three_sqrt_of_res_root=[]
                        if 0:
                                #a test that has prooved useful to me during debugging
                                print(resolvent_output[0])
                                print("resolvent_output[0] \n line 423 \n ")
                                print(resolvent_output[0][0]*abs(resolvent_output[0][1])**2)
                                print("resolvent_output[0][0]*abs(resolvent_output[0][1])**2 \n line 425 \n ")
                                print((gamma/(alpha*8.0))**2)
                                print("(gamma/alpha*8)**2 --4th --- 427")
                        
                                input()

                        #"resolvent_output[0]"   =    [r (real root), a'+ib' (one of the two complex root )]
                        #taking the square roots of these and adding with +/- signs in front of these, means  --> +/- sqrt(r) {+2a or 2ib or -2a or -2ib} with a+ib=sqrt(a'+ib')
                        #generally the solutions are
                        #+/- sqrt(r)+ {+2a or 2ib or -2a or -2ib}

                        real_or_complex=[]
                        
                        for root_res in resolvent_output[0]:
                                if isinstance(root_res,float):
                                        real_or_complex.append('r')
                                        sqrt_root_res=find_roots([0,0,1,0,-1.0*root_res])
                                        # the result of this can be any type of solutions of polynomial 2nd degree
                                        three_sqrt_of_res_root.append(sqrt_root_res[0][0])

                                elif isinstance (root_res,complex):
                                        real_or_complex.append('c')

                                        theta_root_res=phase(root_res)
                                        abs_root_res=abs(root_res)

                                        complex_sqrt=rect(np.sqrt(abs_root_res),theta_root_res/2.0)
                                        three_sqrt_of_res_root.append(complex_sqrt) 
                                else:
                                        print ('big pb non float and non complex!!! line 395')
                                        raise('big pb non float and non complex!!!line 395')
                                
                        
                        # multiplicity is 2,1, there are still 2 roots here actually :(

                        all_the_roots=[three_sqrt_of_res_root[0],-three_sqrt_of_res_root[0]] # these are the two +/- real root of resolvent
                        crude_roots=[[A+el+2*three_sqrt_of_res_root[1].real,\
                                      A+el-2*three_sqrt_of_res_root[1].real,\
                                      A+el+complex(0,1)*2*three_sqrt_of_res_root[1].imag,\
                                      A+el-complex(0,1)*2*three_sqrt_of_res_root[1].imag] \
                                     for el in all_the_roots]
                        # we take into account complex roots of the complex root of resolvent with + or - in front +/-z1 +/-z2
                        # we still have to compute products of roots
                        # to choose the correct solutions

                        product_of_resolvent_roots=[[el*abs(three_sqrt_of_res_root[1])**2,\
                                                     el*abs(three_sqrt_of_res_root[1])**2,\
                                                     -el*abs(three_sqrt_of_res_root[1])**2,\
                                                     -el*abs(three_sqrt_of_res_root[1])**2]
                                                    for el in all_the_roots]
                        this_product=-1.0*gamma/(alpha*8.0)
                        crude_roots=crude_roots[0]+crude_roots[1]
                        product_of_resolvent_roots=product_of_resolvent_roots[0]+product_of_resolvent_roots[1]




                        #differences_between_what_it_should_be_and_what_it_is
                        DBWISBAWII=[this_product-el for el in product_of_resolvent_roots]
                        

                        #let us take the 4 smallest differences

                        good_root_indices=np.abs(np.array(DBWISBAWII)).argsort()[:4]
                        crude_roots=[crude_roots[el] for el in good_root_indices] # A has already been added
                        DBWISBAWII=[DBWISBAWII[el] for el in good_root_indices]
                        
                        return crude_roots,'four roots, and the following nbs should be ##very small## '+str(DBWISBAWII)
                        


                
                elif 'three d' in resolvent_output[1]:
                        # there are three distinct real roots to the resolvent
                        epsilons=[[a,b,c] for a in [-1,1] for b in [-1,1] for c in [-1,1]]

                        three_sqrt_of_res_root=[]
                        for root_res in resolvent_output[0]:
                                sqrt_root_res=find_roots([0,0,1,0,-root_res])
                                three_sqrt_of_res_root.append(sqrt_root_res[0][0])
                                
                        crude_roots=[]

                        #let us take the 4 smallest differences
                        #differences_between_what_it_should_be_and_what_it_is
                        DBWISBAWII=[]
                        
                        
                        for ind_eps in range(len(epsilons)):
                                
                                DBWISBAWII.append(abs(THE_product-np.product(np.multiply(epsilons[ind_eps],three_sqrt_of_res_root))))
                                
                                crude_roots.append(A+np.sum(np.multiply(epsilons[ind_eps],three_sqrt_of_res_root)))
                        good_indices=np.array(DBWISBAWII).argsort()[:4]
                        good_roots=[crude_roots[index_bidule] for index_bidule in good_indices]
                        DBWISBAWII=[DBWISBAWII[el] for el in good_indices]


                        
                        return good_roots,'four roots, and the following nbs should be ##very small## '+str(DBWISBAWII)

                        
                else:
                        raise("strange case 4th degree equation -- case 1")
        
                        
                        
        raise("strange case 4th degree equation -the other case: case 2")

if __name__ == "__main__":
        #main()


        #polynme de degré4 
        #polynome=[1,-4,6,-4,1] # (x-1)^4
        #polynome=[1,0,-2,0,1] # (x^2-1)^2
        #polynome=[1,0,-3,0,2] # (x^2-1)(x^2-2) (+/-1 and +/- sqrt(2) i.e. almost +/- 1.41
        #polynome=[1,0,0,0,1] # roots are +/-exp(i*pi/4) and +/-exp(-i*pi/4)

        #polynome=[0,1,8,46,68] # roots are -2, -(3+5i) and the conjuagte
        #polynome=[1,-3,-2,12,-8] # roots are 1,2,+2,-2 ---> NOT DONE AS YET / test passed ---> WHOOOPIE!
        #polynome=[1,-12,86,-268,533] # roots are 2+3i and 4+5i conjugates
        #polynome=[1,-3,-9,63,-52] # roots: -2i-3 and conjugate and also +1 and -4 
        #polynome=[1,2,-13,-14,24] # roots:1,-2,3,-4
        #polynome=[1,-8,21.5,-22,0.75*8.75] # 4 real roots:1+/-0.5 and 3+/-0.5
        polynome=[1,- 25.15 ,195.665 , - 474.791 , + 355.224] # expand (x-1.6)(x-1.9)(x-10.25)(x-11.4) --> mathworld!
        
        #polynome=[0,0,0,0,0] # nul polynome
        #polynome=[0.0,0.0,0.0,0.0,0.0] # nul polynome
        #polynome=[0.0,0.0,1,0,1] # x^2+1
        #polynome=[0.0,0.0,2,2*(np.sqrt(2)-1),-2*np.sqrt(2)] # 2(x+sqrt(2)) (x-1)
        #polynome=[0.0,0.0,1,2,2] # (x+(1+i)) (x+(1-i))
        #polynome=[0.0,1,-1,1,-1] # (x-1)(x^2+1)
        #polynome=[0.0,1,np.sqrt(3)+1,np.sqrt(3)-2,-2*np.sqrt(3)] # -1,2,sqrt(3)
        #polynome=[0.0,1,-5,3,9] # 1 et -3 double
        #polynome=[0.0,2,-10,6,18] # 1 et -3 double
        #polynome=[1,4,6,4,1] # -1 multiplicity 4
        #polynome=[0.0,1,0,-3,-2] # (x-2)(x+1)^2
        #polynome=[0.0,1,-3,0,4] # (x-2)^2.(x+1)
        #polynome=[0.0,1,5,9,5] # (x+1)(x+2-i)(x+2+i)
        
        print(find_roots(polynome))# list de max 4 float
        print("find_roots(polynome)")# list de max 4 float
