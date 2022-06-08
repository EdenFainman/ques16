import math

import sympy as sp
from sympy.utilities.lambdify import lambdify

def newtonRaphson_calc(f,f_d,f_dShow,start_point ,end_point,epsilon,i):
    """
    calculate the roots with the newton raphson method
    :param f: normal function
    :param f_d: derivative function
    :param f_dShow: function as it looks
    :param start_point: first point
    :param end_point: second point
    :param epsilon: minimum error
    :param i: 1 for derivative 0 for normal function
    """
    difference = 1
    guess = start_point
    iteration = 1
    if i==0:
        while(difference>=epsilon):
            normalF=f(guess)
            der=f_d(guess)
            nextGuess=guess-(normalF/der)
            difference=abs(nextGuess-guess)
            print(f"iteration: {iteration} Xr: {guess} f(x) {normalF} f'(x): {der}")
            lastGuess=guess
            guess=nextGuess
            iteration+=1
        print(f"The root is {lastGuess}\n")

    else:
        f_d2 = sp.diff(f_dShow,x)  # Derivative
        f_d2 = lambdify(x, f_d2)  # create the polinom to function
        while (difference >= epsilon):
            normalF = f_d(guess)
            der = f_d2(guess)
            nextGuess = guess - (normalF / der)
            difference = abs(nextGuess - guess)
            print(f"iteration: {iteration} Xr: {guess} f(x) {normalF} f'(x): {der}")
            lastGuess = guess
            guess = nextGuess
            iteration += 1

        if abs(f(lastGuess))<=epsilon:
            print(f"The root is {lastGuess}\n")
        else:
            print("After testing it can be seen that this root is not suitable because f(x) doesn't equal 0\n")
def newtonRaphson(f,start_point,end_point,epsilon):
    """
    calculates derivative and send the the function to be checked every segment
    :param f: normal function
    :param start_point: first point
    :param end_point: second point
    :param epsilon: minimum error
    """

    f_d = f.diff(x)  # Derivative
    f_dShow=f_d
    f_d = lambdify(x, f_d)  # create the polinom to function
    f = lambdify(x, f)
    range_check_newtonRaphson(f,f_d,f_dShow,start_point,end_point,epsilon,0)
    range_check_newtonRaphson(f,f_d,f_dShow,start_point,end_point,epsilon,1)

def range_check_newtonRaphson(f,f_d,f_dShow,start_point,end_point,epsilon,i):
    """
    do the calculation of the method every segment that there is a root
    :param f: normal function
    :param f_d: derivative function
    :param f_dShow: function as it looks
    :param start_point: first point
    :param end_point: second point
    :param epsilon: minimum error
    :param i: 1 for derivative 0 for normal function
    """
    a=start_point
    if i==0:
        while (a<end_point):
            if (f(a)*f(a+0.1)<0 ):
                newtonRaphson_calc(f,f_d,f_dShow,a,a+0.1,epsilon,0)
                a=a+0.1
            else:
                a=a+0.1
        return
    else:
        while (a<end_point):
            if (f_d(a)*f_d(a+0.1)<0 ):
                newtonRaphson_calc(f,f_d,f_dShow,a,a+0.1,epsilon,1)
                a=a+0.1
            else:
                a=a+0.1
    return



def bisection_calc(f,f_d,start_point ,end_point,epsilon,i,n):
    """
        calculate the roots with the bisection method
        :param f: normal function
        :param f_d: derivative function
        :param start_point: first point
        :param end_point: second point
        :param epsilon: minimum error
        :param i: 1 for derivative 0 for normal function
        :param n: the max iteration
        """
    a_n = start_point
    b_n = end_point
    k=1
    if i==0:
        while (abs(a_n-b_n)>epsilon):
            m_n = (a_n + b_n)/2  # the middle
            f_m_n = f(m_n)
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0:
                print("Found exact solution.")
                return m_n
            elif k >= n:
                print(k, " Is the max iteration according to error calculation\n")
                print("Bisection method fails,please select another methode ")
                return None
            print("Iteration:",k,",a=",a_n,",b=",b_n,",c=",m_n,",f(a)=",f(a_n),",f(b)=",f(b_n))
            k=k+1
        print("The root is : ",(a_n + b_n)/2,"\n")
    else:
        while (abs(a_n-b_n)>epsilon):
            m_n = (a_n + b_n)/2  # the middle
            f_m_n = f_d(m_n)
            if f_d(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f_d(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0:
                print("Found exact solution.")
                return m_n
            elif k>=n:
                print(k," Is the max iteration according to error calculation\n")
                print("Bisection method fails,please select another methode ")
                return None
            print("Iteration:",k,",a=",a_n,",b=",b_n,",c=",m_n,",f(a)=",f(a_n),",f(b)=",f(b_n))
            k=k+1

        if abs(f(a_n+b_n))<=epsilon:
            print("The root is : ",(a_n + b_n)/2,"\n")
        else:
            print("After testing it can be seen that this root is not suitable because f(x) doesn't equal 0\n")
            return

def secannt(f,start_point,end_point,epsilon):
    """
       calculates derivative and send the the function to be checked every segment

       :param f: normal function
       :param start_point: first point
       :param end_point: second point
       :param epsilon: minimum error
       """
    f_d = sp.diff(f,x) #Derivative
    f_d = lambdify(x, f_d) #create the polinom to function
    f = lambdify(x, f)
    range_check_secant(f,f_d,start_point,end_point,epsilon,1)
    range_check_secant(f,f_d,start_point,end_point,epsilon,0)

def range_check_secant(f,f_d,start_point,end_point,epsilon,i):
    """
    do the calculation of the method every segment that there is a root

    :param f: normal function
    :param f_d: derivative function
    :param f_dShow: function as it looks
    :param start_point: first point
    :param end_point: second point
    :param epsilon: minimum error
    :param i: 1 for derivative 0 for normal function
    """
    a=start_point
    if i==0:
        while (a<end_point):
            if (f(a)*f(a+0.1)<0  ):
                secant_calc(f,f_d,a,a+0.1,epsilon,0)
                a=a+0.1
            else:
                a=a+0.1
        return
    else:
        while (a<end_point):
            if (f_d(a)*f_d(a+0.1)<0 ):
                secant_calc(f,f_d,a,a+0.1,epsilon,1)
                a=a+0.1
            else:
                a=a+0.1
    return

def secant_calc(f,f_d,start_point ,end_point,epsilon,i):
    """
        calculate the roots with the secant method
        :param f: normal function
        :param f_d: derivative function
        :param f_dShow: function as it looks
        :param start_point: first point
        :param end_point: second point
        :param epsilon: minimum error
        :param i: 1 for derivative 0 for normal function
        """
    difference = 1
    firstGuess = start_point
    secondGuess=end_point
    iteration = 1
    if i==0:
        while(difference>=epsilon):
            nextGuess=(firstGuess*f(secondGuess)-secondGuess*f(firstGuess))/(f(secondGuess)-f(firstGuess))
            difference=abs(nextGuess-secondGuess)
            print(f"iteration: {iteration} Xr: {firstGuess} Xr+1: {secondGuess} ")
            lastGuess=secondGuess
            secondGuess=nextGuess
            iteration+=1
        print(f"The root is {lastGuess}\n")
    else:
        while(difference>=epsilon):
            nextGuess=(firstGuess*f_d(secondGuess)-secondGuess*f_d(firstGuess))/(f_d(secondGuess)-f_d(firstGuess))
            difference=abs(nextGuess-secondGuess)
            print(f"iteration: {iteration} Xr: {firstGuess} Xr+1: {secondGuess} ")
            lastGuess=secondGuess
            secondGuess=nextGuess
            iteration+=1

        if abs(f(lastGuess))<=epsilon:
            print(f"The root is {lastGuess}\n")
        else:
            print("After testing it can be seen that this root is not suitable because f(x) doesn't equal 0\n")


def calch(b,a,n):
    h= ((b-a))/n
    return h

def xarray(a,h,n):
    xarr=[]
    xarr.append(a)
    for i in range(a+1,n+1):
        xarr.append(xarr[i-1]+h)
    return xarr

def yarray(axarr,n,func):
    yarr=[]
    fd=lambdify(x,func)
    for i in range(0,n+1):
        yarr.append(fd(float(axarr[i])))
    return yarr


def simpsonMethod(yarr,h):

    num=yarr[0]
    for i in range(1,len(yarr)):
        if(i%2==0):
            num=num+2*yarr[i]
        else:
            num=num+4*yarr[i]
    num=num+yarr[len(yarr)-1]
    num=(1/3)*h*num
    return num

def range_check_bisection(f,f_d,start_point,end_point,epsilon,i,n):
    """
      calculates derivative and send the the function to be checked every segment
      :param f: normal function
      :param f_d: derivative function
      :param start_point: first point
      :param end_point: second point
      :param epsilon: minimum error
      :param i: if the function is derivative
      :param n: the max iteration
      """
    a=start_point
    if i==0:
        while (a<end_point):
            if ((f(a)>=0 and f(a+0.1)<=0 )or((f(a)<=0 and f(a+0.1)>=0 ))):
                bisection_calc(f,f_d,a,a+0.1,epsilon,0,n)
                a=a+0.1
            else:
                a=a+0.1
        return
    else:
        while (a<end_point):
            if ((f_d(a)>=0 and f_d(a+0.1)<=0 )or((f_d(a)<=0 and f_d(a+0.1)>=0 ))):
                bisection_calc(f,f_d,a,a+0.1,epsilon,1,n)
                a=a+0.1
            else:
                a=a+0.1
    return


def bisection(f,start_point,end_point,epsilon):
    """
          calculates derivative and send the the function to be checked every segment

          :param f: normal function
          :param start_point: first point
          :param end_point: second point
          :param epsilon: minimum error
          """
    a=end_point-start_point
    n = -((math.log(10**-10/a))/math.log(2))
    n=int(n+1)
    f_d = f.diff(x) #Derivative
    f_d = lambdify(x, f_d) #create the polinom to function
    f = lambdify(x, f)
    range_check_bisection(f,f_d,start_point,end_point,epsilon,1,n)
    range_check_bisection(f,f_d,start_point,end_point,epsilon,0,n)



epsilon=0.0001
x=sp.symbols('x')
f=((x**2*(pow(math.e,((-1)*(x**2)+5*x-3))))*(3*x-5))
start_point = 0
end_point = 3
print(f)

print("Bisection Method:")
bisection(f,start_point,end_point,epsilon)
print("Secannt Method")
secannt(f,start_point,end_point,epsilon)

arrx=(xarray(0,calch(0.5,0,1),10))

print("Simpson Method:")
print(str(simpsonMethod(yarray(arrx,10,f),calch(0.5,0.1,10))))







