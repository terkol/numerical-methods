import numpy as np

# More stable method for computing e**x
def myexp(x):
    k = round(x/np.log(2))  
    r = x-k*np.log(2)
    return (2**k)*np.exp(r)
 
if __name__ == "__main__":
    print(myexp(20))  
    print(np.exp(20))  
