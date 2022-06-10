from numpy import zeros,indices,einsum,array
from numpy.linalg import inv,det
from itertools import product

from scipy.fft import fftn, ifftn
from scipy.constants import e, pi, epsilon_0, m_e, h, eV
from scipy.sparse.linalg import eigsh, LinearOperator as Op
from scipy.signal import fftconvolve

def bulk(latvec,Z,N1=8,N2=4):
    """
    latvec=[v1,v2,v3] are the three lattice vectors.
    Z is (effective) charge.
    N1 is the number of sine functions used to approximate wave function.
    N2 is the degree in convolution
    """
    # reciprocal lattice vector with out 2pi
    rvec = inv(latvec).T
    coef0 = h**2 / 2 / m_e
    coef1 = Z * e**2 / (4 * pi**2 * epsilon_0 * det(latvec) )
    n = indices((2*N1+1,)*3) - N1

    pad = (2*N1+2*N2+1, )*3
    ker = zeros((2*N2+1,)*3)
    for i,j,k in product(range(-N2,N2+1),repeat=3):
        if i==0 and j==0 and k==0:
            continue
        ker[i+N2,j+N2,k+N2] = 1./sum(((i,j,k) @ rvec)**2)
    ker = fftn(ker,pad)

    def ans(n0,Neig):
        n0 = array(n0).reshape((3,1,1,1))
        def H(phi):
            phi = phi.reshape((2*N1+1,) * 3)
            k = einsum("ijkl,im->jklm", n+n0, rvec)
            H0 = coef0 * einsum("ijkl,ijkl,ijk->ijk", k,k,phi)
            H1 = -coef1 * ifftn(fftn(phi,pad)*ker)[N2:-N2,N2:-N2,N2:-N2]
            return (H0+H1.real).flatten()
        H = Op(matvec=H,shape=((2*N1+1)**3,)*2,dtype=float)
        E = eigsh(H,Neig,which="SA",return_eigenvectors=False) / eV
        return E
    return ans

# cell size of silicon
a = 5.43095e-10 / 4
latvec = [[a,a,a],[-a,-a,a],[a,-a,-a]]
silicon = bulk(latvec, Z=4,N1=12,N2=5)
print(silicon((.5,.3,.2), 20))

