#This is merely a draft

def create_public_data(M, num_senders, r):
    st_bd = M.sturm_bound()
    if r < st_bd:
        raise ValueError("%s is less than the Sturm bound %s"%(r,st_bd))
    B = M.basis()
    R = M.base_ring()
    basis_coeffs = []
    coeff = []
    for ix in range(0,num_senders):
        while True:
            coeff = [R.random_element() for b in B]
            fe_coeff_r = sum([coeff[ix]*B[ix].coefficients(r+1)[r] for ix in range(0,len(B))])
            if fe_coeff_r != R.zero():
                break
        basis_coeffs.append(coeff)
    return (basis_coeffs, B, r)

def reconstruct(sender, pub):
    basis_coeffs = pub[0]
    basis = pub[1]
    r = pub[2]
    if sender < 0 or sender >= len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)

    return 0
