#This is merely a draft
#This is broken, because we should never publish the entire modular forms!
def create_public_data_broken(M, num_senders, r):
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

def create_public_data(M, num_senders, r):
    st_bd = M.sturm_bound()
    if r < st_bd:
        raise ValueError("%s is less than the Sturm bound %s"%(r,st_bd))
    B = M.basis()
    R = M.base_ring()
    basis_coeffs = []
    coeff = []
    if B[0].coefficients(r+1)[r] == R.zero():
        raise ValueError("Invalid choice of r: %s"%r)
    for ix in range(0,num_senders):
        while True:
            coeff = [R.random_element() for ix in range(1,len(B))]
            fe_coeff_r = sum([coeff[ix-1]*B[ix].coefficients(r+1)[r] for ix in range(1,len(B))])
            if fe_coeff_r != R.zero():
                break
        basis_coeffs.append(coeff)
    return (basis_coeffs, B, r)

def send_to_auditor_broken(vote, sender, auditor, pub): 
    if not (vote == 1 or vote == 2):
        raise ValueError("Vote %s not supported yet."%vote)
    basis_coeffs = pub[0][sender]
    basis = pub[1]
    r = pub[2]
    M = basis[0].parent()
    R = M.base_ring()
    if sender < 0 or sender > len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)
    if auditor < 0 or auditor > M.sturm_bound():
        raise ValueError("Auditor id %s is invalid."%auditor)
    coeff_r = sum([basis_coeffs[ix]*basis[ix].coefficients(r+1)[r] for ix in range(0,len(basis))])
    coeff_audit = sum([basis_coeffs[ix]*basis[ix].coefficients(auditor+1)[auditor] for ix in range(0,len(basis))])
    vote_ff = R(vote)
    return vote_ff*coeff_r^(-1)*coeff_audit

def send_to_auditor(vote, sender, auditor, pub): 
    if not (vote == 1 or vote == 2):
        raise ValueError("Vote %s not supported yet."%vote)
    basis_coeffs = pub[0][sender]
    basis = pub[1]
    r = pub[2]
    M = basis[0].parent()
    R = M.base_ring()
    if sender < 0 or sender > len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)
    if auditor < 0 or auditor > M.sturm_bound():
        raise ValueError("Auditor id %s is invalid."%auditor)
    coeff_r = sum([basis_coeffs[ix-1]*basis[ix].coefficients(r+1)[r] for ix in range(1,len(basis))])
    coeff_audit = sum([basis_coeffs[ix-1]*basis[ix].coefficients(auditor+1)[auditor] for ix in range(1,len(basis))])
    vote_coeff_r = basis[0].coefficients(r+1)[r]
    vote_ff = R(vote)
    return (-coeff_r+vote_ff)*vote_coeff_r^(-1)*basis[0].coefficients[auditor+1][auditor]+coeff_audit

def reconstruct(sender, from_auditors, pub):
    basis_coeffs = pub[0]
    basis = pub[1]
    r = pub[2]
    M = basis[0].parent()
    if sender < 0 or sender >= len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)
    if len(from_auditors) != M.sturm_bound():
        raise ValueError("Incorrect number of auditors.")
    #pointless, but serves as a reminder that fsigma is constructed from
    #what the auditors sent
    fsigma = [from_auditor for from_auditor in from_auditors]
    #now we perform Gaußian elimination
    return 0
