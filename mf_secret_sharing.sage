#This is merely a draft

#in this version we only publish part of the modular forms, the remaining
#coefficient (for the modular forms basis) is used in the voting
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

def send_to_auditor(vote, sender, auditor, pub): 
    basis_coeffs = pub[0][sender]
    basis = pub[1]
    r = pub[2]
    M = basis[0].parent()
    R = M.base_ring()
    if sender < 0 or sender > len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)
    if auditor < 0 or auditor > M.sturm_bound():
        raise ValueError("Auditor id %s is invalid."%auditor)
    coeff_r = sum([basis_coeffs[ix-1]*basis[ix].coefficients([r])[0] for ix in range(1,len(basis))])
    coeff_audit = sum([basis_coeffs[ix-1]*basis[ix].coefficients([auditor])[0] for ix in range(1,len(basis))])
    vote_coeff_r = basis[0].coefficients([r])[0]
    vote_ff = R(vote)
    #coefficient at prec `auditor` of
    #(-c(c1*f1+c2*f2;r) + vote) * c(f0;r)^(-1) * f0 + c1*f1+c2*f2
    #note that it at precision r is exactly `vote`
    #if coeff_r is equal to vote_ff, we might have a poopy situation
    return (-coeff_r+vote_ff)*vote_coeff_r^(-1)*basis[0].coefficients([auditor])[0]+coeff_audit

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
    #let's do it by hand first
    return 0

R = GF(7)
M = ModularForms(29, 2, base_ring=R)
pub = ([[R(1),R(3)], [R(4), R(6)], [R(3), R(3)]], M.basis(), 100)
votes = [0, 1, 1]
auditors = [R.zero() for ix in range(0,M.sturm_bound()+1)]
for sender in range(0,3):
    for auditor in range(0,len(auditors)):
        auditors[auditor] += send_to_auditor(votes[sender], sender, auditor, pub)

sender_recv=[[],[],[]]
for auditor in range(0,len(auditors)):
    for sender in range(0,3):
        #need id because order is important, and comms might be async
        sender_recv[sender].append((auditor,auditors[auditor]))
