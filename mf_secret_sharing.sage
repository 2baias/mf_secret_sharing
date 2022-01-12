#MIT License
#
#Copyright (c) 2022 Tobias Magnusson
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

#in this version we only publish part of the modular forms, the remaining
#coefficient (for the modular forms basis) is used in the voting
def create_public_data(M, num_senders, r):
    st_bd = M.sturm_bound()
    if r < st_bd:
        raise ValueError("%s is less than the Sturm bound %s"%(r,st_bd))
    B = M.basis()
    R = M.base_ring()
    if num_senders >= R.cardinality():
        raise ValueError("%s cannot contain all different vote sums."%R)
    basis_coeffs = []
    coeff = []
    if B[0].coefficients(r+1)[r] == R.zero():
        raise ValueError("Invalid choice of r: %s"%r)
    for ix in range(0,num_senders):
        while True:
            coeff = [R.random_element() for ix in range(1,len(B))]
            fe_coeff_r = sum([coeff[ix-1]*B[ix].coefficients([r])[0] for ix in range(1,len(B))])
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
    R = M.base_ring()
    if sender < 0 or sender >= len(basis_coeffs):
        raise ValueError("Sender id %s is invalid."%sender)
    if len(from_auditors) != (M.sturm_bound()+1):
        raise ValueError("Incorrect number of auditors.")
    #create fsigma, recall that it isn't necessarily received in order
    fsigma = [R.zero() for from_auditor in from_auditors]
    for aud_pair in from_auditors:
        fsigma[aud_pair[0]] = aud_pair[1]
    cols = [b.coefficients(range(0,M.sturm_bound()+1)) for b in basis]
    cols.append(fsigma)
    coeff_matrix = Matrix(cols).transpose()
    if coeff_matrix.rank() != M.dimension():
        raise RuntimeError("Could not reconstruct modular form.")
    coeff_matrix = coeff_matrix.echelon_form()
    coeff_mfs = [coeff_matrix[ix,M.dimension()] for ix in range(0,len(basis))]
    #extract `r`th coefficient of every modular form
    basis_coeff_r = [b.coefficients([r])[0] for b in basis]
    #multiply with coefficients computed via row reduction
    return sum([coeff_mfs[ix]*basis_coeff_r[ix] for ix in range(0,len(basis))])

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

vote_sum = [R.zero() for ix in range(0,3)]
for sender in range(0,3):
    vote_sum[sender] = reconstruct(sender, sender_recv[sender], pub)
